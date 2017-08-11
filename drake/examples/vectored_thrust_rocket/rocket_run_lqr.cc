#include <memory>

#include <gflags/gflags.h>

#include "drake/common/find_resource.h"
#include "drake/lcm/drake_lcm.h"
#include "drake/multibody/joints/floating_base_types.h"
#include "drake/multibody/parsers/urdf_parser.h"
#include "drake/multibody/rigid_body_plant/drake_visualizer.h"
#include "drake/multibody/rigid_body_plant/rigid_body_plant.h"
#include "drake/multibody/rigid_body_tree.h"
#include "drake/systems/analysis/simulator.h"
#include "drake/systems/controllers/linear_quadratic_regulator.h"
#include "drake/systems/framework/diagram.h"
#include "drake/systems/framework/diagram_builder.h"

namespace drake {
namespace examples {
namespace rocket {
namespace {

// Simple example which simulates the rocket,
// started near the upright, with an
// LQR controller designed to stabilize the unstable fixed point
// Run drake-visualizer to see the animated result.

DEFINE_double(realtime_factor, 1.0,
              "Playback speed.  See documentation for "
              "Simulator::set_target_realtime_rate() for details.");

std::unique_ptr<systems::AffineSystem<double>> GenerateBalancingLQRController(
    const systems::RigidBodyPlant<double>& rocket) {
  auto context = rocket.CreateDefaultContext();

  // Set nominal torques to zero, and force to
  // the total mass of the rocket.
  Eigen::Matrix<double, 3, 1> u_nominal;
  u_nominal << rocket.get_rigid_body_tree().getMass() * 9.81, 0.0, 0.0;
  context->FixInputPort(0, u_nominal);

  // Set nominal state to the upright fixed point.
  auto x_nominal = context->get_mutable_continuous_state_vector();
  DRAKE_ASSERT(x_nominal != nullptr);
  x_nominal->SetZero();
  x_nominal->SetAtIndex(3, 1.0);  // w_w for quaternion.

  // Setup LQR Cost matrices.
  DRAKE_ASSERT(x_nominal->size() == 17);
  Eigen::Matrix<double, 17, 17> Q = Eigen::Matrix<double, 17, 17>::Identity();
  Q(0, 0) = 10;   // x
  Q(1, 1) = 10;   // y
  Q(2, 2) = 10;   // z
  Q(3, 3) = 100;  // w_w
  Q(4, 4) = 100;  // w_x
  Q(5, 5) = 100;  // w_y
  Q(6, 6) = 100;  // w_z
  Q(7, 7) = 5;    // gimbal_x
  Q(8, 8) = 5;    // gimbal_y

  Eigen::Vector2d R = Eigen::Vector2d::Constant(1.0);

  return systems::controllers::LinearQuadraticRegulator(rocket, *context, Q, R);
}

int do_main(int argc, char* argv[]) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  lcm::DrakeLcm lcm;
  systems::RigidBodyPlant<double>* plant = nullptr;
  systems::DiagramBuilder<double> builder;

  {
    auto tree = std::make_unique<RigidBodyTree<double>>();
    parsers::urdf::AddModelInstanceFromUrdfFileToWorld(
        FindResourceOrThrow(
            "drake/examples/vectored_thrust_rocket/Rocket.urdf"),
        multibody::joints::kQuaternion, tree.get());
    plant = builder.AddSystem(
        std::make_unique<systems::RigidBodyPlant<double>>(std::move(tree)));
  }

  const RigidBodyTree<double>& tree = plant->get_rigid_body_tree();
  auto publisher = builder.AddSystem<systems::DrakeVisualizer>(tree, &lcm);
  builder.Connect(plant->get_output_port(0), publisher->get_input_port(0));

  auto controller = builder.AddSystem(GenerateBalancingLQRController(*plant));
  controller->set_name("controller");
  builder.Connect(plant->get_output_port(0), controller->get_input_port());
  builder.Connect(controller->get_output_port(), plant->get_input_port(0));

  auto diagram = builder.Build();
  systems::Simulator<double> simulator(*diagram);
  systems::Context<double>& rocket_context =
      diagram->GetMutableSubsystemContext(*plant,
                                          simulator.get_mutable_context());

  // Force initial state to be not at the fixed point
  auto x0 = rocket_context.get_mutable_continuous_state_vector();
  DRAKE_DEMAND(x0 != nullptr);
  x0->SetAtIndex(0, 1.0);
  x0->SetAtIndex(3, 0.8);
  x0->SetAtIndex(4, 0.1);
  x0->SetAtIndex(5, 0.1);
  x0->SetAtIndex(6, 0.1);

  simulator.set_target_realtime_rate(FLAGS_realtime_factor);
  simulator.Initialize();
  simulator.StepTo(10);
  return 0;
}

}  // namespace
}  // namespace rocket
}  // namespace examples
}  // namespace drake

int main(int argc, char* argv[]) {
  return drake::examples::rocket::do_main(argc, argv);
}
