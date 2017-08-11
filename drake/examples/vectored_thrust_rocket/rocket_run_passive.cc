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
#include "drake/systems/framework/diagram.h"
#include "drake/systems/framework/diagram_builder.h"

namespace drake {
namespace examples {
namespace rocket {
namespace {

// Simple example which simulates the (passive) rocket.  Run drake-visualizer
// to see the animated result.

DEFINE_double(realtime_factor, 1.0,
              "Playback speed.  See documentation for "
              "Simulator::set_target_realtime_rate() for details.");

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
  auto diagram = builder.Build();

  systems::Simulator<double> simulator(*diagram);

  systems::Context<double>& rocket_context =
      diagram->GetMutableSubsystemContext(*plant,
                                          simulator.get_mutable_context());

  printf("Num inputs: %d\n", (int) plant->get_input_port(0).size());

  Eigen::Matrix<double, 3, 1> u0;
  u0 << 0.0, 0.0, 0.0; // Order: force, gimbal_x, gimbal_y.
  rocket_context.FixInputPort(0, u0);

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
