#include <iostream>

#include "drake/Path.h"
#include "drake/systems/LCMSystem.h"
#include "drake/systems/Simulation.h"
#include "drake/systems/cascade_system.h"
#include "drake/systems/plants/BotVisualizer.h"
#include "drake/systems/plants/RigidBodySystem.h"
#include "drake/util/Polynomial.h"
#include "drake/util/drakeAppUtil.h"

using Drake::RigidBodySystem;
using Drake::BotVisualizer;
using Eigen::VectorXd;

int main(int argc, char* argv[]) {
  // Initializes LCM.
  std::shared_ptr<lcm::LCM> lcm = std::make_shared<lcm::LCM>();

  // Instantiates a rigid body system and adds the valkyrie
  // left hand to it.
  auto rigid_body_sys = std::allocate_shared<RigidBodySystem>(
      Eigen::aligned_allocator<RigidBodySystem>());

  rigid_body_sys->addRobotFromFile(
      Drake::getDrakePath() + "/examples/Valkyrie/urdf/urdf/left_hand_test_bench.urdf",
      DrakeJoint::FIXED);

  // Obtains a reference to the rigid body tree within the rigid body system.
  const auto& tree = rigid_body_sys->getRigidBodyTree();

  // Adds the ground.
  /*
  {
    double box_width = 3;
    double box_depth = 0.2;
    DrakeShapes::Box geom(Eigen::Vector3d(box_width, box_width, box_depth));
    Eigen::Isometry3d T_element_to_link = Eigen::Isometry3d::Identity();
    T_element_to_link.translation() << 0, 0,
        -box_depth / 2.0;  // top of the box is at z = 0
    RigidBody& world = tree->world();
    Eigen::Vector4d color;
    color << 0.9297, 0.7930, 0.6758,
        1;  // was hex2dec({'ee','cb','ad'})'/256 in matlab
    world.addVisualElement(
        DrakeShapes::VisualElement(geom, T_element_to_link, color));
    tree->addCollisionElement(
        RigidBody::CollisionElement(geom, T_element_to_link, &world), world,
        "terrain");
    tree->updateStaticCollisionElements();
  }
*/

  // Sets the stiffness of the ground.
  {
    rigid_body_sys->penetration_stiffness = 3000.0;
    rigid_body_sys->penetration_damping = 0;
  }

  // Instantiates additional systems and cascades them with the rigid body
  // system.
  auto visualizer =
      std::make_shared<BotVisualizer<RigidBodySystem::StateVector>>(lcm, tree);

  auto sys = cascade(rigid_body_sys, visualizer);

  // Obtains an initial state of the simulation.
  VectorXd x0 = VectorXd::Zero(rigid_body_sys->getNumStates());
  x0.head(tree->number_of_positions()) = tree->getZeroConfiguration();

  // Specifies the simulation options.
  Drake::SimulationOptions options;
  options.realtime_factor = 0;  // As fast as possible.
  options.initial_step_size = 0.001;

  // Prevents exception from being thrown when simulation runs slower than real
  // time, which it most likely will given the small step size.
  options.warn_real_time_violation = true;

  // Instantates a variable that specifies the duration of the simulation.
  // The default value is 5 seconds.
  double duration = 5.0;

  // Searches through the command line looking for a "--duration" flag followed
  // by a floating point number that specifies a custom duration.
  for (int ii = 1; ii < argc; ++ii) {
    if (std::string(argv[ii]) == "--duration") {
      if (++ii == argc) {
        throw std::runtime_error(
            "ERROR: Command line option \"--duration\" is not followed by a "
            "value!");
      }
      duration = atof(argv[ii]);
    }
  }

  // Starts the simulation.
  const double kStartTime = 0;
  Drake::simulate(*sys.get(), kStartTime, duration, x0, options);
}
