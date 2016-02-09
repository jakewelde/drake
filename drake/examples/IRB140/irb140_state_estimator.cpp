
#include "RigidBodySystem.h"
#include "LinearSystem.h"
#include "BotVisualizer.h"
#include "drakeAppUtil.h"
#include "drake/examples/IRB140/estimator.hpp"

using namespace std;
using namespace Eigen;
using namespace Drake;

int main(int argc, char* argv[]) {
  if (argc > 1) {
    std::cerr << "Usage: " << argv[0] << " <no args>" << std::endl;
    return 1;
  }

  // todo: consider moving this logic into the RigidBodySystem class so it can be reused
  DrakeJoint::FloatingBaseType floating_base_type = DrakeJoint::FIXED;

  // spawn arm itself
  auto rigid_body_sys = make_shared<RigidBodySystem>("urdf/irb_140.urdf",floating_base_type);
  auto const & tree = rigid_body_sys->getRigidBodyTree();
  // set up its simu options, which will ultimately be used by estimator as
  // approximated dynamics
  SimulationOptions options = default_simulation_options;
  rigid_body_sys->penetration_stiffness = 5000.0;
  rigid_body_sys->penetration_damping = rigid_body_sys->penetration_stiffness/10.0;
  options.initial_step_size = 5e-3;
  options.timeout_seconds = numeric_limits<double>::infinity();

  auto rigid_body_sys_estimator = make_shared<IRB140EstimatorSystem>(rigid_body_sys);

  shared_ptr<lcm::LCM> lcm = make_shared<lcm::LCM>();
  auto visualizer = make_shared<BotVisualizer<RigidBodySystem::StateVector>>(lcm,tree);
  auto sys = cascade(rigid_body_sys_estimator, visualizer);

  

  VectorXd x0(rigid_body_sys_estimator->getNumStates());
  x0.head(tree->num_positions) = tree->getZeroConfiguration();

  runLCM(sys,lcm,0,std::numeric_limits<double>::infinity(),x0,options);
//  simulate(*sys,0,std::numeric_limits<double>::infinity(),x0,options);

  return 0;
}
