
#include <stdexcept>
#include "drake/examples/IRB140/estimator.hpp"
#include "urdfParsingUtil.h"

using namespace std;
using namespace Eigen;
using namespace Drake;
using namespace tinyxml2;

void IRB140EstimatorSystem::setupSubscriptions(){
  //lcm->subscribe("SCAN", &IRB140EstimatorSystem::handlePointCloud, this);
  lcm->subscribe("SCAN", &IRB140EstimatorSystem::handlePlanarLidarMsg, this);
  lcm->subscribe("PRE_SPINDLE_TO_POST_SPINDLE", &IRB140EstimatorSystem::handleSpindleFrameMsg, this);
}

void IRB140EstimatorSystem::handlePlanarLidarMsg(const lcm::ReceiveBuffer* rbuf,
                           const std::string& chan,
                           const bot_core::planar_lidar_t* msg){
  printf("Received scan on channel %s\n", chan.c_str());
  // transform according 
}

void IRB140EstimatorSystem::handlePointCloudMsg(const lcm::ReceiveBuffer* rbuf,
                           const std::string& chan,
                           const drake::lcmt_point_cloud* msg){
  printf("Received scan pts on channel %s\n", chan.c_str());
  // todo: transform them all by the lidar frame
}

void IRB140EstimatorSystem::handleSpindleFrameMsg(const lcm::ReceiveBuffer* rbuf,
                           const std::string& chan,
                           const bot_core::rigid_transform_t* msg){
  printf("Received transform on channel %s\n", chan.c_str());
  cout << msg->trans << "," << msg->quat << endl;
  // todo: transform them all by the lidar frame
}


IRB140EstimatorSystem::StateVector<double> IRB140EstimatorSystem::dynamics(const double& t, const IRB140EstimatorSystem::StateVector<double>& x, const IRB140EstimatorSystem::InputVector<double>& u) const {
  using namespace std;
  using namespace Eigen;

  /* Once we figure out how dynamic-size input frames work, we'll
  get this info in from the input from from runLCM. However, for now,
  I'll query LCM manually to get pointcloud info. */
  lcm->handle();

  //cout << x << endl;
  StateVector<double> dot(this->getNumStates());
  return dot;
}

DRAKERBSYSTEM_EXPORT IRB140EstimatorSystem::StateVector<double> Drake::getInitialState(const IRB140EstimatorSystem& sys) {

  VectorXd x0(sys.getNumStates());
  default_random_engine generator;
  x0 << getInitialState(sys.plant);
  return x0;
}