#ifndef DRAKE_IRB140_ESTIMATOR_H
#define DRAKE_IRB140_ESTIMATOR_H

#include "drake/systems/System.h"
#include "drake/systems/LCMSystem.h"
#include <lcm/lcm-cpp.hpp>
#include "drake/solvers/Optimization.h"
#include "drake/systems/plants/RigidBodySystem.h"
#include "drake/systems/plants/RigidBodyFrame.h"
#include "KinematicsCache.h"
#include "drake/drakeRBSystem_export.h"
#include "lcmtypes/drake/lcmt_robot_state.hpp"
#include "bot_core/planar_lidar_t.hpp"
#include "bot_core/rigid_transform_t.hpp"
#include "lcmtypes/drake/lcmt_point_cloud.hpp"

namespace tinyxml2 {
  class XMLElement;
}

namespace Drake {


template <typename ScalarType = double>
class PlanarLidarMsg {
public:
  typedef bot_core::planar_lidar_t LCMMessageType;
  static std::string channel() { return "SCAN"; };

  PlanarLidarMsg(void) {};
  template <typename Derived>
  PlanarLidarMsg(const Eigen::MatrixBase<Derived>& x) {};

  template <typename Derived>
  PlanarLidarMsg& operator=(const Eigen::MatrixBase<Derived>& x) {
    printf("WTF is this method even\n");
    return *this;
  }

  template <typename Derived>
  int size(void){
    return ranges.rows();
  }

  friend Eigen::Vector3d toEigen(const PlanarLidarMsg<ScalarType>& scan) {
    Eigen::VectorXd x;
    for (int i=0; i<scan.nranges; i++)
      x << scan.ranges[i];
    return x;
  }

  friend std::string getCoordinateName(const PlanarLidarMsg<ScalarType>& scan, unsigned int index) {
    switch (index) {
      case 0: return "nranges";
      case 1: return "ranges";
      case 2: return "nintensities";
      case 3: return "intensities";
      case 4: return "rad0";
      case 5: return "radstep";
    }
    return "error";
  }
  const static int RowsAtCompileTime = Eigen::Dynamic;

  std::vector<ScalarType> ranges;
};

bool decode(const bot_core::planar_lidar_t& msg, double& t, PlanarLidarMsg<double>& x) {
  t = double(msg.utime)/1000.0;
  printf("UNHANDLED DECODE.\n");
  return true;
}


template <typename ScalarType = double>
class PointCloudMsg {
public:
  typedef drake::lcmt_point_cloud LCMMessageType;
  static std::string channel() { return "SCAN"; };

  PointCloudMsg(void) {};
  template <typename Derived>
  PointCloudMsg(const Eigen::MatrixBase<Derived>& x) {};

  template <typename Derived>
  PointCloudMsg& operator=(const Eigen::MatrixBase<Derived>& x) {
    printf("WTF is this method even\n");
    return *this;
  }

  template <typename Derived>
  int size(void){
    return points.rows();
  }

  friend Eigen::Matrix3Xd toEigen(const PointCloudMsg<ScalarType>& scan) {
    Eigen::Matrix3Xd x(scan.npoints);
    for (int i=0; i<scan.npoints; i++)
      for (int j=0; j<3; j++)
        x(i, j) = scan.points[i][j];

    return x;
  }

  friend std::string getCoordinateName(const PointCloudMsg<ScalarType>& scan, unsigned int index) {
    switch (index) {
      case 0: return "npoints";
      case 1: return "points";
    }
    return "error";
  }
  const static int RowsAtCompileTime = Eigen::Dynamic;

  Eigen::Matrix3Xd points;
};

bool decode(const drake::lcmt_point_cloud& msg, double& t, PlanarLidarMsg<double>& x) {
  t = double(msg.timestamp)/1000.0;
  printf("UNHANDLED DECODE.\n");
  return true;
}

/** IRB140EstimatorSystem
 * @brief Given a rigid body system describing the IRB140, estimates its states given
 *   inputs to the system and sensor data (for now only actual being piped in).
 *
 *   u---->[   sys   ]------> yout
 *   |      /|\    |           |
 *   |       |     |           |
 *   |       ---x---           \/
 *   |----------------------->[ estimator ]--> x_estimated
 *
 */
  class IRB140EstimatorSystem {
  public:
    template <typename ScalarType> using RigidBodySystemInputVector = Eigen::Matrix<ScalarType,Eigen::Dynamic,1>;
    template <typename ScalarType> using InputVector = RigidBodySystemInputVector<ScalarType>;
    template <typename ScalarType> using StateVector = Eigen::Matrix<ScalarType,Eigen::Dynamic,1>;
    template <typename ScalarType> using OutputVector = Eigen::Matrix<ScalarType,Eigen::Dynamic,1>;

    IRB140EstimatorSystem(const std::shared_ptr<RigidBodySystem>& plant_in, 
                          const std::shared_ptr<RigidBodyFrame>& lidarFrame_in) : 
                                            plant(plant_in), lidarFrame(lidarFrame_in)  { 
      lcm = std::unique_ptr<lcm::LCM>(new lcm::LCM());

      if (lcm->good()){
        this->setupSubscriptions();
        printf("Init with LCM!\n"); 
      } else {
        printf("Could not init LCM!\n");
      }
    };
    virtual ~IRB140EstimatorSystem() {};

    const std::shared_ptr<RigidBodySystem>& getRigidBodySystem(void) { return plant; }
    void setupSubscriptions();
    size_t getNumStates() const { return plant->getNumStates(); }
    size_t getNumInputs() const { return plant->getNumInputs() + plant->getNumOutputs(); }
    size_t getNumOutputs() const { return plant->getNumStates(); }

    /** dynamics
     * wooooooo fun part
     */
    StateVector<double> dynamics(const double& t, const StateVector<double>& x, const InputVector<double>& u) const;

    template <typename ScalarType>
    OutputVector<ScalarType> output(const ScalarType& t, const StateVector<ScalarType>& x, const InputVector<ScalarType>& u) const {
      return x;
    }

    bool isTimeVarying() const  { return false; }
    bool isDirectFeedthrough() const { return false; }

    friend DRAKERBSYSTEM_EXPORT StateVector<double> getInitialState(const IRB140EstimatorSystem& sys);

    void handlePlanarLidarMsg(const lcm::ReceiveBuffer* rbuf,
                              const std::string& chan,
                              const bot_core::planar_lidar_t* msg);
    void handlePointCloudMsg( const lcm::ReceiveBuffer* rbuf,
                              const std::string& chan,
                              const drake::lcmt_point_cloud* msg);
    void handleSpindleFrameMsg(const lcm::ReceiveBuffer* rbuf,
                               const std::string& chan,
                               const bot_core::rigid_transform_t* msg);

  private:
    std::shared_ptr<RigidBodySystem> plant;
    std::unique_ptr<lcm::LCM> lcm;

    bool haveLidarPts = false;
    Eigen::Matrix3Xd latestLidarPts; // todo: replace with rigidbodyframe
    RigidBodyFrame lidarFrame;

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  };
} // end namespace Drake

#endif
