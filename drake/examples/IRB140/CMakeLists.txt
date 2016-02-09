include_directories(${PROJECT_SOURCE_DIR}/core)
include_directories(${PROJECT_SOURCE_DIR}/solvers)
include_directories(${PROJECT_SOURCE_DIR}/systems)
include_directories(${PROJECT_SOURCE_DIR}/systems/frames)
include_directories(${PROJECT_SOURCE_DIR}/systems/plants)
include_directories(${PROJECT_SOURCE_DIR}/systems/plants/joints)
include_directories(${PROJECT_SOURCE_DIR}/systems/controllers)
include_directories(${PROJECT_SOURCE_DIR}/util)

if (LCM_FOUND)
  add_library_with_exports(LIB_NAME drakeIRB140EstimatorSystem SOURCE_FILES estimator.cpp)
  target_link_libraries(drakeIRB140EstimatorSystem drakeRBSystem drakeRBM drakeLCMSystem drakeOptimization)

  include_directories(${PROJECT_SOURCE_DIR}/pod-build/lcmgen)
  include_directories(${PROJECT_SOURCE_DIR}/../build/include/lcmtypes)
  include_directories(${PROJECT_SOURCE_DIR}/../../build/include/lcmtypes)
  include_directories(${PROJECT_SOURCE_DIR}/../../../build/include/lcmtypes)
  add_executable(irb140StateEstimator irb140_state_estimator.cpp)
  target_link_libraries(irb140StateEstimator drakeIRB140EstimatorSystem drakeRBSystem drakeLCMSystem)


endif (LCM_FOUND)
