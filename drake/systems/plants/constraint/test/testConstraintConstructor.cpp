#include <iostream>
#include "RigidBodyTree.h"

#include "../RigidBodyConstraint.h"

using namespace std;

int main()
{
  RigidBodyTree * model = new RigidBodyTree("../../../../examples/Atlas/urdf/atlas_minimal_contact.urdf");
  if(!model)
  {
    cerr << "ERROR: Failed to load model"<<endl;
    return -1;
  }
  Eigen::Vector2d tspan;
  tspan<< 0,1;
  Eigen::Vector3d kc1_lb, kc1_ub;
  kc1_lb<<0,0,0;
  kc1_ub<<0,0,0;
  WorldCoMConstraint* kc1 = new WorldCoMConstraint(model,kc1_lb,kc1_ub,tspan);
  printf("construct WorldCoMConstraint\n");
  delete kc1;
  printf("delete WorldCoMConstraint\n");
  delete model;
  return 0;
}
