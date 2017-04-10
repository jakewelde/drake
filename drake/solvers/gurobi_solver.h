#pragma once

#include <string>

#include "drake/common/drake_copyable.h"
#include "drake/solvers/mathematical_program.h"

namespace drake {
namespace solvers {

class GurobiSolver : public MathematicalProgramSolverInterface {
 public:
  DRAKE_NO_COPY_NO_MOVE_NO_ASSIGN(GurobiSolver)

  GurobiSolver() = default;
  ~GurobiSolver() override = default;

  // This solver is implemented in various pieces depending on if
  // Gurobi was available during compilation.
  bool available() const override;

  // This passes the model and model environment to a
  // an external 
  typedef void (*mipSolCallbackFunction)(const MathematicalProgram&, void *, Eigen::VectorXd&, VectorXDecisionVariable&);
  void addMIPNodeCallback(mipSolCallbackFunction fnc, void * usrdata){
  	mip_node_callback_ = fnc;
  	mip_node_callback_usrdata_ = usrdata;
  }

  SolutionResult Solve(MathematicalProgram& prog) const override;

  SolverType solver_type() const override { return SolverType::kGurobi; }

  std::string SolverName() const override { return "Gurobi"; }

 private:
  mipSolCallbackFunction mip_node_callback_ = NULL;
  void * mip_node_callback_usrdata_ = NULL;
};

}  // end namespace solvers
}  // end namespace drake
