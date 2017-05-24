#pragma once

#include <string>

#include <Eigen/Core>

#include "drake/common/drake_copyable.h"
#include "drake/solvers/mathematical_program.h"

namespace drake {
namespace solvers {

class MosekSolver : public MathematicalProgramSolverInterface {
 public:
  DRAKE_NO_COPY_NO_MOVE_NO_ASSIGN(MosekSolver)

  MosekSolver() : MathematicalProgramSolverInterface(SolverType::kMosek) {}

  /**
   * Defined true if Mosek was included during compilation, false otherwise.
   */
  bool available() const override;

  /**
   * Sets the branching priority for a set of variables.
   */

  void set_branch_priority(const VariableRefList& vars, int priority);

  SolutionResult Solve(MathematicalProgram& prog) const override;
 
 private:
  std::list<std::pair<drake::symbolic::Variable, int>>
    branch_priority_settings_;
};

}  // namespace solvers
}  // namespace drake
