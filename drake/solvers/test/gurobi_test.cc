#include "gtest/gtest.h"

#include "drake/multibody/parsers/urdf_parser.h"
#include "drake/solvers/mathematical_program.h"
#include "drake/solvers/gurobi_solver.h"
#include "drake/common/eigen_matrix_compare.h"
#include "drake/common/eigen_types.h"

using Eigen::VectorXd;

namespace drake {
namespace solvers {

GTEST_TEST(GurobiTest, ) {
  MathematicalProgram prog;
  auto B_knapsack = prog.NewBinaryVariables(10, "B_knapsack");

  VectorXd knapsack_costs(10);
  VectorXd knapsack_obj_coeffs(10);
  knapsack_costs << 16, 16,  8,  8, 4, 4, 2, 2, 1, 1;
  knapsack_obj_coeffs << 32, 32, 15, 15, 6, 6, 1, 1, 1, 1;
  double total_budget = 33;
  prog.AddLinearConstraint(knapsack_costs.transpose(), 0, total_budget, {B_knapsack});
  prog.AddLinearCost(-1.0*knapsack_obj_coeffs, {B_knapsack});

  GurobiSolver gurobi_solver;

  prog.SetSolverOption("GUROBI", "PoolSolutions", 1024);
  prog.SetSolverOption("GUROBI", "PoolGap", 0.10);
  prog.SetSolverOption("GUROBI", "PoolSearchMode", 2);

  prog.SetSolverOption("GUROBI", "LogToConsole", 1);
  prog.SetSolverOption("GUROBI", "OutputFlag", 1);

  gurobi_solver.Solve(prog);

  printf("done!\n");
}

}  // namespace solvers
}  // namespace drake

