#include "gtest/gtest.h"

#include "drake/multibody/parsers/urdf_parser.h"
#include "drake/solvers/mathematical_program.h"
#include "drake/solvers/gurobi_solver.h"
#include "drake/common/eigen_matrix_compare.h"
#include "drake/common/eigen_types.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;

namespace drake {
namespace solvers {

GTEST_TEST(GurobiTest, StarBattle) {
  MathematicalProgram prog;
  auto B = prog.NewBinaryVariables(10, 10, "B");

  VectorXd constr = VectorXd::Ones(10);
  for (int i=0; i<10; i++){
    prog.AddLinearEqualityConstraint(constr.transpose(), 2, {B.row(i).transpose()});
    prog.AddLinearEqualityConstraint(constr.transpose(), 2, {B.col(i)});
  }

  MatrixXd board(10, 10);
  board << 0, 0, 1, 1, 1, 1, 1, 1, 2, 2,
           0, 1, 1, 0, 3, 3, 3, 3, 3, 2,
           0, 0, 0, 0, 3, 4, 2, 2, 2, 2,
           3, 3, 3, 3, 3, 4, 4, 5, 5, 5,
           3, 6, 6, 6, 6, 6, 4, 7, 7, 5,
           6, 6, 8, 8, 8, 8, 4, 7, 5, 5,
           6, 8, 8, 9, 9, 8, 4, 7, 7, 7,
           6, 8, 6, 6, 9, 8, 4, 4, 9, 7,
           6, 8, 8, 6, 9, 8, 9, 9, 9, 7,
           6, 6, 6, 6, 9, 9, 9, 7, 7, 7;

  for (int i=0; i<=9; i++){
    VariableListRef vars;
    for (int k=0; k<10; k++){
      for (int l=0; l<10; l++){
        if (board(k, l) == i){
          vars.push_back(B.block<1,1>(k, l));
        }
      }
    }
    prog.AddLinearEqualityConstraint(VectorXd::Ones(vars.size()).transpose(), 2, vars);
  }

  for (int k=0; k<9; k++){
    for (int l=0; l<9; l++){
        prog.AddLinearConstraint(VectorXd::Ones(4).transpose(), 0, 1, {B.block<2, 1>(k, l), B.block<2, 1>(k, l+1)});
    }
  }
  GurobiSolver gurobi_solver;

  prog.SetSolverOption("GUROBI", "PoolSolutions", 1024);
  prog.SetSolverOption("GUROBI", "PoolGap", 0.0);
  prog.SetSolverOption("GUROBI", "PoolSearchMode", 2);

  gurobi_solver.Solve(prog); //prog.Solve();

  std::cout << "Sol: " << std::endl;
  std::cout << prog.GetSolution(B) << std::endl;
}

GTEST_TEST(GurobiTest, Knapsack) {
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

