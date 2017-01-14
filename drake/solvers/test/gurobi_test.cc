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
/*
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
*/

GTEST_TEST(GurobiTest, Knapsack) {
  MathematicalProgram prog;
  auto B_knapsack = prog.NewBinaryVariables(10, "B_knapsack");

  VectorXd knapsack_costs(10);
  VectorXd knapsack_obj_coeffs(10);
  knapsack_costs << 16, 16,  8,  8, 4, 4, 2, 2, 1, 1;
  knapsack_obj_coeffs << 32, 32, 15, 15, 6, 6, 1, 1, 1, 1;
  double total_budget = 33;
  prog.AddLinearConstraint(knapsack_costs.transpose(), -std::numeric_limits<double>::infinity(), total_budget, {B_knapsack});
  prog.AddLinearCost(-1.0*knapsack_obj_coeffs, {B_knapsack});

  GurobiSolver gurobi_solver;

  prog.SetSolverOption("GUROBI", "PoolSolutions", 1024);
  prog.SetSolverOption("GUROBI", "PoolGap", 0.10);
  prog.SetSolverOption("GUROBI", "PoolSearchMode", 2);

  //prog.SetSolverOption("GUROBI", "LogToConsole", 1);
  //prog.SetSolverOption("GUROBI", "OutputFlag", 1);

  gurobi_solver.Solve(prog);

  // should find 21 solutions
  EXPECT_EQ(prog.get_num_solutions(), 21);

  // Can verify values TODO(gizatt)

  printf("done!\n");
}

extern "C" {
#include "gurobi_c.h"
}
GTEST_TEST(GurobiTest, KnapsackDirectGurobiInterface) {
  GRBenv   *env   = NULL;
  GRBenv   *menv  = NULL;
  GRBmodel *model = NULL;
  int       error = 0;
  char      buffer[1000];
  int e, status, nSolutions, prlen;
  double objval, *cval = NULL;
  int *cind = NULL;

  /* Sample data */
  const int groundSetSize = 10;
  double objCoef[10] =
    {32, 32, 15, 15, 6, 6, 1, 1, 1, 1};
  double knapsackCoef[10] =
    {16, 16,  8,  8, 4, 4, 2, 2, 1, 1};
  double Budget = 33;

  /* Create environment */
  error = GRBloadenv(&env, "poolsearch_c.log");
  if (error) goto QUIT;

  /* Create initial model */
  error = GRBnewmodel(env, &model, "poolsearch_c", groundSetSize, NULL,
                      NULL, NULL, NULL, NULL);
  if (error) goto QUIT;

  /* get model environment */
  menv = GRBgetenv(model);
  if (!menv) {
    fprintf(stderr, "Error: could not get model environment\n");
    goto QUIT;
  }

  /* set objective function */
  error = GRBsetdblattrarray(model, "Obj", 0, groundSetSize, objCoef);
  if (error) goto QUIT;

  /* set variable types and names */
  for (e = 0; e < groundSetSize; e++) {
    sprintf(buffer, "El%d", e);
    error = GRBsetcharattrelement(model, "VType", e, GRB_BINARY);
    if (error) goto QUIT;

    error = GRBsetstrattrelement(model, "VarName", e, buffer);
    if (error) goto QUIT;
  }

  /* Make space for constraint data */
  cind = (int *) malloc(sizeof(int) * groundSetSize);
  if (!cind) goto QUIT;
  for (e = 0; e < groundSetSize; e++)
    cind[e] = e;

  /* Constraint: limit total number of elements to be picked to be at most
   * Budget */
  sprintf (buffer, "Budget");
  error = GRBaddconstr(model, groundSetSize, cind, knapsackCoef,
                       GRB_LESS_EQUAL, Budget, buffer);
  if (error) goto QUIT;

  /* set global sense for ALL objectives */
  error = GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE, GRB_MAXIMIZE);
  if (error) goto QUIT;

  /* Limit how many solutions to collect */
  error = GRBsetintparam(menv, GRB_INT_PAR_POOLSOLUTIONS, 1024);
  if (error) goto QUIT;

  /* Limit the search space by setting a gap for the worst possible solution that will be accepted */
  error = GRBsetdblparam(menv, GRB_DBL_PAR_POOLGAP, 0.10);
  if (error) goto QUIT;

  /* do a systematic search for the k-best solutions */
  error = GRBsetintparam(menv, GRB_INT_PAR_POOLSEARCHMODE, 2);
  if (error) goto QUIT;

  /* save problem */
  error = GRBwrite(model, "poolsearch_c.lp");
  if (error) goto QUIT;
  error = GRBwrite(model, "poolsearch_c.mps");
  if (error) goto QUIT;

  /* Optimize */
  error = GRBoptimize(model);
  if (error) goto QUIT;

  /* Status checking */
  error = GRBgetintattr(model, "Status", &status);
  if (error) goto QUIT;

  if (status == GRB_INF_OR_UNBD ||
      status == GRB_INFEASIBLE  ||
      status == GRB_UNBOUNDED     ) {
    printf("The model cannot be solved "
           "because it is infeasible or unbounded\n");
    goto QUIT;
  }
  if (status != GRB_OPTIMAL) {
    printf("Optimization was stopped with status %d\n", status);
    goto QUIT;
  }

  /* make space for optimal solution */
  cval = (double *) malloc(sizeof(double) * groundSetSize);
  if (!cval) goto QUIT;

  /* Print best selected set */
  error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, groundSetSize, cval);
  if (error) goto QUIT;

  printf("Selected elements in best solution:\n\t");
  for (e = 0; e < groundSetSize; e++) {
    if (cval[e] < .9) continue;
    printf("El%d ", e);
  }

  /* print number of solutions stored */
  error = GRBgetintattr(model, GRB_INT_ATTR_SOLCOUNT, &nSolutions);
  if (error) goto QUIT;
  printf("\nNumber of solutions found: %d\nValues:", nSolutions);

  /* print objective values of alternative solutions */
  prlen = 0;
  for (e = 0; e < nSolutions; e++) {
    error = GRBsetintparam(menv, GRB_INT_PAR_SOLUTIONNUMBER, e);
    if (error) goto QUIT;

    error = GRBgetdblattr(model, GRB_DBL_ATTR_POOLOBJVAL, &objval);
    if (error) goto QUIT;

    prlen += printf(" %g", objval);
    if (prlen >= 75 && e+1 < nSolutions) {
      prlen = printf("\n    ");
    }
  }
  printf("\n");

  /* print fourth best set if available */
  if (nSolutions >= 4) {
    error = GRBsetintparam(menv, GRB_INT_PAR_SOLUTIONNUMBER, 3);
    if (error) goto QUIT;

    /* get the solution vector */
    error = GRBgetdblattrarray(model, GRB_DBL_ATTR_XN, 0, groundSetSize, cval);
    if (error) goto QUIT;

    printf("Selected elements in fourth best solution:\n\t");
    for (e = 0; e < groundSetSize; e++) {
      if (cval[e] < .9) continue;
      printf("El%d ", e);
    }
    printf("\n");
  }

  QUIT:
    if (model != NULL) GRBfreemodel(model);
    if (env != NULL)   GRBfreeenv(env);
    if (cind != NULL)  free(cind);
    if (cval != NULL)  free(cval);
}



}  // namespace solvers
}  // namespace drake

