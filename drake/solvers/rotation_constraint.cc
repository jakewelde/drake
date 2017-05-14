/* clang-format off */
#include "drake/solvers/rotation_constraint.h"
#include "drake/solvers/rotation_constraint_internal.h"
/* clang-format on */

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>

#include "drake/common/symbolic_expression.h"
#include "drake/math/cross_product.h"

using std::numeric_limits;
using drake::symbolic::Expression;
using Eigen::VectorXd;
using Eigen::RowVectorXd;
using Eigen::MatrixXi;

namespace drake {
namespace solvers {

namespace internal {
  MatrixXi CalculateReflectedGrayCodes(int k){
    int num_codes = k == 0 ? 0 : 1 << k;
    MatrixXi return_codes(num_codes, k);
    return_codes.setZero();
    for (int i = 0; i < num_codes; i++) {
      for (int j = 0; j < k; j++) {
        // The jth digit cycles
        // 2^j times off, then 2^j times on.
        // The jth digit goes through a quarter-period
        // (2^j / 2) off first.
        int quarter_period_time = 1 << (j);
        int half_period_time = 2 << (j);
        int full_period_time = 4 << (j);
        if ((i + quarter_period_time) % full_period_time >= half_period_time) {
          return_codes(i, j) = 1;
        }
      }
    }
    return return_codes;
  }
}

/* 
   Encode a Special Ordered Set of the 2nd type...
      i.e., only 2 adjacent members of lambda should
      be nonzero at a time.

   Given a set of expressions V expressions Lambda,
   and a set of N binary variables Bin, and a V x N
   Gray Coding that maps activatings of Bin to nonzero-ness
   of Lambda.

   Binary variable i being inactive means that 

*/
void AddSOS2EncodingConstraints(MathematicalProgram* prog,
    Eigen::Matrix<Expression, Eigen::Dynamic, 1> & lambda,
    MatrixXDecisionVariable& bin,
    MatrixXi coding) {
  for (int i = 0; i < bin.rows(); i++) {
    Expression lb = coding(0, i) * lambda(0);
    Expression ub = coding(0, i) * lambda(0);
    for (int v = 1; v < lambda.rows()-1; v++) {
      lb += fmin(coding(v, i), coding(v-1, i))*lambda(v);
      ub += fmax(coding(v, i), coding(v-1, i))*lambda(v);
    }
    lb += coding(lambda.rows()-2, i)*lambda(lambda.rows()-1);
    ub += coding(lambda.rows()-2, i)*lambda(lambda.rows()-1);
    prog->AddLinearConstraint(lb <= bin(i));
    prog->AddLinearConstraint(ub >= bin(i));
  }
}

void Add2DLogarithmicMcCormickEnvelope(MathematicalProgram* prog, 
                              drake::symbolic::Variable& w,
                              drake::symbolic::Variable& x, 
                              drake::symbolic::Variable& y, 
                              std::string corename,
                              double xL, 
                              double xH, 
                              double yL, 
                              double yH, 
                              int M_x,
                              int M_y){
  auto x_bounds = VectorXd::LinSpaced(M_x+1, xL, xH);
  auto y_bounds = VectorXd::LinSpaced(M_y+1, yL, yH);

  // Disjunction of envelope regions in convex hull form:
  MatrixXDecisionVariable lambda = prog->NewContinuousVariables(M_x+1, M_y+1, corename + "_lambda");
  // Lambdas are affine coefficients.
  prog->AddBoundingBoxConstraint(0., 1., lambda);
  prog->AddLinearConstraint(
    RowVectorXd::Ones(M_y+1)*((RowVectorXd::Ones(M_x+1) * lambda).transpose())
    == VectorXd::Ones(1));
  
  // Lambdas dictate affine combination of bounds to generate
  // values for x and y; and can be used similarly to construct
  // the approximate value for w = x * y.
  Expression x_affine_combination = 0.0;
  Expression y_affine_combination = 0.0;
  Expression w_affine_combination = 0.0;
  for (int i = 0; i < M_x+1; i++){
    for (int j = 0; j < M_y+1; j++){
      x_affine_combination += lambda(i, j) * x_bounds(i);
      y_affine_combination += lambda(i, j) * y_bounds(j);
      w_affine_combination += lambda(i, j) * x_bounds(i) * y_bounds(j);
    }
  }
  prog->AddLinearConstraint(x == x_affine_combination);
  prog->AddLinearConstraint(y == y_affine_combination);
  prog->AddLinearConstraint(w == w_affine_combination);


  // For each X bin...
  if (M_x > 1) {
    // Form an expression that corresponds to what y is selected (by inspecting the lambdas
    // in this column).
    Eigen::Matrix<Expression, Eigen::Dynamic, 1> y_col = RowVectorXd::Ones(M_y+1) * lambda.transpose();
    int k = ceil(log2(M_x));
    auto coding = internal::CalculateReflectedGrayCodes(k);
    coding.conservativeResize(M_x, k);
    // And generate a corresponding number of binary variables indicating membership in
    // those regions.
    auto bin = prog->NewBinaryVariables(k, 1, corename + "_x_bin");
    AddSOS2EncodingConstraints(prog, y_col, bin, coding);
  }

  // For each Y bin...
  if (M_y > 1) {
    // Form an expression that corresponds to what x is selected (by inspecting the lambdas
    // in this column).
    Eigen::Matrix<Expression, Eigen::Dynamic, 1> x_col = RowVectorXd::Ones(M_x+1) * lambda;
    int k = ceil(log2(M_y));
    auto coding = internal::CalculateReflectedGrayCodes(k);
    coding.conservativeResize(M_y, k);
    // And generate a corresponding number of binary variables indicating membership in
    // those regions.
    auto bin = prog->NewBinaryVariables(k, 1, corename + "_y_bin");
    AddSOS2EncodingConstraints(prog, x_col, bin, coding);
  }
}

MatrixDecisionVariable<3, 3> NewRotationMatrixVars(MathematicalProgram* prog,
                                                   const std::string& name) {
  MatrixDecisionVariable<3, 3> R = prog->NewContinuousVariables<3, 3>(name);

  // Forall i,j, -1 <= R(i,j) <=1.
  prog->AddBoundingBoxConstraint(-1, 1, R);

  // -1 <= trace(R) <= 3.
  // Proof sketch:
  //   orthonormal => |lambda_i|=1.
  //   R is real => eigenvalues either real or appear in complex conj pairs.
  //   Case 1: All real (lambda_i \in {-1,1}).
  //     det(R)=lambda_1*lambda_2*lambda_3=1 => lambda_1=lambda_2, lambda_3=1.
  //   Case 2: Two imaginary, pick conj(lambda_1) = lambda_2.
  //    => lambda_1*lambda_2 = 1.  =>  lambda_3 = 1.
  //    and also => lambda_1 + lambda_2 = 2*Re(lambda_1) \in [-2,2].
  prog->AddLinearConstraint(Eigen::RowVector3d::Ones(), -1, 3, R.diagonal());
  return R;
}

void AddBoundingBoxConstraintsImpliedByRollPitchYawLimits(
    MathematicalProgram* prog,
    const Eigen::Ref<const MatrixDecisionVariable<3, 3>>& R,
    RollPitchYawLimits limits) {
  // Based on the RPY to Rotation Matrix conversion:
  // [ cp*cy, cy*sp*sr - cr*sy, sr*sy + cr*cy*sp]
  // [ cp*sy, cr*cy + sp*sr*sy, cr*sp*sy - cy*sr]
  // [   -sp,            cp*sr,            cp*cr]
  // where cz = cos(z) and sz = sin(z), and using
  //  kRoll_NegPI_2_to_PI_2 = 1 << 1,   // => cos(r)>=0
  //  kRoll_0_to_PI = 1 << 2,           // => sin(r)>=0
  //  kPitch_NegPI_2_to_PI_2 = 1 << 3,  // => cos(p)>=0
  //  kPitch_0_to_PI = 1 << 4,          // => sin(p)>=0
  //  kYaw_NegPI_2_to_PI_2 = 1 << 5,    // => cos(y)>=0
  //  kYaw_0_to_PI = 1 << 6,            // => sin(y)>=0

  if ((limits & kPitch_NegPI_2_to_PI_2) && (limits & kYaw_NegPI_2_to_PI_2))
    prog->AddBoundingBoxConstraint(0, 1, R(0, 0));

  if ((limits & kPitch_NegPI_2_to_PI_2) && (limits & kYaw_0_to_PI))
    prog->AddBoundingBoxConstraint(0, 1, R(1, 0));

  if (limits & kPitch_0_to_PI) prog->AddBoundingBoxConstraint(-1, 0, R(2, 0));

  if ((limits & kRoll_NegPI_2_to_PI_2) && (limits & kYaw_NegPI_2_to_PI_2) &&
      (limits & kPitch_0_to_PI) && (limits & kRoll_0_to_PI) &&
      (limits & kYaw_0_to_PI))
    prog->AddBoundingBoxConstraint(0, 1, R(1, 1));

  if ((limits & kPitch_NegPI_2_to_PI_2) && (limits & kRoll_0_to_PI))
    prog->AddBoundingBoxConstraint(0, 1, R(2, 1));

  if ((limits & kRoll_0_to_PI) && (limits & kYaw_0_to_PI) &&
      (limits & kRoll_NegPI_2_to_PI_2) && (limits & kYaw_NegPI_2_to_PI_2) &&
      (limits & kPitch_0_to_PI))
    prog->AddBoundingBoxConstraint(0, 1, R(0, 2));

  if ((limits & kPitch_NegPI_2_to_PI_2) && (limits & kRoll_NegPI_2_to_PI_2))
    prog->AddBoundingBoxConstraint(0, 1, R(2, 2));
}

void AddBoundingBoxConstraintsImpliedByRollPitchYawLimitsToBinary(
    MathematicalProgram* prog,
    const Eigen::Ref<const MatrixDecisionVariable<3, 3>>& B,
    RollPitchYawLimits limits) {
  if ((limits & kPitch_NegPI_2_to_PI_2) && (limits & kYaw_NegPI_2_to_PI_2))
    prog->AddBoundingBoxConstraint(1, 1, B(0, 0));

  if ((limits & kPitch_NegPI_2_to_PI_2) && (limits & kYaw_0_to_PI))
    prog->AddBoundingBoxConstraint(1, 1, B(1, 0));

  if (limits & kPitch_0_to_PI) prog->AddBoundingBoxConstraint(0, 0, B(2, 0));

  if ((limits & kRoll_NegPI_2_to_PI_2) && (limits & kYaw_NegPI_2_to_PI_2) &&
      (limits & kPitch_0_to_PI) && (limits & kRoll_0_to_PI) &&
      (limits & kYaw_0_to_PI))
    prog->AddBoundingBoxConstraint(1, 1, B(1, 1));

  if ((limits & kPitch_NegPI_2_to_PI_2) && (limits & kRoll_0_to_PI))
    prog->AddBoundingBoxConstraint(1, 1, B(2, 1));

  if ((limits & kRoll_0_to_PI) && (limits & kYaw_0_to_PI) &&
      (limits & kRoll_NegPI_2_to_PI_2) && (limits & kYaw_NegPI_2_to_PI_2) &&
      (limits & kPitch_0_to_PI))
    prog->AddBoundingBoxConstraint(1, 1, B(0, 2));

  if ((limits & kPitch_NegPI_2_to_PI_2) && (limits & kRoll_NegPI_2_to_PI_2))
    prog->AddBoundingBoxConstraint(1, 1, B(2, 2));
}

void AddRotationMatrixSpectrahedralSdpConstraint(
    MathematicalProgram* prog,
    const Eigen::Ref<const MatrixDecisionVariable<3, 3>>& R) {
  // TODO(russt): Clean this up using symbolic expressions!
  Eigen::Matrix4d F0 = Eigen::Matrix4d::Identity();
  Eigen::Matrix4d F11 = Eigen::Matrix4d::Zero();
  F11(0, 0) = -1;
  F11(1, 1) = 1;
  F11(2, 2) = 1;
  F11(3, 3) = -1;
  Eigen::Matrix4d F21 = Eigen::Matrix4d::Zero();
  F21(0, 2) = -1;
  F21(1, 3) = 1;
  F21(2, 0) = -1;
  F21(3, 1) = 1;
  Eigen::Matrix4d F31 = Eigen::Matrix4d::Zero();
  F31(0, 1) = 1;
  F31(1, 0) = 1;
  F31(2, 3) = 1;
  F31(3, 2) = 1;
  Eigen::Matrix4d F12 = Eigen::Matrix4d::Zero();
  F12(0, 2) = 1;
  F12(1, 3) = 1;
  F12(2, 0) = 1;
  F12(3, 1) = 1;
  Eigen::Matrix4d F22 = Eigen::Matrix4d::Zero();
  F22(0, 0) = -1;
  F22(1, 1) = -1;
  F22(2, 2) = 1;
  F22(3, 3) = 1;
  Eigen::Matrix4d F32 = Eigen::Matrix4d::Zero();
  F32(0, 3) = 1;
  F32(1, 2) = -1;
  F32(2, 1) = -1;
  F32(3, 0) = 1;
  Eigen::Matrix4d F13 = Eigen::Matrix4d::Zero();
  F13(0, 1) = 1;
  F13(1, 0) = 1;
  F13(2, 3) = -1;
  F13(3, 2) = -1;
  Eigen::Matrix4d F23 = Eigen::Matrix4d::Zero();
  F23(0, 3) = 1;
  F23(1, 2) = 1;
  F23(2, 1) = 1;
  F23(3, 0) = 1;
  Eigen::Matrix4d F33 = Eigen::Matrix4d::Zero();
  F33(0, 0) = 1;
  F33(1, 1) = -1;
  F33(2, 2) = 1;
  F33(3, 3) = -1;

  prog->AddLinearMatrixInequalityConstraint(
      {F0, F11, F21, F31, F12, F22, F32, F13, F23, F33},
      {R.col(0), R.col(1), R.col(2)});
}

namespace {

void AddOrthogonalConstraint(
    MathematicalProgram* prog,
    const Eigen::Ref<const VectorDecisionVariable<3>>& v1,
    const Eigen::Ref<const VectorDecisionVariable<3>>& v2) {
  // We do this by introducing
  //   |v1+v2|^2 = v1'v1 + 2v1'v2 + v2'v2 <= 2
  //   |v1-v2|^2 = v1'v1 - 2v1'v2 + v2'v2 <= 2
  // This is tight when v1'v1 = 1 and v2'v2 = 1.

  // TODO(russt): Consider generalizing this to |v1+alpha*v2|^2 <= 1+alpha^2,
  // for any real-valued alpha.  When |R1|<|R2|<=1 or |R2|<|R1|<=1,
  // different alphas represent different constraints.

  Eigen::Matrix<double, 5, 6> A;
  Eigen::Matrix<double, 5, 1> b;

  // |v1+v2|^2 <= 2
  // Implemented as a Lorenz cone using z = [ sqrt(2); v1+v2 ].
  Vector4<symbolic::Expression> z;
  z << std::sqrt(2), v1 + v2;
  prog->AddLorentzConeConstraint(z);

  // |v1-v2|^2 <= 2
  // Implemented as a Lorenz cone using z = [ sqrt(2); v1-v2 ].
  z.tail<3>() = v1 - v2;
  prog->AddLorentzConeConstraint(z);
}

}  // namespace

void AddRotationMatrixOrthonormalSocpConstraint(
    MathematicalProgram* prog,
    const Eigen::Ref<const MatrixDecisionVariable<3, 3>>& R) {
  // All columns should be unit length (but we can only write Ri'Ri<=1),
  // implemented as a rotated Lorenz cone with z = Ax+b = [1;1;R.col(i)].
  Eigen::Matrix<double, 5, 3> A = Eigen::Matrix<double, 5, 3>::Zero();
  A.bottomRows<3>() = Eigen::Matrix3d::Identity();
  Eigen::Matrix<double, 5, 1> b;
  b << 1, 1, 0, 0, 0;
  for (int i = 0; i < 3; i++) {
    prog->AddRotatedLorentzConeConstraint(A, b, R.col(i));
    prog->AddRotatedLorentzConeConstraint(A, b, R.row(i).transpose());
  }

  AddOrthogonalConstraint(prog, R.col(0), R.col(1));  // R0'*R1 = 0.
  AddOrthogonalConstraint(prog, R.col(1), R.col(2));  // R1'*R2 = 0.
  AddOrthogonalConstraint(prog, R.col(0), R.col(2));  // R0'*R2 = 0.

  // Same for the rows
  AddOrthogonalConstraint(prog, R.row(0).transpose(), R.row(1).transpose());
  AddOrthogonalConstraint(prog, R.row(1).transpose(), R.row(2).transpose());
  AddOrthogonalConstraint(prog, R.row(0).transpose(), R.row(2).transpose());
}

namespace {

// Decodes the discretization of the axes.
// For compactness, this method is referred to as phi(i) in the documentation
// below.  The implementation must give a valid number even for i<0 and
// i>num_binary_variables_per_half_axis.
double EnvelopeMinValue(int i, int num_binary_variables_per_half_axis) {
  return static_cast<double>(i) / num_binary_variables_per_half_axis;
}

// Given (an integer enumeration of) the orthant, takes a vector in the
// positive orthant into that orthant by flipping the signs of the individual
// elements.
Eigen::Vector3d FlipVector(const Eigen::Ref<const Eigen::Vector3d>& vpos,
                           int orthant) {
  DRAKE_ASSERT(vpos(0) >= 0 && vpos(1) >= 0 && vpos(2) >= 0);
  DRAKE_DEMAND(orthant >= 0 && orthant <= 7);
  Eigen::Vector3d v = vpos;
  if (orthant & (1 << 2)) v(0) = -v(0);
  if (orthant & (1 << 1)) v(1) = -v(1);
  if (orthant & 1) v(2) = -v(2);
  return v;
}

// Given (an integer enumeration of) the orthant, return a vector c with
// c(i) = a(i) if element i is positive in the indicated orthant, otherwise
// c(i) = b(i).
template <typename Derived>
Eigen::Matrix<Derived, 3, 1> PickPermutation(
    const Eigen::Matrix<Derived, 3, 1>& a,
    const Eigen::Matrix<Derived, 3, 1>& b, int orthant) {
  DRAKE_DEMAND(orthant >= 0 && orthant <= 7);
  Eigen::Matrix<Derived, 3, 1> c = a;
  if (orthant & (1 << 2)) c(0) = b(0);
  if (orthant & (1 << 1)) c(1) = b(1);
  if (orthant & 1) c(2) = b(2);
  return c;
}

// Given two coordinates, find the (positive) third coordinate on that
// intersects with the unit circle.
double Intercept(double x, double y) {
  DRAKE_ASSERT(x * x + y * y <= 1);
  return std::sqrt(1 - x * x - y * y);
}

}  // namespace

namespace internal {
// Given an axis-aligned box in the first orthant, computes and returns all the
// intersecting points between the edges of the box and the unit sphere.
// @param bmin  The vertex of the box closest to the origin.
// @param bmax  The vertex of the box farthest from the origin.
std::vector<Eigen::Vector3d> ComputeBoxEdgesAndSphereIntersection(
    const Eigen::Vector3d& bmin, const Eigen::Vector3d& bmax) {
  // Assumes the positive orthant (and bmax > bmin).
  DRAKE_ASSERT(bmin(0) >= 0 && bmin(1) >= 0 && bmin(2) >= 0);
  DRAKE_ASSERT(bmax(0) > bmin(0) && bmax(1) > bmin(1) && bmax(2) > bmin(2));

  // Assumes the unit circle intersects the box.
  DRAKE_ASSERT(bmin.lpNorm<2>() <= 1);
  DRAKE_ASSERT(bmax.lpNorm<2>() >= 1);

  std::vector<Eigen::Vector3d> intersections;

  if (bmin.lpNorm<2>() == 1) {
    // Then only the min corner intersects.
    intersections.push_back(bmin);
    return intersections;
  }

  if (bmax.lpNorm<2>() == 1) {
    // Then only the max corner intersects.
    intersections.push_back(bmax);
    return intersections;
  }

  // The box has at most 12 edges, each edge can intersect with the unit sphere
  // at most once, since the box is in the first orthant.
  intersections.reserve(12);

  // 1. Loop through each vertex of the box, add it to intersections if
  // the vertex is on the sphere.
  for (int i = 0; i < 8; ++i) {
    Eigen::Vector3d vertex{};
    for (int axis = 0; axis < 3; ++axis) {
      vertex(axis) = i & (1 << axis) ? bmin(axis) : bmax(axis);
    }
    if (vertex.norm() == 1) {
      intersections.push_back(vertex);
    }
  }

  // 2. Loop through each edge, find the intersection between each edge and the
  // unit sphere, if one exists.
  for (int axis = 0; axis < 3; ++axis) {
    // axis = 0 means edges along x axis;
    // axis = 1 means edges along y axis;
    // axis = 2 means edges along z axis;
    int fixed_axis1 = (axis + 1) % 3;
    int fixed_axis2 = (axis + 2) % 3;
    // 4 edges along each axis;

    // First finds the two end points on the edge.
    Eigen::Vector3d pt_closer, pt_farther;
    pt_closer(axis) = bmin(axis);
    pt_farther(axis) = bmax(axis);
    std::array<double, 2> fixed_axis1_val = {
        {bmin(fixed_axis1), bmax(fixed_axis1)}};
    std::array<double, 2> fixed_axis2_val = {
        {bmin(fixed_axis2), bmax(fixed_axis2)}};
    for (double val1 : fixed_axis1_val) {
      pt_closer(fixed_axis1) = val1;
      pt_farther(fixed_axis1) = pt_closer(fixed_axis1);
      for (double val2 : fixed_axis2_val) {
        pt_closer(fixed_axis2) = val2;
        pt_farther(fixed_axis2) = pt_closer(fixed_axis2);

        // Determines if there is an intersecting point between the edge and the
        // sphere. If the intersecting point is not the vertex of the box, then
        // push this intersecting point to intersections directly.
        if (pt_closer.norm() < 1 && pt_farther.norm() > 1) {
          Eigen::Vector3d pt_intersect{};
          pt_intersect(fixed_axis1) = pt_closer(fixed_axis1);
          pt_intersect(fixed_axis2) = pt_closer(fixed_axis2);
          pt_intersect(axis) =
              Intercept(pt_intersect(fixed_axis1), pt_intersect(fixed_axis2));
          intersections.push_back(pt_intersect);
        }
      }
    }
  }
  return intersections;
}

/**
 * Compute the outward unit length normal of the triangle, with the three
 * vertices being `pt0`, `pt1` and `pt2`.
 * @param pt0 A vertex of the triangle, in the first orthant (+++).
 * @param pt1 A vertex of the triangle, in the first orthant (+++).
 * @param pt2 A vertex of the triangle, in the first orthant (+++).
 * @param n The unit length normal vector of the triangle, pointing outward from
 * the origin.
 * @param d The intersecpt of the plane. Namely nᵀ * x = d for any point x on
 * the triangle.
 */
void ComputeTriangleOutwardNormal(const Eigen::Vector3d& pt0,
                                  const Eigen::Vector3d& pt1,
                                  const Eigen::Vector3d& pt2,
                                  Eigen::Vector3d* n, double* d) {
  DRAKE_DEMAND((pt0.array() >= 0).all());
  DRAKE_DEMAND((pt1.array() >= 0).all());
  DRAKE_DEMAND((pt2.array() >= 0).all());
  *n = (pt2 - pt0).cross(pt1 - pt0);
  // If the three points are almost colinear, then throw an error.
  double n_norm = n->norm();
  if (n_norm < 1E-3) {
    throw std::runtime_error("The points are almost colinear.");
  }
  *n = (*n) / n_norm;
  if (n->sum() < 0) {
    (*n) *= -1;
  }
  *d = pt0.dot(*n);
  DRAKE_DEMAND((n->array() >= 0).all());
}

/**
 * For the vertices in `pts`, determine if these vertices are co-planar. If they
 * are, then compute that plane nᵀ * x = d.
 * @param pts The vertices to be checked.
 * @param n The unit length normal vector of the plane, points outward from the
 * origin. If the vertices are not co-planar, leave `n` to 0.
 * @param d The intersecpt of the plane. If the vertices are not co-planar, set
 * this to 0.
 * @return If the vertices are co-planar, set this to true. Otherwise set to
 * false.
 */
bool AreAllVerticesCoPlanar(const std::vector<Eigen::Vector3d>& pts,
                            Eigen::Vector3d* n, double* d) {
  DRAKE_DEMAND(pts.size() >= 3);
  ComputeTriangleOutwardNormal(pts[0], pts[1], pts[2], n, d);
  // Determine if the other vertices are on the plane nᵀ * x = d.
  bool pts_on_plane = true;
  for (int i = 3; i < static_cast<int>(pts.size()); ++i) {
    if (std::abs(n->dot(pts[i]) - *d) > 1E-10) {
      pts_on_plane = false;
      n->setZero();
      *d = 0;
      break;
    }
  }
  return pts_on_plane;
}

/*
 * For the intersection region between the surface of the unit sphere, and the
 * interior of a box aligned with the axes, use a half space relaxation for
 * the intersection region as nᵀ * v >= d
 * @param[in] pts. The vertices containing the intersecting points between edges
 * of the box and the surface of the unit sphere.
 * @param[out] n. The unit length normal vector of the halfspace, pointing
 * outward.
 * @param[out] d. The intercept of the halfspace.
 */
void ComputeHalfSpaceRelaxationForBoxSphereIntersection(
    const std::vector<Eigen::Vector3d>& pts, Eigen::Vector3d* n, double* d) {
  DRAKE_DEMAND(pts.size() >= 3);
  // We first prove that for a given normal vector n, and ANY unit
  // length vector v within the intersection region between the
  // surface of the unit sphere and the interior of the axis-aligned
  // box, the minimal of nᵀ * v, always occurs at one of the vertex of
  // the intersection region, if the box and the vector n are in the
  // same orthant. Namely min nᵀ * v = min(nᵀ * pts.col(i))
  // To see this, for any vector v in the intersection region, suppose
  // it is on an arc, aligned with one axis. Without loss of
  // generality we assume the aligned axis is x axis, namely
  // v(0) = t, box_min(0) <= t <= box_max(0)
  // and v(1)² + v(2)² = 1 - t², with the bounds
  // box_min(1) <= v(1) <= box_max(1)
  // box_min(2) <= v(2) <= box_max(2)
  // And the inner product nᵀ * v =
  // n(0) * t + s * (n(1) * cos(α) + n(2) * sin(α))
  // where we define s = sqrt(1 - t²)
  // Using the property of trigonometric function, we know that
  // the minimal of (n(1) * cos(α) + n(2) * sin(α)) is obtained at
  // the boundary of α. Thus we know that the minimal of nᵀ * v is
  // always obtained at one of the vertex pts.col(i).

  // To find the tightest bound d satisfying nᵀ * v >= d for all
  // vector v in the intersection region, we use the fact that for
  // a given normal vector n, the minimal of nᵀ * v is always obtained
  // at one of the vertices pts.col(i), and formulate the following
  // SOCP to find the normal vector n
  // max d
  // s.t d <= nᵀ * pts.col(i)
  //     nᵀ * n <= 1

  // First compute a plane coinciding with a triangle, formed by 3 vertices
  // in the intersection region. If all the vertices are on that plane, then the
  // normal of the plane is n, and we do not need to run the optimization.
  // If there are only 3 vertices in the intersection region, then the normal
  // vector n is the normal of the triangle, formed by these three vertices.

  bool pts_on_plane = AreAllVerticesCoPlanar(pts, n, d);
  if (pts_on_plane) {
    return;
  }

  // If there are more than 3 vertices in the intersection region, and these
  // vertices are not co-planar, then we find the normal vector `n` through an
  // optimization, whose formulation is mentioned above.
  MathematicalProgram prog_normal;
  auto n_var = prog_normal.NewContinuousVariables<3>();
  auto d_var = prog_normal.NewContinuousVariables<1>();
  prog_normal.AddLinearCost(-d_var(0));
  for (const auto& pt : pts) {
    prog_normal.AddLinearConstraint(n_var.dot(pt) >= d_var(0));
  }

  // TODO(hongkai.dai): This optimization is expensive, especially if we have
  // multiple rotation matrices, all relaxed with the same number of binary
  // variables per half axis, the result `n` and `d` are the same. Should
  // consider hard-coding the result, to avoid repeated computation.

  Vector4<symbolic::Expression> lorentz_cone_vars;
  lorentz_cone_vars << 1, n_var;
  prog_normal.AddLorentzConeConstraint(lorentz_cone_vars);
  prog_normal.Solve();
  *n = prog_normal.GetSolution(n_var);
  *d = prog_normal.GetSolution(d_var(0));

  DRAKE_DEMAND((*n)(0) > 0 && (*n)(1) > 0 && (*n)(2) > 0);
  DRAKE_DEMAND(*d > 0 && *d < 1);
}

/**
 * For the intersection region between the surface of the unit sphere, and the
 * interior of a box aligned with the axes, relax this nonconvex intersection
 * region to its convex hull. This convex hull has some planar facets (formed
 * by the triangles connecting the vertices of the intersection region). This
 * function computes these planar facets. It is guaranteed that any point x on
 * the intersection region, satisfies A * x <= b.
 * @param[in] pts The vertices of the intersection region. Same as the `pts` in
 * ComputeHalfSpaceRelaxationForBoxSphereIntersection()
 * @param[out] A The rows of A are the normal vector of facets. Each row of A is
 * a unit length vector.
 * @param b b(i) is the interscept of the i'th facet.
 * @pre pts[i] are all in the first orthant, namely (pts[i].array() >=0).all()
 * should be true.
 */
void ComputeInnerFacetsForBoxSphereIntersection(
    const std::vector<Eigen::Vector3d>& pts,
    Eigen::Matrix<double, Eigen::Dynamic, 3>* A, Eigen::VectorXd* b) {
  for (const auto& pt : pts) {
    DRAKE_DEMAND((pt.array() >= 0).all());
  }
  A->resize(0, 3);
  b->resize(0);
  // Loop through each triangle, formed by connecting the vertices of the
  // intersection region. We write the plane coinciding with the triangle as
  // cᵀ * x >= d. If all the vertices of the intersection region satisfies
  // cᵀ * pts[i] >= d, then we know the intersection region satisfies
  // cᵀ * x >= d for all x being a point in the intersection region. Here we
  // use the proof in the ComputeHalfSpaceRelaxationForBoxSphereIntersection(),
  // that the minimal value of cᵀ * x over all x inside the intersection
  // region, occurs at one of the vertex of the intersection region.
  for (int i = 0; i < static_cast<int>(pts.size()); ++i) {
    for (int j = i + 1; j < static_cast<int>(pts.size()); ++j) {
      for (int k = j + 1; k < static_cast<int>(pts.size()); ++k) {
        // First compute the triangle formed by vertices pts[i], pts[j] and
        // pts[k].
        Eigen::Vector3d c;
        double d;
        ComputeTriangleOutwardNormal(pts[i], pts[j], pts[k], &c, &d);
        // A halfspace cᵀ * x >= d is valid, if all vertices pts[l] satisfy
        // cᵀ * pts[l] >= d.
        bool is_valid_halfspace = true;
        // Now check if the other vertices pts[l] satisfies cᵀ * pts[l] >= d.
        for (int l = 0; l < static_cast<int>(pts.size()); ++l) {
          if ((l != i) && (l != j) && (l != k)) {
            if (c.dot(pts[l]) < d - 1E-10) {
              is_valid_halfspace = false;
              break;
            }
          }
        }
        // If all vertices pts[l] satisfy cᵀ * pts[l] >= d, then add this
        // constraint to A * x <= b
        if (is_valid_halfspace) {
          A->conservativeResize(A->rows() + 1, Eigen::NoChange);
          b->conservativeResize(b->rows() + 1, Eigen::NoChange);
          A->row(A->rows() - 1) = -c.transpose();
          (*b)(b->rows() - 1) = -d;
        }
      }
    }
  }
}

/**
 * Add the McCormick constraint for a unit length vector `v`, if it lies within
 * a box box_min <= |v| <= box_max. where |v| is the elementwise absolute value
 * of v.
 * @param prog The optimization program.
 * @param v The unit length vector that will be constrained.
 * @param box_min The vertex closest to the origin.
 * @param box_max The vertex farthest from the origin.
 * @param this_cpos this_cpos(i) = 1 => box_min(i) <= v(i) <= box_max(i)
 * @param this_cneg this_cneg(i) = 1 => box_min(i) <= -v(i) <= box_max(i)
 * @param v1 Another vector that should be perpendicular to v
 * @param v2 Another vector that should be perpendicular to v, and v2 = cross(v, v1)
 */
void AddMcCormickVectorConstraintsWithBox(
    MathematicalProgram* prog, const VectorDecisionVariable<3>& v,
    const Eigen::Ref<const Eigen::Vector3d>& box_min,
    const Eigen::Ref<const Eigen::Vector3d>& box_max,
    const VectorDecisionVariable<3>& this_cpos,
    const VectorDecisionVariable<3>& this_cneg,
    const VectorDecisionVariable<3>& v1,
    const VectorDecisionVariable<3>& v2) {
  DRAKE_DEMAND((box_min.array() >= 0).all());
  DRAKE_DEMAND((box_max.array() >= box_min.array()).all());
  const double box_min_norm = box_min.lpNorm<2>();
  const double box_max_norm = box_max.lpNorm<2>();
  if (box_min_norm <= 1.0 + 2 * numeric_limits<double>::epsilon() &&
      box_max_norm >= 1.0 - 2 * numeric_limits<double>::epsilon()) {
    // The box intersects with the surface of the unit sphere.
    // Two possible cases
    // 1. If the box bmin <= x <= bmax intersects with the surface of the
    // unit sphere at a unique point (either bmin or bmax),
    // 2. Otherwise, there is a region of intersection.

    // We choose the error as 2 * eps here. The reason is that
    // if x.norm() == 1, then another vector y which is different from
    // x by eps (y(i) = x(i) + eps), the norm of y is at most 1 + 2 * eps.
    if (std::abs(box_min_norm - 1.0) <
        2 * numeric_limits<double>::epsilon() ||
        std::abs(box_max_norm - 1.0) <
            2 * numeric_limits<double>::epsilon()) {
      // If box_min or box_max is on the sphere, then denote the point on
      // the sphere as u, we have the following condition
      // if c[xi](0) = 1 and c[yi](1) == 1 and c[zi](2) == 1, then
      //     v = u
      //     vᵀ * v1 = 0
      //     vᵀ * v2 = 0
      //     v.cross(v1) = v2
      // Translate this to constraint, we have
      //   2 * (c[xi](0) + c[yi](1) + c[zi](2)) - 6
      //       <= v - u <= -2 * (c[xi](0) + c[yi](1) + c[zi](2)) + 6
      //
      //   c[xi](0) + c[yi](1) + c[zi](2) - 3
      //       <= vᵀ * v1 <= 3 - (c[xi](0) + c[yi](1) + c[zi](2))
      //
      //   c[xi](0) + c[yi](1) + c[zi](2) - 3
      //       <= vᵀ * v2 <= 3 - (c[xi](0) + c[yi](1) + c[zi](2))
      //
      //   2 * c[xi](0) + c[yi](1) + c[zi](2) - 6
      //       <= v.cross(v1) - v2 <= 6 - 2 * (c[xi](0) + c[yi](1) +
      //       c[zi](2))
      Eigen::Vector3d unique_intersection;  // `u` in the documentation
      // above
      if (std::abs(box_min_norm - 1.0) <
          2 * numeric_limits<double>::epsilon()) {
        unique_intersection = box_min / box_min_norm;
      } else {
        unique_intersection = box_max / box_max_norm;
      }
      Eigen::Vector3d orthant_u;
      VectorDecisionVariable<3> orthant_c;
      for (int o = 0; o < 8; o++) {  // iterate over orthants
        orthant_u = FlipVector(unique_intersection, o);
        orthant_c = PickPermutation(this_cpos, this_cneg, o);

        // TODO(hongkai.dai): remove this for loop when we can handle
        // Eigen::Array of symbolic formulae.
        Expression orthant_c_sum = orthant_c.cast<Expression>().sum();
        for (int i = 0; i < 3; ++i) {
          prog->AddLinearConstraint(v(i) - orthant_u(i) <=
              6 - 2 * orthant_c_sum);
          prog->AddLinearConstraint(v(i) - orthant_u(i) >=
              2 * orthant_c_sum - 6);
        }
        const Expression v_dot_v1 = orthant_u.dot(v1);
        const Expression v_dot_v2 = orthant_u.dot(v2);
        prog->AddLinearConstraint(v_dot_v1 <= 3 - orthant_c_sum);
        prog->AddLinearConstraint(orthant_c_sum - 3 <= v_dot_v1);
        prog->AddLinearConstraint(v_dot_v2 <= 3 - orthant_c_sum);
        prog->AddLinearConstraint(orthant_c_sum - 3 <= v_dot_v2);
        const Vector3<Expression> v_cross_v1 = orthant_u.cross(v1);
        for (int i = 0; i < 3; ++i) {
          prog->AddLinearConstraint(v_cross_v1(i) - v2(i) <=
              6 - 2 * orthant_c_sum);
          prog->AddLinearConstraint(v_cross_v1(i) - v2(i) >=
              2 * orthant_c_sum - 6);
        }
      }
    } else {
      // Find the intercepts of the unit sphere with the box, then find
      // the tightest linear constraint of the form:
      //    d <= n'*v
      // that puts v inside (but as close as possible to) the unit circle.
      auto pts = internal::ComputeBoxEdgesAndSphereIntersection(box_min,
                                                                box_max);
      DRAKE_DEMAND(pts.size() >= 3);

      double d(0);
      Eigen::Vector3d normal{};
      internal::ComputeHalfSpaceRelaxationForBoxSphereIntersection(
          pts, &normal, &d);

      Eigen::VectorXd b(0);
      Eigen::Matrix<double, Eigen::Dynamic, 3> A(0, 3);

      internal::ComputeInnerFacetsForBoxSphereIntersection(pts, &A, &b);

      // theta is the maximal angle between v and normal, where v is an
      // intersecting point between the box and the sphere.
      double cos_theta = d;
      const double theta = std::acos(cos_theta);

      Eigen::Matrix<double, 1, 6> a;
      Eigen::Matrix<double, 3, 9> A_cross;

      Eigen::RowVector3d orthant_normal;
      VectorDecisionVariable<3> orthant_c;
      for (int o = 0; o < 8; o++) {  // iterate over orthants
        orthant_normal = FlipVector(normal, o).transpose();
        orthant_c = PickPermutation(this_cpos, this_cneg, o);

        for (int i = 0; i < A.rows(); ++i) {
          // Add the constraint that A * v <= b, representing the inner
          // facets f the convex hull, obtained from the vertices of the
          // intersection region.
          // This constraint is only active if the box is active.
          // We impose the constraint
          // A.row(i) * v - b(i) <= 3 - 3 * b(i) + (b(i) - 1) * (c[xi](0)
          // + c[yi](1) + c[zi](2))
          // Or in words
          // If c[xi](0) = 1 and c[yi](1) = 1 and c[zi](1) = 1
          //   A.row(i) * v <= b(i)
          // Otherwise
          //   A.row(i) * v -b(i) is not constrained
          Eigen::Vector3d orthant_a =
              -FlipVector(-A.row(i).transpose(), o);
          prog->AddLinearConstraint(
              orthant_a.dot(v) - b(i) <=
                  3 - 3 * b(i) +
                      (b(i) - 1) *
                          orthant_c.cast<symbolic::Expression>().sum());
        }

        // Max vector norm constraint: -1 <= normal'*x <= 1.
        // No need to restrict to this orthant, but also no need to apply
        // the same constraint twice (would be the same for opposite
        // orthants), so skip all of the -x orthants.
        if (o % 2 == 0)
          prog->AddLinearConstraint(orthant_normal, -1, 1, v);

        // Dot-product constraint: ideally v.dot(v1) = v.dot(v2) = 0.
        // The cone of (unit) vectors within theta of the normal vector
        // defines a band of admissible vectors v1 and v2 which are
        // orthogonal to v.  They must satisfy the constraint:
        //    -sin(theta) <= normal.dot(vi) <= sin(theta)
        // Proof sketch:
        //   v is within theta of normal.
        //   => vi must be within theta of a vector orthogonal to
        //      the normal.
        //   => vi must be pi/2 +/- theta from the normal.
        //   => |normal||vi| cos(pi/2 + theta) <= normal.dot(vi) <=
        //                |normal||vi| cos(pi/2 - theta).
        // Since normal and vi are both unit length,
        //     -sin(theta) <= normal.dot(vi) <= sin(theta).
        // To activate this only when this box is active, we use
        //   -sin(theta)-6+2*c[xi](0)+2*c[yi](1)+2*c[zi](2) <=
        //   normal.dot(vi)
        //     normal.dot(vi) <=
        //     sin(theta)+6-2*c[xi](0)-2*c[yi](1)-2*c[zi](2).
        // Note: (An alternative tighter, but SOCP constraint)
        //   v, v1, v2 forms an orthornormal basis. So n'*v is the
        //   projection of n in the v direction, same for n'*v1, n'*v2.
        //   Thus
        //     (n'*v)^2 + (n'*v1)^2 + (n'*v2)^2 = n'*n
        //   which translates to "The norm of a vector is equal to the
        //   sum of squares of the vector projected onto each axes of an
        //   orthornormal basis".
        //   This equation is the same as
        //     (n'*v1)^2 + (n'*v2)^2 <= sin(theta)^2
        //   So actually instead of imposing
        //     -sin(theta)<=n'*vi <=sin(theta),
        //   we can impose a tighter Lorentz cone constraint
        //     [|sin(theta)|, n'*v1, n'*v2] is in the Lorentz cone.
        a << orthant_normal, -2, -2, -2;
        prog->AddLinearConstraint(a, -sin(theta) - 6, 1, {v1, orthant_c});
        prog->AddLinearConstraint(a, -sin(theta) - 6, 1, {v2, orthant_c});

        a.tail<3>() << 2, 2, 2;
        prog->AddLinearConstraint(a, -1, sin(theta) + 6, {v1, orthant_c});
        prog->AddLinearConstraint(a, -1, sin(theta) + 6, {v2, orthant_c});

        // Cross-product constraint: ideally v2 = v.cross(v1).
        // Since v is within theta of normal, we will prove that
        // |v2 - normal.cross(v1)| <= 2 * sin(theta / 2)
        // Notice that (v2 - normal.cross(v1))ᵀ * (v2 - normal.cross(v1))
        // = v2ᵀ * v2 + (normal.cross(v1))ᵀ * (normal.cross(v1)) -
        //      2 * v2ᵀ * (normal.cross(v1))
        // <= 1 + 1 - 2 * cos(theta)
        // = (2 * sin(theta / 2))²
        // Thus we get |v2 - normal.cross(v1)| <= 2 * sin(theta / 2)
        // Here we consider to use an elementwise linear constraint
        // -2*sin(theta / 2) <=  v2 - normal.cross(v1) <= 2*sin(theta / 2)
        // Since 0<=theta<=pi/2, this should be enough to rule out the
        // det(R)=-1 case (the shortest projection of a line across the
        // circle onto a single axis has length 2sqrt(3)/3 > 1.15), and
        // can be significantly tighter.

        // To activate this only when the box is active, the complete
        // constraints are
        //  -2*sin(theta/2)-6+2(cxi+cyi+czi) <= v2-normal.cross(v1)
        //    v2-normal.cross(v1) <= 2*sin(theta/2)+6-2(cxi+cyi+czi)
        // Note: Again this constraint could be tighter as a Lorenz cone
        // constraint of the form:
        //   |v2 - normal.cross(v1)| <= 2*sin(theta/2).
        A_cross << Eigen::Matrix3d::Identity(),
            -math::VectorToSkewSymmetric(orthant_normal),
            Eigen::Matrix3d::Constant(-2);
        prog->AddLinearConstraint(
            A_cross, Eigen::Vector3d::Constant(-2 * sin(theta / 2) - 6),
            Eigen::Vector3d::Constant(2), {v2, v1, orthant_c});

        A_cross.rightCols<3>() = Eigen::Matrix3d::Constant(2.0);
        prog->AddLinearConstraint(
            A_cross, Eigen::Vector3d::Constant(-2),
            Eigen::Vector3d::Constant(2 * sin(theta / 2) + 6),
            {v2, v1, orthant_c});
      }
    }
  } else {
    // This box does not intersect with the surface of the sphere.
    for (int o = 0; o < 8; ++o) {  // iterate over orthants
      prog->AddLinearConstraint(PickPermutation(this_cpos, this_cneg, o)
                                    .cast<Expression>()
                                    .sum(),
                                0.0, 2.0);
    }
  }
}
}  // namespace internal

namespace {

void AddMcCormickVectorConstraints(
    MathematicalProgram* prog, const VectorDecisionVariable<3>& v,
    const std::vector<MatrixDecisionVariable<3, 1>>& cpos,
    const std::vector<MatrixDecisionVariable<3, 1>>& cneg,
    const std::vector<VectorDecisionVariable<3>>& bpos,
    const std::vector<VectorDecisionVariable<3>>& bneg,
    const VectorDecisionVariable<3>& v1, const VectorDecisionVariable<3>& v2) {
  const int N = cpos.size();  // number of discretization points.

  // Iterate through regions.
  Eigen::Vector3d box_min, box_max;
  for (int xi = 0; xi < N; xi++) {
    box_min(0) = EnvelopeMinValue(xi, N);
    box_max(0) = EnvelopeMinValue(xi + 1, N);
    for (int yi = 0; yi < N; yi++) {
      box_min(1) = EnvelopeMinValue(yi, N);
      box_max(1) = EnvelopeMinValue(yi + 1, N);
      for (int zi = 0; zi < N; zi++) {
        box_min(2) = EnvelopeMinValue(zi, N);
        box_max(2) = EnvelopeMinValue(zi + 1, N);

        VectorDecisionVariable<3> this_cpos, this_cneg;
        this_cpos << cpos[xi](0), cpos[yi](1), cpos[zi](2);
        this_cneg << cneg[xi](0), cneg[yi](1), cneg[zi](2);

        VectorDecisionVariable<3> this_bpos, this_bneg;
        this_bpos << bpos[xi](0), bpos[yi](1), bpos[zi](2);
        this_bneg << bneg[xi](0), bneg[yi](1), bneg[zi](2);

        internal::AddMcCormickVectorConstraintsWithBox(prog, v, box_min, box_max, this_cpos, this_cneg, v1, v2);
        if (xi == 0 && yi == 0 && zi == 0) {
          // If we know which orthant the point `v` is in, then we can constrain
          // `v` to be in the McCormick envelope of that orthant surface, this
          // can potentially rule out the R = 0 case.
          internal::AddMcCormickVectorConstraintsWithBox(prog,
                                                         v,
                                                         box_min,
                                                         Eigen::Vector3d::Ones(),
                                                         this_bpos,
                                                         this_bneg,
                                                         v1,
                                                         v2);
        }
      }
    }
  }
}

/**
 * Add the constraint that vector R.col(i) and R.col(j) are not in the
 * same or opposite orthants. This constraint should be satisfied since
 * R.col(i) should be perpendicular to R.col(j). For example, if both
 * R.col(i) and R.col(j) are in the first orthant (+++), their inner product
 * has to be non-negative. If the inner product of two first orthant vectors
 * is exactly zero, then both vectors has to be on the boundaries of the first
 * orthant. But we can then assign the vector to a different orthant. The same
 * proof applies to the opposite orthant case.
 * To impose the constraint that R.col(0) and R.col(1) are not both in the first
 * orthant, we consider the constraint
 * Bpos0.col(0).sum() + Bpos0.col(1).sum() <= 5.
 * Namely, not all 6 entries in Bpos0.col(0) and Bpos0.col(1) can be 1 at the
 * same time, which is another way of saying R.col(0) and R.col(1) cannot be
 * both in the first orthant.
 * Similarly we can impose the constraint on the other orthant.
 * @param prog Add the constraint to this mathematical program.
 * @param Bpos0 Defined in AddRotationMatrixMcCormickEnvelopeMilpConstraints(),
 * Bpos0(i,j) = 1 => R(i, j) >= 0.
 * @param Bneg0 Defined in AddRotationMatrixMcCormickEnvelopeMilpConstraints(),
 * Bneg0(i,j) = 1 => R(i, j) <= 0.
 */
void AddNotInSameOrOppositeOrthantConstraint(
    MathematicalProgram* prog,
    const Eigen::Ref<const MatrixDecisionVariable<3, 3>>& Bpos0,
    const Eigen::Ref<const MatrixDecisionVariable<3, 3>>& Bneg0) {
  const std::array<std::pair<int, int>, 3> column_idx = {
      {{0, 1}, {0, 2}, {1, 2}}};
  for (const auto& column_pair : column_idx) {
    const int col_idx0 = column_pair.first;
    const int col_idx1 = column_pair.second;
    for (int o = 0; o < 8; ++o) {
      // To enforce that R.col(i) and R.col(j) are not simultaneously in the
      // o'th orthant, we will impose the constraint
      // vars_same_orthant.sum() < = 5. The variables in vars_same_orthant
      // depend on the orthant number o.
      // To enforce that R.col(i) and R.col(j) are not in the opposite
      // orthants, we will impose the constraint
      // vars_oppo_orthant.sum() <= 5. The variables in vars_oppo_orthant
      // depnd on the orthant number o.
      Eigen::Matrix<symbolic::Variable, 6, 1> vars_same_orthant;
      Eigen::Matrix<symbolic::Variable, 6, 1> vars_oppo_orthant;
      for (int axis = 0; axis < 3; ++axis) {
        // axis chooses x, y, or z axis.
        if (o & (1 << axis)) {
          // If the orthant has positive value along the `axis`, then
          // `vars_same_orthant` choose the positive component Bpos0.
          vars_same_orthant(2 * axis)     = Bpos0(axis, col_idx0);
          vars_same_orthant(2 * axis + 1) = Bpos0(axis, col_idx1);
          vars_oppo_orthant(2 * axis)     = Bpos0(axis, col_idx0);
          vars_oppo_orthant(2 * axis + 1) = Bneg0(axis, col_idx1);
        } else {
          // If the orthant has negative value along the `axis`, then
          // `vars_same_orthant` choose the negative component Bneg0.
          vars_same_orthant(2 * axis)     = Bneg0(axis, col_idx0);
          vars_same_orthant(2 * axis + 1) = Bneg0(axis, col_idx1);
          vars_oppo_orthant(2 * axis)     = Bneg0(axis, col_idx0);
          vars_oppo_orthant(2 * axis + 1) = Bpos0(axis, col_idx1);
        }
      }
      prog->AddLinearConstraint(Eigen::Matrix<double, 1, 6>::Ones(), 0, 5,
                                vars_same_orthant);
      prog->AddLinearConstraint(Eigen::Matrix<double, 1, 6>::Ones(), 0, 5,
                                vars_oppo_orthant);
    }
  }
}
}  // namespace

AddRotationMatrixMcCormickEnvelopeReturnType
AddRotationMatrixMcCormickEnvelopeMilpConstraints(
    MathematicalProgram* prog,
    const Eigen::Ref<const MatrixDecisionVariable<3, 3>>& R,
    int num_binary_vars_per_half_axis, RollPitchYawLimits limits) {
  DRAKE_DEMAND(num_binary_vars_per_half_axis >= 1);

  // Use a simple lambda to make the constraints more readable below.
  // Note that
  //  forall k>=0, 0<=phi(k), and
  //  forall k<=num_binary_vars_per_half_axis, phi(k)<=1.
  auto phi = [&](int k) -> double {
    return EnvelopeMinValue(k, num_binary_vars_per_half_axis);
  };

  // Creates binary decision variables which discretize each axis.
  //   BRpos[k](i,j) = 1 => R(i,j) >= phi(k)
  //   BRneg[k](i,j) = 1 => R(i,j) <= -phi(k)
  //
  // For convenience, we introduce additional (continuous) variables to
  // represent the individual sections of the real line
  //   CRpos[k](i,j) = BRpos[k](i,j) if k=N-1, otherwise
  //   CRpos[k](i,j) = BRpos[k](i,j) - BRpos[k+1](i,j)
  // This is useful only because the *number of decision variables* that we
  // pass into the constraints changes for the k=N-1 case.  Otherwise we
  // could do a simple substitution everywhere.
  // TODO(russt): Use symbolic constraints and remove these decision variables!
  std::vector<MatrixDecisionVariable<3, 3>> BRpos, BRneg;
  std::vector<MatrixDecisionVariable<3, 3>> CRpos, CRneg;
  for (int k = 0; k < num_binary_vars_per_half_axis; k++) {
    BRpos.push_back(
        prog->NewBinaryVariables<3, 3>("BRpos" + std::to_string(k)));
    BRneg.push_back(
        prog->NewBinaryVariables<3, 3>("BRneg" + std::to_string(k)));
    CRpos.push_back(
        prog->NewContinuousVariables<3, 3>("CRpos" + std::to_string(k)));
    CRneg.push_back(
        prog->NewContinuousVariables<3, 3>("CRneg" + std::to_string(k)));
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < num_binary_vars_per_half_axis; k++) {
        // R(i,j) > phi(k) => BRpos[k](i,j) = 1
        // R(i,j) < phi(k) => BRpos[k](i,j) = 0
        // R(i,j) = phi(k) => BRpos[k](i,j) = 0 or 1
        // Since -s1 <= R(i, j) - phi(k) <= s2,
        // where s1 = 1 + phi(k), s2 = 1 - phi(k). The point
        // [R(i,j) - phi(k), BRpos[k](i,j)] has to lie within the convex hull,
        // whose vertices are (-s1, 0), (0, 0), (s2, 1), (0, 1). By computing
        // the edges of this convex hull, we get
        // -s1 + s1*BRpos[k](i,j) <= R(i,j)-phi(k) <= s2 * BRpos[k](i,j)
        double s1 = 1 + phi(k);
        double s2 = 1 - phi(k);
        prog->AddLinearConstraint(R(i, j) - phi(k) >=
                                  -s1 + s1 * BRpos[k](i, j));
        prog->AddLinearConstraint(R(i, j) - phi(k) <= s2 * BRpos[k](i, j));

        // -R(i,j) > phi(k) => BRneg[k](i,j) = 1
        // -R(i,j) < phi(k) => BRneg[k](i,j) = 0
        // -R(i,j) = phi(k) => BRneg[k](i,j) = 0 or 1
        // Since -s2 <= R(i, j) + phi(k) <= s1,
        // where s1 = 1 + phi(k), s2 = 1 - phi(k). The point
        // [R(i,j) + phi(k), BRneg[k](i,j)] has to lie within the convex hull
        // whose vertices are (-s2, 1), (0, 0), (s1, 0), (0, 1). By computing
        // the edges of the convex hull, we get
        // -s2 * BRneg[k](i,j) <= R(i,j)+phi(k) <= s1-s1*BRneg[k](i,j)
        prog->AddLinearConstraint(R(i, j) + phi(k) <= s1 - s1 * BRneg[k](i, j));
        prog->AddLinearConstraint(R(i, j) + phi(k) >= -s2 * BRneg[k](i, j));

        if (k == num_binary_vars_per_half_axis - 1) {
          //   CRpos[k](i,j) = BRpos[k](i,j)
          prog->AddLinearConstraint(CRpos[k](i, j) == BRpos[k](i, j));

          //   CRneg[k](i,j) = BRneg[k](i,j)
          prog->AddLinearConstraint(CRneg[k](i, j) == BRneg[k](i, j));
        } else {
          //   CRpos[k](i,j) = BRpos[k](i,j) - BRpos[k+1](i,j)
          prog->AddLinearConstraint(CRpos[k](i, j) ==
                                    BRpos[k](i, j) - BRpos[k + 1](i, j));
          //   CRneg[k](i,j) = BRneg[k](i,j) - BRneg[k+1](i,j)
          prog->AddLinearConstraint(CRneg[k](i, j) ==
                                    BRneg[k](i, j) - BRneg[k + 1](i, j));
        }
      }
      // R(i,j) has to pick a side, either non-positive or non-negative.
      prog->AddLinearConstraint(BRpos[0](i, j) + BRneg[0](i, j) == 1); 
      // Any indicator variable of one sign should force all
      // indicator variables of the opposite sign to be off.
      for (int k = 0; k < num_binary_vars_per_half_axis; k++) {
        for (int k_other = 0; k_other < num_binary_vars_per_half_axis; k_other++) {
          if (k == k_other && k == 0) continue; // Skip the case constrained above.
          prog->AddLinearConstraint(BRpos[k](i, j) + BRneg[k](i, j) <= 1);
        }
      }
      // for debugging: constrain to positive orthant.
      //      prog->AddBoundingBoxConstraint(1,1,{BRpos[0].block<1,1>(i,j)});
    }
  }

  // Add constraint that no two rows (or two columns) can lie in the same
  // orthant (or opposite orthant).
  AddNotInSameOrOppositeOrthantConstraint(prog, BRpos[0], BRneg[0]);
  AddNotInSameOrOppositeOrthantConstraint(prog, BRpos[0].transpose(),
                                          BRneg[0].transpose());

  // Add angle limit constraints.
  // Bounding box will turn on/off an orthant.  It's sufficient to add the
  // constraints only to the positive orthant.
  AddBoundingBoxConstraintsImpliedByRollPitchYawLimitsToBinary(prog, BRpos[0],
                                                               limits);

  // Add constraints to the column and row vectors.
  std::vector<VectorDecisionVariable<3>> cpos(num_binary_vars_per_half_axis),
      cneg(num_binary_vars_per_half_axis);
  std::vector<VectorDecisionVariable<3>> bpos(num_binary_vars_per_half_axis),
      bneg(num_binary_vars_per_half_axis);
  for (int i = 0; i < 3; i++) {
    // Make lists of the decision variables in terms of column vectors and row
    // vectors to facilitate the calls below.
    // TODO(russt): Consider reorganizing the original CRpos/CRneg variables to
    // avoid this (albeit minor) cost?
    for (int k = 0; k < num_binary_vars_per_half_axis; k++) {
      cpos[k] = CRpos[k].col(i);
      cneg[k] = CRneg[k].col(i);
      bpos[k] = BRpos[k].col(i);
      bneg[k] = BRneg[k].col(i);
    }
    AddMcCormickVectorConstraints(prog, R.col(i), cpos, cneg, bpos, bneg,
                                  R.col((i + 1) % 3), R.col((i + 2) % 3));

    for (int k = 0; k < num_binary_vars_per_half_axis; k++) {
      cpos[k] = CRpos[k].row(i).transpose();
      cneg[k] = CRneg[k].row(i).transpose();
      bpos[k] = BRpos[k].row(i).transpose();
      bneg[k] = BRneg[k].row(i).transpose();
    }
    AddMcCormickVectorConstraints(prog, R.row(i).transpose(), cpos, cneg, bpos, bneg,
                                  R.row((i + 1) % 3).transpose(),
                                  R.row((i + 2) % 3).transpose());
  }
  return make_tuple(CRpos, CRneg, BRpos, BRneg);
}

symbolic::Variable GetOrConstructBilinearApproxTerm(
    MathematicalProgram* prog,
    std::map<std::pair<size_t, size_t>, symbolic::Variable> bilinear_approx_terms,
    symbolic::Variable x, symbolic::Variable y,
    int num_bins_per_bilinear_term){

  std::pair<size_t, size_t> bilinear_term_ids(x.get_hash(), y.get_hash());
  auto iter = bilinear_approx_terms.find(bilinear_term_ids);
  if (iter == bilinear_approx_terms.end()){
    // If we didn't find this term already, go make a new approximation of
    // it, and store it in our bilinear term storage.
    auto w = prog->NewContinuousVariables<1>("w_" + x.get_name() + "_" + y.get_name());
    prog->AddBoundingBoxConstraint(-1., 1., w);
    Add2DLogarithmicMcCormickEnvelope(prog, w[0], x, y, w[0].get_name(), -1., 1., -1., 1.,
    num_bins_per_bilinear_term, num_bins_per_bilinear_term);
    bilinear_approx_terms[bilinear_term_ids] = w[0];
    iter = bilinear_approx_terms.find(bilinear_term_ids);
    DRAKE_ASSERT(iter != bilinear_approx_terms.end());
  }
  return iter->second;
}

void
AddRotationMatrix2DLogMcCormickEnvelopeMilpConstraints(
    MathematicalProgram* prog,
    const Eigen::Ref<const MatrixDecisionVariable<3, 3>>& R,
    int num_bins_per_bilinear_term) {
  DRAKE_DEMAND(num_bins_per_bilinear_term >= 1);

  // Keep track of bilinear approximations we make with a dictionary,
  // so we don't repeat
  std::map<
    std::pair<size_t, size_t>,
    symbolic::Variable
  > bilinear_approx_terms;

  // Add constraints that fall out of R^T R = I
  for (int i = 0; i < 3; i++) {
    for (int j = i; j < 3; j ++) {
      VectorDecisionVariable<3> r_a = R.col(i);
      VectorDecisionVariable<3> r_b = R.col(j);
      double equals = i == j ? 1. : 0.;

      Expression sum_of_bilinear_terms = 0.0;
      for (int k = 0; k < 3; k ++) {
        sum_of_bilinear_terms += GetOrConstructBilinearApproxTerm(prog, bilinear_approx_terms, r_a[k], r_b[k], num_bins_per_bilinear_term);
      }
      prog->AddLinearConstraint(sum_of_bilinear_terms == equals);
    }
  }

  // And add the additional constraint that r_1 cross r_2 = r_
  Eigen::Matrix3i skew_sym_index_and_sign;
  skew_sym_index_and_sign << 0, -3, 2, 3, 0, -1, -2, 1, 0; // Indices are 1-indexed.

  for (int i = 0; i < 3; i++) { // i^{th} row of skew sym matrix times r_2 = r_{3, i}
    Expression lhs = 0.0;
    for (int j = 0; j < 3; j++) {
      if (i != j) { // diagonal terms are all 0
        int ind = skew_sym_index_and_sign(i, j);
        DRAKE_ASSERT(ind != 0);
        auto new_expr = GetOrConstructBilinearApproxTerm(prog, bilinear_approx_terms, 
            R(abs(ind)-1, 0), 
            R(j, 1), num_bins_per_bilinear_term);
        std::cout << "Adding " << (ind / abs(ind)) * new_expr << std::endl;
        lhs += (ind / abs(ind)) * new_expr;
        std::cout << "Added!" << std::endl;
      }
    }
    std::cout << "Cross prod constraint says " << lhs << " == " << R(i, 2) << std::endl;
    prog->AddLinearConstraint(lhs == R(i, 2));
  }
  
}

}  // namespace solvers
}  // namespace drake
