#ifndef CONSTR_VIOL_H
#define CONSTR_VIOL_H

/**
  * Violations of (mainly functional) constraints
  */

#include <cmath>

#include "constr_functional.h"
#include "constr_general.h"

namespace mp {

/// Compute result of the max constraint.
template <class VarVec>
double ComputeValue(const MaxConstraint& con, const VarVec& x) {
  double result = -INFINITY;
  for (auto i: con.GetArguments()) {
    if (result < x[i])
      result = x[i];
  }
  return result;
}

/// Compute result of the min constraint.
template <class VarVec>
double ComputeValue(const MinConstraint& con, const VarVec& x) {
  double result = INFINITY;
  for (auto i: con.GetArguments()) {
    if (result > x[i])
      result = x[i];
  }
  return result;
}

/// Compute result of the abs constraint.
template <class VarVec>
double ComputeValue(const AbsConstraint& con, const VarVec& x) {
  return std::fabs(x[con.GetArguments()[0]]);
}

/// Compute result of the and constraint.
template <class VarVec>
double ComputeValue(const AndConstraint& con, const VarVec& x) {
  for (auto i: con.GetArguments()) {
    if (x[i] < 0.5)
      return 0.0;
  }
  return 1.0;
}

/// Compute result of the or constraint.
template <class VarVec>
double ComputeValue(const OrConstraint& con, const VarVec& x) {
  for (auto i: con.GetArguments()) {
    if (x[i] >= 0.5)
      return 1.0;
  }
  return 0.0;
}

/// Compute result of the not constraint.
template <class VarVec>
double ComputeValue(const NotConstraint& con, const VarVec& x) {
  return x[con.GetArguments()[0]] < 0.5;
}


}  // namespace mp

#endif // CONSTR_VIOL_H
