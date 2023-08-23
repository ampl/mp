#ifndef CONSTR_EVAL_H
#define CONSTR_EVAL_H

/**
  * Evaluations and violations
  * of (mainly functional) constraints
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

/// Compute result of the div constraint.
template <class VarVec>
double ComputeValue(const DivConstraint& con, const VarVec& x) {
  auto divt = x[con.GetArguments()[0]];
  auto divs = x[con.GetArguments()[1]];
  return 0.0==std::fabs(divs)
      ? (divt>=0.0 ? INFINITY : -INFINITY)
      : divt / divs;
}

/// Compute result of the IfThen constraint.
template <class VarVec>
double ComputeValue(const IfThenConstraint& con, const VarVec& x) {
  auto i = x[con.GetArguments()[0]];
  return x[
      con.GetArguments()[(i>=0.5) ? 1 : 2]];
}

/// Compute result of the Implication constraint.
template <class VarVec>
double ComputeValue(const ImplicationConstraint& con, const VarVec& x) {
  auto i = x[con.GetArguments()[0]];
  auto c1 = x[con.GetArguments()[1]];
  auto c2 = x[con.GetArguments()[2]];
  return (i>=0.5 && c1>=0.5) || (i<0.5 && c2>=0.5);
}

/// Compute result of the AllDiff constraint.
template <class VarVec>
double ComputeValue(const AllDiffConstraint& con, const VarVec& x) {
  const auto& args = con.GetArguments();
  for (auto i=args.size(); i--; ) {
    for (auto j=i; j--; ) {         // Should be integer vars.
      if (std::round(x[args[i]]) == std::round(x[args[j]]))
        return 0.0;
    }
  }
  return 1.0;
}

/// Compute result of the NumberofConst constraint.
/// Currently assumes integer variables
/// (otherwise the result may be larger.)
template <class VarVec>
double ComputeValue(const NumberofConstConstraint& con, const VarVec& x) {
  double result = 0.0;
  auto k = std::round(con.GetParameters()[0]);
  for (auto v: con.GetArguments()) {
    if (std::round(x[v]) == k)
      ++result;
  }
  return result;
}

/// Compute result of the NumberofVar constraint.
/// Currently assumes integer variables
/// (otherwise the result may be larger.)
template <class VarVec>
double ComputeValue(const NumberofVarConstraint& con, const VarVec& x) {
  double result = 0.0;
  const auto& args = con.GetArguments();
  auto k = std::round(x[args[0]]);
  for (auto i=args.size(); --i; ) {
    if (std::round(x[args[i]]) == k)
      ++result;
  }
  return result;
}

/// Compute result of the count constraint.
template <class VarVec>
double ComputeValue(const CountConstraint& con, const VarVec& x) {
  double result = 0.0;
  for (auto v: con.GetArguments()) {
    if (x[v] >= 0.5)
      ++result;
  }
  return result;
}

/// Compute result of the exp constraint.
template <class VarVec>
double ComputeValue(const ExpConstraint& con, const VarVec& x) {
  return std::exp(x[con.GetArguments()[0]]);
}

/// Compute result of the expA constraint.
template <class VarVec>
double ComputeValue(const ExpAConstraint& con, const VarVec& x) {
  return std::pow(
        con.GetParameters()[0], x[con.GetArguments()[0]]);
}

/// Compute result of the log constraint.
template <class VarVec>
double ComputeValue(const LogConstraint& con, const VarVec& x) {
  return std::log(x[con.GetArguments()[0]]);
}

/// Compute result of the logA constraint.
template <class VarVec>
double ComputeValue(const LogAConstraint& con, const VarVec& x) {
  return
      std::log(x[con.GetArguments()[0]])
      / std::log(con.GetParameters()[0]);
}

/// Compute result of the pow constraint.
template <class VarVec>
double ComputeValue(const PowConstraint& con, const VarVec& x) {
  return std::pow(
        x[con.GetArguments()[0]], con.GetParameters()[0]);
}

/// Compute result of the sin constraint.
template <class VarVec>
double ComputeValue(const SinConstraint& con, const VarVec& x) {
  return std::sin(x[con.GetArguments()[0]]);
}

/// Compute result of the cos constraint.
template <class VarVec>
double ComputeValue(const CosConstraint& con, const VarVec& x) {
  return std::cos(x[con.GetArguments()[0]]);
}

/// Compute result of the tan constraint.
template <class VarVec>
double ComputeValue(const TanConstraint& con, const VarVec& x) {
  return std::tan(x[con.GetArguments()[0]]);
}

/// Compute result of the asin constraint.
template <class VarVec>
double ComputeValue(const AsinConstraint& con, const VarVec& x) {
  return std::asin(x[con.GetArguments()[0]]);
}

/// Compute result of the acos constraint.
template <class VarVec>
double ComputeValue(const AcosConstraint& con, const VarVec& x) {
  return std::acos(x[con.GetArguments()[0]]);
}

/// Compute result of the atan constraint.
template <class VarVec>
double ComputeValue(const AtanConstraint& con, const VarVec& x) {
  return std::atan(x[con.GetArguments()[0]]);
}

/// Compute result of the sinh constraint.
template <class VarVec>
double ComputeValue(const SinhConstraint& con, const VarVec& x) {
  return std::sinh(x[con.GetArguments()[0]]);
}

/// Compute result of the cosh constraint.
template <class VarVec>
double ComputeValue(const CoshConstraint& con, const VarVec& x) {
  return std::cosh(x[con.GetArguments()[0]]);
}

/// Compute result of the tanh constraint.
template <class VarVec>
double ComputeValue(const TanhConstraint& con, const VarVec& x) {
  return std::tanh(x[con.GetArguments()[0]]);
}

/// Compute result of a conditional constraint.
/// Actually, for violation itself, we could do this:
/// If the subconstr is violated but should hold,
/// return the exact gap, despite this is a logical constraint.
/// If the subconstr holds but should not, return 0.0.
template <class Con, class VarVec>
double ComputeValue(
    const ConditionalConstraint<Con>& con, const VarVec& x) {
  auto viol = con.GetConstraint().ComputeViolation(x);
  bool ccon_valid = viol<=0.0;
  bool has_arg = x[con.GetResultVar()] >= 0.5;
  switch (con.GetContext()) {
  case Context::CTX_MIX:
    return has_arg == ccon_valid;
  case Context::CTX_POS:
    return has_arg < ccon_valid;
  case Context::CTX_NEG:
    return has_arg > ccon_valid;
  default:
    return INFINITY;
  }
}


}  // namespace mp

#endif // CONSTR_EVAL_H
