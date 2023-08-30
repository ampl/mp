#ifndef CONSTR_EVAL_H
#define CONSTR_EVAL_H

/**
  * Evaluations and violations
  * of (mainly functional) constraints.
  *
  * For most evaluators, it's enough to supply
  * a simple vector x.
  * For some, it needs to be an object with
  * extra API.
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
template <class VarInfo>
double ComputeValue(const NumberofConstConstraint& con, const VarInfo& x) {
  double result = 0.0;
  auto k = con.GetParameters()[0];
  for (auto v: con.GetArguments()) {
    if ((x.is_var_int(v)
         && std::round(x[v]) == k)
        || std::fabs(x[v] - k) <= x.feastol())
      ++result;
  }
  return result;
}

/// Compute result of the NumberofVar constraint.
template <class VarVec>
double ComputeValue(const NumberofVarConstraint& con, const VarVec& x) {
  double result = 0.0;
  const auto& args = con.GetArguments();
  auto k = x[args[0]];
  for (auto i=args.size(); --i; ) {
    auto v = args[i];
    if ((x.is_var_int(v)
         && std::round(x[v]) == k)
        || std::fabs(x[v] - k) <= x.feastol())
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

/// Compute result of the asinh constraint.
template <class VarVec>
double ComputeValue(const AsinhConstraint& con, const VarVec& x) {
  return std::asinh(x[con.GetArguments()[0]]);
}

/// Compute result of the acosh constraint.
template <class VarVec>
double ComputeValue(const AcoshConstraint& con, const VarVec& x) {
  return std::acosh(x[con.GetArguments()[0]]);
}

/// Compute result of the atanh constraint.
template <class VarVec>
double ComputeValue(const AtanhConstraint& con, const VarVec& x) {
  return std::atanh(x[con.GetArguments()[0]]);
}

/// Compute result of the LinearFuncCon constraint.
template <class VarVec>
double ComputeValue(
    const LinearFunctionalConstraint& con, const VarVec& x) {
  return con.GetAffineExpr().ComputeValue(x);
}

/// Compute result of the QuadrFuncCon constraint.
template <class VarVec>
double ComputeValue(
    const QuadraticFunctionalConstraint& con, const VarVec& x) {
  return con.GetQuadExpr().ComputeValue(x);
}



/// Compute result of a conditional constraint.
/// Just return bool(viol(subcon) <= 0).
/// This is not used to compute violation.
template <class Con, class VarVec>
double ComputeValue(
    const ConditionalConstraint<Con>& con, const VarVec& x) {
  auto viol = con.GetConstraint().ComputeViolation(x);
  bool ccon_valid = viol<=0.0;
  return ccon_valid;
}

/// Compute violation of the QuadraticCone constraint.
template <class VarVec>
double ComputeViolation(
    const QuadraticConeConstraint& con, const VarVec& x) {
  const auto& args = con.GetArguments();
  const auto& params = con.GetParameters();
  assert(args.size()==params.size());
  double sum = 0.0;
  for (auto i=args.size(); --i; )
    sum += std::pow( params[i]*x[args[i]], 2.0 );
  return std::sqrt(sum) - params[0]*x[args[0]];
}

/// Compute violation of the RotatedQuadraticCone constraint.
template <class VarVec>
double ComputeViolation(
    const RotatedQuadraticConeConstraint& con, const VarVec& x) {
  const auto& args = con.GetArguments();
  const auto& params = con.GetParameters();
  assert(args.size()==params.size());
  double sum = 0.0;
  for (auto i=args.size(); --i>1; )
    sum += std::pow( params[i]*x[args[i]], 2.0 );
  return sum
      - 2.0 * params[0]*x[args[0]] * params[1]*x[args[1]];
}

/// Compute violation of the ExponentialCone constraint.
template <class VarVec>        // ax >= by exp(cz / (by))
double ComputeViolation(       // where ax, by >= 0
    const ExponentialConeConstraint& con, const VarVec& x) {
  const auto& args = con.GetArguments();
  const auto& params = con.GetParameters();
  auto ax = params[0]*x[args[0]];
  auto by = params[1]*x[args[1]];
  if (0.0==std::fabs(by))
    return -ax;
  auto cz = params[2]*x[args[2]];
  return by * std::exp(cz / by) - ax;
}

/// Compute result of the PL constraint.
template <class VarVec>
double ComputeValue(const PLConstraint& con, const VarVec& x) {
  const auto& plp = con.GetParameters().GetPLPoints();
  assert(!plp.empty());
  auto x0 = x[con.GetArguments()[0]];        // position
  if (x0<plp.x_.front())
    return plp.y_.front()
        - plp.PreSlope()*(plp.x_.front() - x0);
  if (x0>plp.x_.back())
    return plp.y_.back()
        + plp.PostSlope()*(x0 - plp.x_.back());
  int i0=0;
  for ( ; x0 > plp.x_[i0]; ++i0) ;
  return plp.x_[i0]==x0
      ? plp.y_[i0]
        : (plp.y_[i0-1]
           + (plp.y_[i0]-plp.y_[i0-1])
        * (x0-plp.x_[i0-1]) / (plp.x_[i0]-plp.x_[i0-1]));
}

/// Should be here,
/// after ComputeViolation() is specialized
/// for some constraints.
template <class Args, class Params,
          class NumOrLogic, class Id>
template <class VarVec>
double
CustomFunctionalConstraint<Args, Params, NumOrLogic, Id>
::ComputeViolation(const VarVec& x) const
{ return mp::ComputeViolation(*this, x); }


}  // namespace mp

#endif // CONSTR_EVAL_H
