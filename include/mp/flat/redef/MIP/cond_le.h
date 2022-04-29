#ifndef COND_LE_H
#define COND_LE_H

#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constraints_std.h"

namespace mp {

/// Underlying algebraic constraint of a conditional constraint,
/// parameterized by comparison type
template <class CondCon, int kind>
using UndAlgCon =
  AlgebraicConstraint<
    typename CondCon::ConType::BodyType, AlgConRhs<kind> >;

/// Specialize for LE
template <class CondCon>
using UndAlgConLE = UndAlgCon<CondCon, -1>;

/// Convert LE / LT constraints for MIP,
/// positive context: res==1 ==> (body) <= rhs-eps
/// @param mc: ModelConverter
/// @param lelt: LE / LTConstraint
/// @param eps: rhs modifier
template <class ModelConverter, class CondCon>
void ConvertLELT_MIP_CtxPos(ModelConverter& mc,
                            const CondCon& lelt, double eps) {
  const auto& con = lelt.GetConstraint();
  const auto res = lelt.GetResultVar();
  if (con.empty()) {
    if (con.rhs() < eps)           // TODO option to switch off
      mc.NarrowVarBounds(res, 0.0, 0.0);
  } else {
    if (mc.is_fixed(res)) {
      if (mc.fixed_value(res)) {       // res==1
        mc.AddConstraint(
              UndAlgConLE<CondCon>{ con.GetBody(), con.rhs()-eps } );
      }
    } else {
      mc.AddConstraint(IndicatorConstraint< UndAlgConLE<CondCon> >(
                              res, 1,
                              { con.GetBody(), con.rhs()-eps }));
    }
  }
}

/// Convert LE / LT constraints for MIP,
/// negative context.
/// resvar==0 --> (body) >= d+eps
/// @param mc: ModelConverter
/// @param lelt: LE / LTConstraint
/// @param eps: rhs modifier
template <class ModelConverter, class CondCon>
void ConvertLELT_MIP_CtxNeg(ModelConverter& mc,
                             const CondCon& lelt, double eps) {
  const auto& con = lelt.GetConstraint();
  const auto res = lelt.GetResultVar();
  if (con.empty()) {
    if (con.rhs() > -eps)
      mc.NarrowVarBounds(res, 1.0, 1.0);
  } else {
    auto con = lelt.GetConstraint();
    con.negate();
    con.set_rhs(con.rhs() - eps);
    if (mc.is_fixed(res)) {
      if (!mc.fixed_value(res)) {            // fixed to 0
        mc.AddConstraint(
              UndAlgConLE<CondCon>{ std::move(con.GetBody()), con.rhs() });
      }
    } else {
      mc.AddConstraint(IndicatorConstraint< UndAlgConLE<CondCon> >(
                              res, 0, { std::move(con.GetBody()), con.rhs() }));
    }
  }
}


/// Converts conditional <= for MIP
template <class ModelConverter, class AlgConBody>
class CondLEConverter_MIP :
    public BasicFuncConstrCvt<
      CondLEConverter_MIP<ModelConverter, AlgConBody>,
      ModelConverter> {
public:
  /// Base class
  using Base = BasicFuncConstrCvt<
    CondLEConverter_MIP<ModelConverter, AlgConBody>,
    ModelConverter>;

  /// Constructor
  CondLEConverter_MIP(ModelConverter& mc) : Base(mc) { }

  /// Underlying algebraic constraint,
  /// parameterized by comparison type
  template <int kind>
  using AlgCon =
    AlgebraicConstraint< AlgConBody, AlgConRhs<kind> >;

  /// Converted item type
  using ItemType = ConditionalConstraint< AlgCon<-1> >;

  /// Convert in positive context
  void ConvertCtxPos(const ItemType& le0c, int ) {
    ConvertLELT_MIP_CtxPos(GetMC(), le0c, 0.0);
  }

  /// Convert in negative context.
  /// resvar==0 --> c'x >(=) d
  void ConvertCtxNeg(const ItemType& le0c, int ) {
    auto bNt = GetMC().ComputeBoundsAndType(le0c.GetArguments());
    double cmpEps = GetMC().ComparisonEps( bNt.get_result_type() );
    ConvertLELT_MIP_CtxNeg(GetMC(), le0c, cmpEps);
  }

  /// Reuse the stored ModelConverter
  using Base::GetMC;
};


/// Specialize for linear constraint
template <class MC>
using CondLinLEConverter_MIP = CondLEConverter_MIP<MC, LinTerms>;

/// Specialize for quadratic constraint
template <class MC>
using CondQuadLEConverter_MIP = CondLEConverter_MIP<MC, QuadAndLinTerms>;

} // namespace mp

#endif // COND_LE_H
