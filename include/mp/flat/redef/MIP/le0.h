#ifndef LE0_H
#define LE0_H

#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constraints_std.h"

namespace mp {

/// Convert LE0 / LT0 constraints for MIP,
/// positive context.
/// @param mc: ModelConverter
/// @param lelt0: LE0 / LT0Constraint
/// @param eps: rhs modifier
template <class ModelConverter, class Con>
void ConvertLELT0_MIP_CtxPos(ModelConverter& mc,
                             const Con& lelt0, double eps) {
  const auto& ae = lelt0.GetArguments();
  const auto res = lelt0.GetResultVar();
  if (ae.is_constant()) {
    if (ae.constant_term()+eps > 0.0)
      mc.NarrowVarBounds(res, 0.0, 0.0);
  } else {
    if (mc.is_fixed(res)) {
      if (mc.fixed_value(res)) {  // fixed to 1
        mc.AddConstraint( ExtractConstraint(lelt0, eps) );
      }
    } else {
      mc.AddConstraint(IndicatorConstraintLinLE(
                              res, 1,
                              ExtractConstraint(lelt0, eps)));
    }
  }
}

/// Convert LE0 / LT0 constraints for MIP,
/// negative context.
/// resvar==0 --> c'x >(=) d
/// @param mc: ModelConverter
/// @param lelt0: LE0 / LT0Constraint
/// @param eps: rhs modifier
template <class ModelConverter, class Con>
void ConvertLELT0_MIP_CtxNeg(ModelConverter& mc,
                             const Con& lelt0, double eps) {
  auto ae = lelt0.GetArguments();
  const auto res = lelt0.GetResultVar();
  if (ae.is_constant()) {
    if (ae.constant_term() <= 0.0)
      mc.NarrowVarBounds(res, 1.0, 1.0);
  } else {
    if (mc.is_fixed(res)) {
      if (!mc.fixed_value(res)) {       // fixed to 0
        mc.AddConstraint(LinConGE(
                                LinTerms(ae),
                                -ae.constant_term()+eps));
      }
    } else {
      ae.negate();
      double d = ae.constant_term() + eps;
      mc.AddConstraint(IndicatorConstraintLinLE(
                              res, 0,
                              { LinTerms(ae), -d }));
    }
  }
}


/// Converts LE0Constraint for MIP
template <class ModelConverter>
class LE0Converter_MIP :
    public BasicFuncConstrCvt<
      LE0Converter_MIP<ModelConverter>, ModelConverter> {
public:
  /// Base class
  using Base = BasicFuncConstrCvt<
    LE0Converter_MIP<ModelConverter>, ModelConverter>;
  /// Constructor
  LE0Converter_MIP(ModelConverter& mc) : Base(mc) { }
  /// Converted item type
  using ItemType = LE0Constraint;

  /// Convert in positive context
  void ConvertCtxPos(const ItemType& le0c, int ) {
    ConvertLELT0_MIP_CtxPos(GetMC(), le0c, 0.0);
  }

  /// Convert in negative context.
  /// resvar==0 --> c'x >(=) d
  void ConvertCtxNeg(const ItemType& le0c, int ) {
    auto bNt = GetMC().ComputeBoundsAndType(le0c.GetArguments());
    double cmpEps = GetMC().ComparisonEps( bNt.get_result_type() );
    ConvertLELT0_MIP_CtxNeg(GetMC(), le0c, cmpEps);
  }

  /// Reuse the stored ModelConverter
  using Base::GetMC;
};

} // namespace mp

#endif // LE0_H
