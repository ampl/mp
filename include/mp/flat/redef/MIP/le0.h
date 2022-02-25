#ifndef LE0_H
#define LE0_H

#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constraints_std.h"

namespace mp {

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
    const auto& ae = le0c.GetArguments();
    const auto res = le0c.GetResultVar();
    if (ae.is_constant()) {
      if (ae.constant_term() > 0.0)
        GetMC().NarrowVarBounds(res, 0.0, 0.0);
    } else {
      if (GetMC().is_fixed(res)) {
        if (GetMC().fixed_value(res)) {  // fixed to 1
          GetMC().AddConstraint( ExtractConstraint(le0c) );
        }
      } else {
        GetMC().AddConstraint(IndicatorConstraintLinLE(
                                res, 1,
                                ExtractConstraint(le0c)));
      }
    }
  }

  /// Convert in negative context
  /// b==0 --> c'x > d
  void ConvertCtxNeg(const ItemType& le0c, int ) {
    auto ae = le0c.GetArguments();
    const auto res = le0c.GetResultVar();
    if (ae.is_constant()) {
      if (ae.constant_term() <= 0.0)
        GetMC().NarrowVarBounds(res, 1.0, 1.0);
    } else {
      if (GetMC().is_fixed(res)) {
        if (!GetMC().fixed_value(res)) {       // fixed to 0
          GetMC().AddConstraint(LinConGE(
                                  LinTerms(ae),
                                  -ae.constant_term()+1));
        }
      } else {
        ae.negate();
        auto bNt = GetMC().ComputeBoundsAndType(ae);
        double cmpEps = GetMC().ComparisonEps( bNt.get_result_type() );
        double d = ae.constant_term() + cmpEps;
        GetMC().AddConstraint(IndicatorConstraintLinLE(
                                res, 0,
                                { LinTerms(ae), -d }));
      }
    }
  }

  /// Reuse the stored ModelConverter
  using Base::GetMC;
};

} // namespace mp

#endif // LE0_H
