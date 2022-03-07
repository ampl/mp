#ifndef LT0_H
#define LT0_H

#include "mp/flat/redef/MIP/le0.h"

namespace mp {

/// Converts LT0Constraint for MIP
template <class ModelConverter>
class LT0Converter_MIP :
    public BasicFuncConstrCvt<
      LT0Converter_MIP<ModelConverter>, ModelConverter> {
public:
  /// Base class
  using Base = BasicFuncConstrCvt<
    LT0Converter_MIP<ModelConverter>, ModelConverter>;
  /// Constructor
  LT0Converter_MIP(ModelConverter& mc) : Base(mc) { }
  /// Converted item type
  using ItemType = LT0Constraint;

  /// Convert in positive context
  void ConvertCtxPos(const ItemType& le0c, int ) {
    auto bNt = GetMC().ComputeBoundsAndType(le0c.GetArguments());
    double cmpEps = GetMC().ComparisonEps( bNt.get_result_type() );
    ConvertLELT0_MIP_CtxPos(GetMC(), le0c, cmpEps);
  }

  /// Convert in negative context.
  /// resvar==0 --> c'x >(=) d
  void ConvertCtxNeg(const ItemType& le0c, int ) {
    ConvertLELT0_MIP_CtxNeg(GetMC(), le0c, 0.0);
  }

  /// Reuse the stored ModelConverter
  using Base::GetMC;
};

} // namespace mp

#endif // LT0_H
