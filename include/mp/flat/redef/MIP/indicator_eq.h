#ifndef INDICATOR_EQ_H
#define INDICATOR_EQ_H

#include "mp/flat/redef/MIP/indicator_le.h"

namespace mp {

/// Convert IndicatorLinLE.
/// b==val ==> c'x==d
template <class ModelConverter>
class IndicatorLinEQConverter_MIP :
    public IndicatorLinLEConverter_MIP<ModelConverter> {
public:
  /// Base class
  using Base = IndicatorLinLEConverter_MIP<ModelConverter>;
  /// Constructor
  IndicatorLinEQConverter_MIP(ModelConverter& mc) : Base(mc) { }
  /// Converted item type
  using ItemType = IndicatorConstraintLinEQ;

  /// Conversion
  void Convert(const ItemType& indc, int ) {
    auto binvar=indc.get_binary_var();
    auto ae = indc.to_lhs_expr();
    auto bnds = GetMC().ComputeBoundsAndType(ae);
    ConvertImplicationLE(binvar, indc.get_binary_value(),
                         bnds, ae);
    ae.negate();
    bnds.NegateBounds();
    ConvertImplicationLE(binvar, indc.get_binary_value(),
                         bnds, std::move(ae));
  }

protected:
  using Base::GetMC;
  using Base::ConvertImplicationLE;
};

} // namespace mp

#endif // INDICATOR_EQ_H
