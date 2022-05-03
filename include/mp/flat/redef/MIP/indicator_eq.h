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
    auto bnds = GetMC().ComputeBoundsAndType(
          indc.get_constraint().GetBody());
    /// Converting b==val ==> c'x==d to
    ///   ==> c'x <= d and
    ///   ==> c'x >= d
    auto con = LinConLE{ indc.get_constraint().GetBody(),
        indc.get_constraint().rhs() };
    ConvertImplicationLE(binvar, indc.get_binary_value(),
                         bnds.ub(), con);
    con.negate();
    bnds.NegateBounds();
    ConvertImplicationLE(binvar, indc.get_binary_value(),
                         bnds.ub(), std::move(con));
  }

protected:
  using Base::GetMC;
  using Base::ConvertImplicationLE;
};

} // namespace mp

#endif // INDICATOR_EQ_H
