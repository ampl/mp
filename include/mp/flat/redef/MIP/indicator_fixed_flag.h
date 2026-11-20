#ifndef INDICATOR_FIXED_FLAG_H
#define INDICATOR_FIXED_FLAG_H

#include "mp/flat/redef/redef_base.h"

namespace mp {

/// Base of the indicator converters:
/// folds an indicator whose binary is fixed.
///
/// The static presolve (PreprocessStaticConstraint())
/// folds such an indicator when it is added.
/// A binary can be fixed later, e.g., by a root logical
/// constraint flattened afterwards. Then the indicator
/// would reach the solver with a constant binary,
/// which some solvers reject once the binary has been made
/// continuous (cvt:pre:continuous_fixed_vars).
/// Hence the item is reconsidered in conversion:
/// IfNeedsConversion() also makes accepted indicators convert,
/// and FlatConverter::ConvertDelayedItems() revisits
/// items behind the conversion cursor.
/// Controlled by cvt:pre:indfixedflag.
template <class ModelConverter>
class BasicIndicatorConverter_MIP :
    public BasicItemConverter<ModelConverter> {
public:
  /// Base class
  using Base = BasicItemConverter<ModelConverter>;
  /// Constructor
  BasicIndicatorConverter_MIP(ModelConverter& mc) : Base(mc) { }

  /// An accepted indicator needs conversion iff its
  /// binary is fixed (and folding is enabled)
  template <class Indicator>
  bool IfNeedsConversion(const Indicator& indc, int ) {
    return GetMC().IfIndicatorFlagFoldable(indc.get_binary_var());
  }

protected:
  using Base::GetMC;

  /// Fold \a indc if its binary is fixed:
  /// add the subconstraint if the binary has the indicator's value,
  /// otherwise drop the indicator.
  /// @return true iff folded
  template <class Indicator>
  bool FoldIfFlagFixed(const Indicator& indc) {
    if (!GetMC().IfIndicatorFlagFoldable(indc.get_binary_var()))
      return false;
    const bool folded = GetMC().PreprocessStaticConstraint(indc);
    assert(folded);
    return folded;
  }
};

} // namespace mp

#endif // INDICATOR_FIXED_FLAG_H
