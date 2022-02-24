#ifndef ALLDIFF_H
#define ALLDIFF_H

/*
 * Convert AllDiff for MIP
 */

#include "mp/flat/redef_base.h"
#include "mp/flat/constraints_std.h"

namespace mp {

/// Converts proper range linear constraints to c'x-slack=ub,
/// otherwise to c'x ? rhs.
template <class ModelConverter>
class AllDiffConverter_MIP :
    public BasicFuncConstrCvt<
      AllDiffConverter_MIP<ModelConverter>, ModelConverter> {
public:
  /// Base class
  using Base = BasicFuncConstrCvt<
    AllDiffConverter_MIP<ModelConverter>, ModelConverter>;
  /// Constructor
  AllDiffConverter_MIP(ModelConverter& mc) : Base(mc) { }
  /// Converted item type
  using ItemType = AllDiffConstraint;

  /// Convert in positive context
  /// TODO add presolve bridge?
  void ConvertCtxPos(const ItemType& alld, int ) {
    const auto& args = alld.GetArguments();
    const auto lba_dbl = GetMC().lb_array(args);
    const auto uba_dbl = GetMC().ub_array(args);
    if (lba_dbl<=GetMC().MinusInfty() ||
        uba_dbl>=GetMC().Infty())
      MP_RAISE("MP2MIP: AllDiff on unbounded variables not implemented");
    if (lba_dbl<std::numeric_limits<int>::min() ||
        uba_dbl>std::numeric_limits<int>::max())
      MP_RAISE("MP2MIP: AllDiff on variables with domain "
               "out of integer range not implemented");
    const int lba = (int)std::round(lba_dbl);
    const int uba = (int)std::round(uba_dbl);
    std::vector<double> coefs(args.size(), 1.0);
    std::vector<int> flags(args.size());        // unary encoding flags
    double rhs = 1.0;
    if (!GetMC().is_fixed(alld.GetResultVar())) {  // implied version: b -> alldiff
      coefs.push_back(args.size()-1.0);
      flags.push_back(alld.GetResultVar());
      rhs = (double)args.size();
    }
    for (int v=lba; v!=uba+1; ++v) {            // for each value in the domain union
      for (size_t ivar = 0; ivar < args.size(); ++ivar) {
        flags[ivar] = GetMC().AssignResultVar2Args(
              EQ0Constraint( { {{1.0}, {args[ivar]}}, -double(v) } ) );
      }
      GetMC().AddConstraint( LinConLE( {coefs, flags}, {rhs} ) );
    }
  }

  /// Reuse the stored ModelConverter
  using Base::GetMC;
};

} // namespace mp

#endif // ALLDIFF_H
