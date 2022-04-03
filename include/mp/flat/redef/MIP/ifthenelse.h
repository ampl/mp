#ifndef IFTHENELSE_H
#define IFTHENELSE_H

#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constraints_std.h"

namespace mp {

/// Converts ITE for MIP
template <class ModelConverter>
class IfThenElseConverter_MIP :
    public BasicFuncConstrCvt<
      IfThenElseConverter_MIP<ModelConverter>, ModelConverter> {
public:
  /// Base class
  using Base = BasicFuncConstrCvt<
    IfThenElseConverter_MIP<ModelConverter>, ModelConverter>;
  /// Constructor
  IfThenElseConverter_MIP(ModelConverter& mc) : Base(mc) { }
  /// Converted item type
  using ItemType = IfThenConstraint;

  /// Convert in any context
  void Convert(const ItemType& itc, int ) {
    assert(!itc.GetContext().IsNone());
    const auto& args = itc.GetArguments();
    if (!GetMC().is_fixed(args[1]) || !GetMC().is_fixed(args[2]))
      MP_RAISE("MP2MIP: IfThen with variable then/else "
               "arguments not implemented");
    else
      ConvertIfThen_constantThenElse(itc);
  }

protected:
  void ConvertIfThen_constantThenElse(const IfThenConstraint& itc) {
    const auto& args = itc.GetArguments();
    assert((GetMC().is_fixed(args[1]) && GetMC().is_fixed(args[2])));
    const double const1 = GetMC().fixed_value(args[1]);
    const double const2 = GetMC().fixed_value(args[2]);
    /// Obtain negation variable via map
    int var_res_lin = GetMC().AssignResultVar2Args(
          LinearFunctionalConstraint(
            { {{const1-const2}, {args[0]}}, const2 } ));
    GetMC().AddConstraint(LinConEQ{
                            { {-1.0, 1.0},
                              {itc.GetResultVar(), var_res_lin} },
                            {0.0}});
  }

  /// Reuse the stored ModelConverter
  using Base::GetMC;
};

} // namespace mp

#endif // IFTHENELSE_H
