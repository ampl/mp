#ifndef IMPL_H
#define IMPL_H

/**
 * Redefinition for MIP of
 *   (cond) ==> (constr1) [ else (constr2) ]
 */

#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constr_std.h"

namespace mp {

/// Convert ImplicationConstraint into AndConstraint
template <class ModelConverter>
class ImplicationConverter_MIP :
    public BasicFuncConstrCvt<
    ImplicationConverter_MIP<ModelConverter>, ModelConverter> {
public:
  /// Base class
  using Base = BasicFuncConstrCvt<
  ImplicationConverter_MIP<ModelConverter>, ModelConverter>;
  /// Constructor
  ImplicationConverter_MIP(ModelConverter& mc) : Base(mc) { }
  /// Converted item type
  using ItemType = ImplicationConstraint;

  /// Convert in any context
  void Convert(const ItemType& cc, int ) {
    const auto& args = cc.GetArguments();
    auto compl_arg0 = GetMC().MakeComplementVar(args[0]);
    /// args[0] ==> args[1]
    auto disj1 = GetMC().AssignResultVar2Args(
          OrConstraint{ {compl_arg0, args[1]} });
    /// !args[0] ==> args[2]
    auto disj2 = GetMC().AssignResultVar2Args(
          OrConstraint{ {args[0], args[2]} });
    if (GetMC().is_fixed(cc.GetResultVar())   // the impl is static
        && 1 == GetMC().fixed_value(cc.GetResultVar())) {
      GetMC().FixAsTrue(disj1);   // XPRESS 9.0 needs this form
      GetMC().FixAsTrue(disj2);
    } else { // Redefine constraint which can be the init expr
      GetMC().RedefineVariable(cc.GetResultVar(),
                               AndConstraint{ {disj1, disj2} });
      GetMC().PropagateResultOfInitExpr(
            cc.GetResultVar(), cc.GetContext());
    }
  }


protected:
  /// Reuse the stored ModelConverter
  using Base::GetMC;
};

} // namespace mp

#endif // IMPL_H
