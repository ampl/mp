#ifndef LOGICAL_NOT_H
#define LOGICAL_NOT_H

#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constr_std.h"

namespace mp {

/// Converts NOT for MIP
template <class ModelConverter>
class NotConverter_MIP :
    public BasicFuncConstrCvt<
      NotConverter_MIP<ModelConverter>, ModelConverter> {
public:
  /// Base class
  using Base = BasicFuncConstrCvt<
    NotConverter_MIP<ModelConverter>, ModelConverter>;
  /// Constructor
  NotConverter_MIP(ModelConverter& mc) : Base(mc) { }
  /// Converted item type
  using ItemType = NotConstraint;

  /// Convert in both contexts (full reification)
  Context Convert(const ItemType& nc, int ) {
    if (GetMC().is_fixed(nc.GetResultVar())) {           // fixed result?
      assert(GetMC().is_fixed(nc.GetArguments()[0]));    // propagated down
      assert(GetMC().fixed_value(nc.GetResultVar())
             == 1.0 - GetMC().fixed_value(nc.GetArguments()[0]));
      return Context::CTX_MIX;
    }
    if (GetMC().is_fixed(nc.GetArguments()[0])) {
      auto resval = 1.0 - GetMC().fixed_value(nc.GetArguments()[0]);
      GetMC().NarrowVarBounds(nc.GetResultVar(), resval, resval);
      return Context::CTX_MIX;
    }
    LinearFunctionalConstraint funccon {{{{-1.0}, {nc.GetArguments()[0]}}, 1.0}};
    funccon.SetContext( GetMC().GetInitExprContext(nc.GetResultVar()) );
    GetMC().RedefineVariable(nc.GetResultVar(), std::move(funccon));
    // Old way, no expression: /// Obtain negation variable via map
    // int var_res_lin = GetMC().AssignResultVar2Args(
    //       LinearFunctionalConstraint(
    //         {{{-1.0}, {nc.GetArguments()[0]}}, 1.0}));
    // GetMC().AddConstraint_AS_ROOT(LinConEQ{    // Could use RedefineVariable()
    //                         { {-1.0, 1.0},
    //                           {nc.GetResultVar(), var_res_lin} },
    //                         {0.0}});
    return Context::CTX_MIX;
  }

  /// Reuse the stored ModelConverter
  using Base::GetMC;
};

} // namespace mp

#endif // LOGICAL_NOT_H
