#ifndef LIN_APPROX_H
#define LIN_APPROX_H

/** Redefinition of smooth nonlinear functions
 *  via piecewise-linear approximation
 */

#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constr_std.h"
#include "mp/flat/redef/MIP/core/lin_approx_core.h"

namespace mp {

/// Pl-approximates FuncCon for MIP
template <class ModelConverter, class FuncCon>
class FuncConConverter_MIP :
    public BasicFuncConstrCvt<
      FuncConConverter_MIP<ModelConverter, FuncCon>,
      ModelConverter> {
public:
  /// Base class
  using Base = BasicFuncConstrCvt<
    FuncConConverter_MIP<ModelConverter, FuncCon>, ModelConverter>;
  /// Constructor
  FuncConConverter_MIP(ModelConverter& mc) : Base(mc) { }
  /// Converted item type
  using ItemType = FuncCon;

  /// Convert in any context
  void Convert(const ItemType& con, int ) {
    assert(!con.GetContext().IsNone());
    assert(1==con.GetArguments().size());          // 1 argument var
    auto x = con.GetArguments()[0];
    auto y = con.GetResultVar();
    LinApproxParams laPrm;
    laPrm.lbx = std::max(GetMC().lb(x), -1e6);
    laPrm.ubx = std::min(GetMC().ub(x), 1e6);
    laPrm.lby = std::max(GetMC().lb(y), -1e6);
    laPrm.uby = std::min(GetMC().ub(y), 1e6);
    GetMC().RedefineVariable(con.GetResultVar(),
          PLConstraint({x}, PLApproximate(con, laPrm)));
    GetMC(). PropagateResultOfInitExpr(     // propagate ctx into new constr
          con.GetResultVar(), con.GetContext());
  }


protected:
  /// Reuse the stored ModelConverter
  using Base::GetMC;
};

} // namespace mp

#endif // LIN_APPROX_H
