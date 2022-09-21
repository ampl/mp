#ifndef LIN_APPROX_H
#define LIN_APPROX_H

/** Redefinition of smooth nonlinear functions
 *  via piecewise-linear approximation
 */

#include "mp/common.h"

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
    PLApproxParams laPrm;
    /// Narrow graph domain to +-1e6
    laPrm.grDom.lbx = std::max(GetMC().lb(x), -1e6);
    laPrm.grDom.ubx = std::min(GetMC().ub(x), 1e6);
    laPrm.grDom.lby = std::max(GetMC().lb(y), -1e6);
    laPrm.grDom.uby = std::min(GetMC().ub(y), 1e6);

    PLApproximate(con, laPrm);
    if (!laPrm.fUsePeriod) {
      GetMC().RedefineVariable(y,
                               PLConstraint({x}, laPrm.plPoints));
      GetMC().PropagateResultOfInitExpr(
            // propagate ctx into new constr
            con.GetResultVar(), con.GetContext());
    } else {
      auto rmd = GetMC().AddVar(
            laPrm.periodRemainderRange.lb,
            laPrm.periodRemainderRange.ub);
      auto factor = GetMC().AddVar(
            laPrm.periodicFactorRange.lb,
            laPrm.periodicFactorRange.ub,
            var::INTEGER);
      GetMC().RedefineVariable(y,
                               PLConstraint({rmd}, laPrm.plPoints));
      GetMC().PropagateResultOfInitExpr(
            // propagate ctx into new constr
            con.GetResultVar(), con.GetContext());
      GetMC().AddConstraint( LinConEQ{  // x = period * factor + rmd
                               { {laPrm.periodLength, 1.0, -1.0},
                                 {factor, rmd, x} },
                               {0.0} } );
    }
  }


protected:
  /// Reuse the stored ModelConverter
  using Base::GetMC;
};


/// Typedef FuncConConverter_MIP_Exp
template <class MC>
using FuncConConverter_MIP_Exp = FuncConConverter_MIP<MC, ExpConstraint>;

/// Typedef FuncConConverter_MIP_Log
template <class MC>
using FuncConConverter_MIP_Log = FuncConConverter_MIP<MC, LogConstraint>;

/// Typedef FuncConConverter_MIP_LogA
template <class MC>
using FuncConConverter_MIP_LogA = FuncConConverter_MIP<MC, LogAConstraint>;

/// Typedef FuncConConverter_MIP_ExpA
template <class MC>
using FuncConConverter_MIP_ExpA = FuncConConverter_MIP<MC, ExpAConstraint>;

/// Typedef FuncConConverter_MIP_Sin
template <class MC>
using FuncConConverter_MIP_Sin = FuncConConverter_MIP<MC, SinConstraint>;


} // namespace mp

#endif // LIN_APPROX_H
