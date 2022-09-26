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
template <class Impl, class ModelConverter, class FuncCon>
class FuncConConverter_MIP_CRTP :
    public BasicFuncConstrCvt<Impl, ModelConverter> {
public:
  /// Base class
  using Base = BasicFuncConstrCvt<Impl, ModelConverter>;
  /// Constructor
  FuncConConverter_MIP_CRTP(ModelConverter& mc) : Base(mc) { }
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


/// Final CRTP class
template <class ModelConverter, class FuncCon>
class FuncConConverter_MIP :
    public FuncConConverter_MIP_CRTP<
      FuncConConverter_MIP<ModelConverter, FuncCon>,
    ModelConverter, FuncCon> {
public:
  /// Typedef Base
  using Base = FuncConConverter_MIP_CRTP<
        FuncConConverter_MIP<ModelConverter, FuncCon>,
      ModelConverter, FuncCon>;
  /// Constructor
  FuncConConverter_MIP(ModelConverter& mc) : Base(mc) { }
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

/// Typedef FuncConConverter_MIP_Cos
template <class MC>
using FuncConConverter_MIP_Cos = FuncConConverter_MIP<MC, CosConstraint>;

/// Typedef FuncConConverter_MIP_Tan
template <class MC>
using FuncConConverter_MIP_Tan = FuncConConverter_MIP<MC, TanConstraint>;

/// Typedef FuncConConverter_MIP_ASin
template <class MC>
using FuncConConverter_MIP_Asin = FuncConConverter_MIP<MC, AsinConstraint>;

/// Typedef FuncConConverter_MIP_ACos
template <class MC>
using FuncConConverter_MIP_Acos = FuncConConverter_MIP<MC, AcosConstraint>;

/// Typedef FuncConConverter_MIP_ATan
template <class MC>
using FuncConConverter_MIP_Atan = FuncConConverter_MIP<MC, AtanConstraint>;

/// Typedef FuncConConverter_MIP_Sinh
template <class MC>
using FuncConConverter_MIP_Sinh = FuncConConverter_MIP<MC, SinhConstraint>;

/// Typedef FuncConConverter_MIP_Cosh
template <class MC>
using FuncConConverter_MIP_Cosh = FuncConConverter_MIP<MC, CoshConstraint>;

/// Typedef FuncConConverter_MIP_Tanh
template <class MC>
using FuncConConverter_MIP_Tanh = FuncConConverter_MIP<MC, TanhConstraint>;

/// Typedef FuncConConverter_MIP_ASinh
template <class MC>
using FuncConConverter_MIP_Asinh = FuncConConverter_MIP<MC, AsinhConstraint>;

/// Typedef FuncConConverter_MIP_ACosh
template <class MC>
using FuncConConverter_MIP_Acosh = FuncConConverter_MIP<MC, AcoshConstraint>;

/// Typedef FuncConConverter_MIP_ATanh
template <class MC>
using FuncConConverter_MIP_Atanh = FuncConConverter_MIP<MC, AtanhConstraint>;


} // namespace mp

#endif // LIN_APPROX_H
