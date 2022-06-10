#ifndef COMPLEMENT_H
#define COMPLEMENT_H

#include <cmath>

#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constr_std.h"

namespace mp {

/// Converts Complementarity for MIP
template <class ModelConverter, class ComplCon>
class ComplementarityConverter_MIP :
    public BasicFuncConstrCvt<
      ComplementarityConverter_MIP<ModelConverter, ComplCon>,
      ModelConverter> {
public:
  /// Base class
  using Base = BasicFuncConstrCvt<
    ComplementarityConverter_MIP<ModelConverter, ComplCon>,
    ModelConverter>;

  /// Constructor
  ComplementarityConverter_MIP(ModelConverter& mc) : Base(mc) { }

  /// Converted item type
  using ItemType = ComplCon;

  /// Convert in any context
  void Convert(const ItemType& cc, int ) {
    const auto& expr = cc.GetExpression();
    auto compl_var = cc.GetVariable();

    double
        var_lb = GetMC().lb(compl_var),
        var_ub = GetMC().ub(compl_var);

    bool
        fin_var_lb = std::isfinite(var_lb),
        fin_var_ub = std::isfinite(var_ub);

    /// Using algebraic expression (expr.body + 0.0)
    auto expr_var = GetMC().Convert2Var(
          AlgebraicExpression<typename ComplCon::ExprType::BodyType>{
            expr.GetBody(), 0.0} );
    double con_rhs = -expr.constant_term();

    if (fin_var_lb && !fin_var_ub) {
      /// res1 = (var <= var_lb)
      auto res_neg_var_lb = GetMC().AssignResultVar2Args(
            CondLinConLE{ { {{1.0}, {compl_var}}, var_lb } });
      /// res2 = (body <= lb)
      auto res_neg_con_lb = GetMC().AssignResultVar2Args(
            CondLinConLE{ { {{1.0}, {expr_var}}, con_rhs } });
      /// res3 = (res1 \/ res2)
      auto res_disj = GetMC().AssignResultVar2Args(
            OrConstraint{ { res_neg_var_lb, res_neg_con_lb } });
      GetMC().FixAsTrue(res_disj);
      /// Add the algebraic constraint via the representing variable
      /// Propagate mixed context (logical constraint would set CTX_NEG)
      GetMC().set_var_lb_context(expr_var, con_rhs, Context::CTX_MIX);
    } else if (fin_var_ub && !fin_var_lb) {
      /// res1 = (var >= var_ub)
      auto res_neg_var_ub = GetMC().AssignResultVar2Args(
            CondLinConGE{ { {{1.0}, {compl_var}}, var_ub } });
      /// res2 = (body >= ub)
      auto res_neg_con_ub = GetMC().AssignResultVar2Args(
            CondLinConGE{ { {{1.0}, {expr_var}}, con_rhs } });
      /// res3 = (res1 \/ res2)
      auto res_disj = GetMC().AssignResultVar2Args(
            OrConstraint{ { res_neg_var_ub, res_neg_con_ub } });
      GetMC().FixAsTrue(res_disj);
      /// Add the algebraic constraint via the representing variable
      /// Propagate mixed context (logical constraint would set CTX_POS)
      GetMC().set_var_ub_context(expr_var, con_rhs, Context::CTX_MIX);
    } else {
      assert(fin_var_lb && fin_var_ub);
      /// res1 = (var <= lb && con >= 0)
      auto res1 = GetMC().AssignResultVar2Args(
            AndConstraint{ {
                GetMC().AssignResultVar2Args(
                             CondLinConLE{ { {{1.0}, {compl_var}}, var_lb } }),
                GetMC().AssignResultVar2Args(
                             CondLinConGE{ { {{1.0}, {expr_var}},  con_rhs } })
                           } });
      /// res2 = (body==0)
      auto res2 = GetMC().AssignResultVar2Args(
            CondLinConEQ{ { {{1.0}, {expr_var}}, con_rhs } });
      /// res3 = (var >= ub && con <= 0)
      auto res3 = GetMC().AssignResultVar2Args(
            AndConstraint{ {
                GetMC().AssignResultVar2Args(
                             CondLinConGE{ { {{1.0}, {compl_var}}, var_ub } }),
                GetMC().AssignResultVar2Args(
                             CondLinConLE{ { {{1.0}, {expr_var}}, con_rhs } })
                           } });
      /// res4 = (res1 \/ res2 \/ res3)
      auto res4 = GetMC().AssignResultVar2Args(
            OrConstraint{ { res1, res2, res3 } });
      GetMC().FixAsTrue(res4);
      /// Not adding any static algebraic constraint
    }
  }

  /// Reuse the stored ModelConverter
  using Base::GetMC;
};


/// Typedef linear compl cvt
template <class MC>
using ComplCvtLin_MIP = ComplementarityConverter_MIP<MC, ComplementarityLinear>;


/// Typedef quadratic compl cvt
template <class MC>
using ComplCvtQuad_MIP = ComplementarityConverter_MIP<MC, ComplementarityQuadratic>;

} // namespace mp

#endif // COMPLEMENT_H
