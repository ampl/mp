#ifndef COMPLEMENT_H
#define COMPLEMENT_H

#include <cmath>

#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constr_std.h"
#include "mp/flat/nl_expr/constr_nl.h"

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

  /// Reuse the stored ModelConverter
  using Base::GetMC;

  /// Skip conversion?
  bool IfDelayConversion(const ItemType& , int ) {
    return
        GetMC().IfWantNLOutput()
        && GetMC().UserAcceptsAndRecommends(
               (const NLComplementarity*)nullptr);
  }

  /// Convert in any context
  Context Convert(const ItemType& cc, int i) {
    switch (GetMC().ComplementarityCvt()) {
    case 0:
      return Convert2Disj(cc, i);
    case 1:
      return Convert2Prod(cc, i);
    case 2:
      return Convert2FB(cc, i);
    case 3:
      return Convert2Min(cc, i);
    default:
      MP_RAISE("Wrong value for cvt:compl");
    }
  }

protected:
  /// Convert into disjunction
  Context Convert2Disj(const ItemType& cc, int ) {
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
        AlgebraicExpression<typename ComplCon::ExprType::BodyType>
        {expr.GetBody(), 0.0} );
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
    return Context::CTX_ROOT;
  }

  /// Convert into product - which solvers do this well?
  /// Baron, according to its docu. But should be var1*var2? No.
  Context Convert2Prod(const ItemType& cc, int ) {
    const auto& expr = cc.GetExpression();
    auto compl_var = cc.GetVariable();

    double
        var_lb = GetMC().lb(compl_var),
        var_ub = GetMC().ub(compl_var);

    bool
        fin_var_lb = std::isfinite(var_lb),
        fin_var_ub = std::isfinite(var_ub);

    /// Using algebraic expression (expr 1:1 here)
    auto expr_var = GetMC().Convert2Var( expr );

    if (fin_var_lb && !fin_var_ub) {
      // (v-lb)*expr = 0
      GetMC().AddConstraint(QuadConEQ{ { {{-var_lb}, {expr_var}},
                                       {{1.0}, {compl_var}, {expr_var} }},
                                      GetMC().ComplementarityCvtTol() });
      // Add the algebraic constraint via the representing variable
      // Propagate mixed context (logical constraint would set CTX_NEG)
      GetMC().set_var_lb_context(expr_var, 0.0, Context::CTX_MIX);
    } else if (fin_var_ub && !fin_var_lb) {
      // (v-ub)*expr = 0
      GetMC().AddConstraint(QuadConEQ{ { {{-var_ub}, {expr_var}},
                                       {{1.0}, {compl_var}, {expr_var} }},
                                      GetMC().ComplementarityCvtTol() });
      // Add the algebraic constraint via the representing variable
      // Propagate mixed context (logical constraint would set CTX_POS)
      GetMC().set_var_ub_context(expr_var, 0.0, Context::CTX_MIX);
    } else {
      assert(fin_var_lb && fin_var_ub);
      // Reduce mixed compl into 2 standard ones
      // le <= expr, le <= 0
      int le = (int)GetMC().AddVar(GetMC().MinusInfty(), 0.0);
      GetMC().AddConstraint(LinConLE{{ {1.0, -1.0}, {le, expr_var} }, 0.0});
      // ue >= expr, ue >= 0
      int ue = (int)GetMC().AddVar(0.0, GetMC().Infty());
      GetMC().AddConstraint(LinConGE{{ {1.0, -1.0}, {ue, expr_var} }, 0.0});
      // (v-lb)*ue = 0
      GetMC().AddConstraint(QuadConEQ{ { {{-var_lb}, {ue}},
                                       {{1.0}, {compl_var}, {ue} }},
                                      GetMC().ComplementarityCvtTol() });
      // (v-ub)*le = 0
      GetMC().AddConstraint(QuadConEQ{ { {{-var_ub}, {le}},
                                       {{1.0}, {compl_var}, {le} }},
                                      GetMC().ComplementarityCvtTol() });
    }
    return Context::CTX_ROOT;
  }


  /// Convert into Fischer-Burmeister function
  Context Convert2FB(const ItemType& cc, int i) {
    auto FBfn =               // @todo NormConstraint could be simpler
        [this](double sv, int v, double bnd, double se, int e, double eps) {
          auto root_from = this->GetMC().AssignResultVar2Args(
              QuadraticFunctionalConstraint( {
                  { { {2*sv*bnd}, {v} },       // lin part
                   { {1.0, 1.0},               // QP part
                    {v, e}, {v, e} } },
                  bnd*bnd + 2*eps } ) );
          auto root = this->GetMC().AssignResultVar2Args(
              PowConstExpConstraint( {root_from}, DblParamArray1{0.5} ) );
          return
              LinConEQ( { {1.0, -sv, -se}, {root, v, e} }, bnd );
        };
    return ConvertWithNCPCon(cc, i, FBfn);
  }


  /// Convert into min(., .)
  Context Convert2Min(const ItemType& cc, int i) {
    auto MINfn =
        [this](double sv, int v, double bnd, double se, int e, double eps) {
          auto a = this->GetMC().AssignResultVar2Args(
              LinearFunctionalConstraint(
                  { { {sv}, {v} }, bnd } ) );
          auto b = this->GetMC().AssignResultVar2Args(
              LinearFunctionalConstraint(
                  { { {se}, {e} }, 0.0 } ) );
          auto minab = this->GetMC().AssignResultVar2Args(
              MinConstraint( {a, b}, {} ) );
          return              // @todo could just set bounds minab==0
              LinConEQ( { {1.0}, {minab} }, 0.0 );
        };
    return ConvertWithNCPCon(cc, i, MINfn);
  }

  /// Convert with a given NCP constraint.
  /// NCPcon(sv, v, bnd, se, e, eps):
  ///   e.g., FB(sv*v+bnd, se*e, eps) == 0
  template <class NCPConFn>
  Context ConvertWithNCPCon(const ItemType& cc, int i, NCPConFn NCPcon) {
    const auto& expr = cc.GetExpression();
    auto compl_var = cc.GetVariable();

    double
        var_lb = GetMC().lb(compl_var),
        var_ub = GetMC().ub(compl_var);

    bool
        fin_var_lb = std::isfinite(var_lb),
        fin_var_ub = std::isfinite(var_ub);

    /// Using algebraic expression (expr 1:1 here)
    auto expr_var = GetMC().Convert2Var( expr );

    if (fin_var_lb && !fin_var_ub) {
      // NCP(v-lb)*expr = 0
      GetMC().AddConstraint(
          NCPcon(1.0, compl_var, -var_lb, 1.0, expr_var,
                GetMC().ComplementarityCvtTol() ) );
      // Add the algebraic constraint via the representing variable.
      // Actually not necessary with a NCP function.
      // Propagate mixed context (logical constraint would set CTX_NEG)
      GetMC().set_var_lb_context(expr_var, 0.0, Context::CTX_MIX);
    } else if (fin_var_ub && !fin_var_lb) {
      // (v-ub)*expr = 0
      GetMC().AddConstraint(
          NCPcon(-1.0, compl_var, var_ub, -1.0, expr_var,
                GetMC().ComplementarityCvtTol() ) );
      // Add the algebraic constraint via the representing variable
      // Actually not necessary with a NCP function.
      // Propagate mixed context (logical constraint would set CTX_POS)
      GetMC().set_var_ub_context(expr_var, 0.0, Context::CTX_MIX);
    } else {
      assert(fin_var_lb && fin_var_ub);
      // Reduce mixed compl into 2 standard ones
      // le <= expr, le <= 0
      int le = (int)GetMC().AddVar(GetMC().MinusInfty(), 0.0);
      GetMC().AddConstraint(LinConLE{{ {1.0, -1.0}, {le, expr_var} }, 0.0});
      // ue >= expr, ue >= 0
      int ue = (int)GetMC().AddVar(0.0, GetMC().Infty());
      GetMC().AddConstraint(LinConGE{{ {1.0, -1.0}, {ue, expr_var} }, 0.0});
      // (v-lb)*ue = 0
      GetMC().AddConstraint(
          NCPcon(1.0, compl_var, -var_lb, 1.0, ue,
                 GetMC().ComplementarityCvtTol() ) );
      // (ub-v)*(-le) = 0
      GetMC().AddConstraint(
          NCPcon(-1.0, compl_var, var_ub, -1.0, le,
                 GetMC().ComplementarityCvtTol() ) );
    }
    return Context::CTX_ROOT;
  }
};


/// Typedef linear compl cvt
template <class MC>
using ComplCvtLin_MIP = ComplementarityConverter_MIP<MC, ComplementarityLinear>;


/// Typedef quadratic compl cvt
template <class MC>
using ComplCvtQuad_MIP = ComplementarityConverter_MIP<MC, ComplementarityQuadratic>;

} // namespace mp

#endif // COMPLEMENT_H
