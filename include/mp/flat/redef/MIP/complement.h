#ifndef COMPLEMENT_H
#define COMPLEMENT_H

#include <cmath>

#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constraints_std.h"

namespace mp {

/// Converts Complementatrity for MIP
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

    using AlgConBodyType = typename
        std::decay< decltype (expr.GetAlgConBody()) >::type;
    using AlgConLE = AlgebraicConstraint<
        AlgConBodyType, AlgConRhs<-1> >;
    using AlgConEQ = AlgebraicConstraint<
        AlgConBodyType, AlgConRhs<0> >;
    using AlgConGE = AlgebraicConstraint<
        AlgConBodyType, AlgConRhs<1> >;

    using CondConLE = ConditionalConstraint< AlgConLE >;
    using CondConEQ = ConditionalConstraint< AlgConEQ >;

    if (fin_var_lb && !fin_var_ub) {
      /// res1 = (var <= var_lb)
      auto res_neg_var_lb = GetMC().AssignResultVar2Args(
            CondLinConLE{ { {{1.0}, {compl_var}}, var_lb } });
      /// res2 = (body <= lb)
      auto res_neg_con_lb = GetMC().AssignResultVar2Args(
            CondConLE{ { expr.GetAlgConBody(), -expr.constant_term() } });
      /// res3 = (res1 \/ res2)
      auto res_disj = GetMC().AssignResultVar2Args(
            OrConstraint{ { res_neg_var_lb, res_neg_con_lb } });
      GetMC().FixAsTrue(res_disj);
      /// Add the algebraic constraint
      GetMC().AddConstraint(
            AlgConGE{ expr.GetAlgConBody(), -expr.constant_term() } );
    } else if (fin_var_ub && !fin_var_lb) {
      /// res1 = (var >= var_ub)
      auto res_neg_var_ub = GetMC().AssignResultVar2Args(
            CondLinConLE{ { {{-1.0}, {compl_var}}, -var_ub } });
      /// res2 = (body >= ub)
      auto expr = cc.GetExpression();                         // copy
      expr.negate();
      auto res_neg_con_ub = GetMC().AssignResultVar2Args(
            CondConLE{ { expr.GetAlgConBody(), -expr.constant_term() } });
      /// res3 = (res1 \/ res2)
      auto res_disj = GetMC().AssignResultVar2Args(
            OrConstraint{ { res_neg_var_ub, res_neg_con_ub } });
      GetMC().FixAsTrue(res_disj);
      /// Add the algebraic constraint
      GetMC().AddConstraint(
            AlgConGE{ expr.GetAlgConBody(), -expr.constant_term() } );
    } else {
      assert(fin_var_lb && fin_var_ub);
      /// res1 = (var <= lb && con >= 0)
      auto expr_neg = expr;
      expr_neg.negate();
      auto res1 = GetMC().AssignResultVar2Args(
            AndConstraint{ {
                GetMC().AssignResultVar2Args(
                             CondLinConLE{ { {{1.0}, {compl_var}}, var_lb } }),
                GetMC().AssignResultVar2Args(        // TODO wrong for QuadCon
                             CondConLE{ { expr_neg.GetAlgConBody(),
                                             -expr_neg.constant_term() } })
                           } });
      /// res2 = (body==0)
      auto res2 = GetMC().AssignResultVar2Args(      // TODO wrong for QuadCon
            CondConEQ{ { expr.GetAlgConBody(), -expr.constant_term() } });
      /// res3 = (var >= ub && con <= 0)
      auto res3 = GetMC().AssignResultVar2Args(
            AndConstraint{ {
                GetMC().AssignResultVar2Args(
                             CondLinConLE{ { {{-1.0}, {compl_var}}, -var_ub } }),
                GetMC().AssignResultVar2Args(        // TODO wrong for QuadCon
                             CondConLE{ { expr.GetAlgConBody(),
                                             -expr.constant_term() } })
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
