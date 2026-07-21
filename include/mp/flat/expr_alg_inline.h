#ifndef EXPR_ALG_INLINE_H
#define EXPR_ALG_INLINE_H

#include <array>
#include <type_traits>

#include "mp/flat/constr_std.h"
#include "mp/flat/bucketaccum.h"
#include "mp/flat/eexpr.h"
#include "mp/valcvt-link.h"

namespace mp {

/// A mix-in class inlining algebraic expressions #266,
/// represented as functional flat constraints
/// (LFC and QFC).
template <class Impl>
class AlgebraicExpressionInliner {
public:
  /// Consider inlining if the corresponding type
  /// of expressions exist and desired for that.
  bool ConsiderInliningAlgExpr(
      int flags) {
    bool result {};
    fLin_ = flags & 2;
    fQuad_ = flags & 4;

    auto nlfc = MPCD( GetNumberOfAddable(
        (LinearFunctionalConstraint*)nullptr) );
    auto nqfc = MPCD( GetNumberOfAddable(
        (QuadraticFunctionalConstraint*)nullptr) );
    if ((nlfc && fLin_) || (nqfc && fQuad_)) {
      result = result || DoConsiderInliningStaticAlgCons();
    }

    if (flags & 16) {
      nlfc = MPCD( GetNumberOfAddable(
          (LinearFunctionalConstraint*)nullptr) );
      nqfc = MPCD( GetNumberOfAddable(
          (QuadraticFunctionalConstraint*)nullptr) );
      if ((nlfc && fLin_) || (nqfc && fQuad_)) {
        result = result || DoConsiderInliningInIndicators();
      }
    }
    return result;
  }

  /// Check if the variable can be eliminated.
  /// For that, it should not have stronger submitted bounds,
  /// and not marked as explicit, e.g., by dvelim=0.
  bool CanBeEliminated(int var) const {
    return
        CanBeEliminated_FastCheck(var)
        &&
        !MPCD( IfSubmittedVarBoundsStrongerThanInitExpr(var) );
  }

  /// Check if the variable can be eliminated.
  /// Only the fast check:
  /// not marked as explicit, e.g., by dvelim=0.
  bool CanBeEliminated_FastCheck(int var) const
  { return !MPCD(IsExplicitDV(var)); }


protected:
  /// @todo Currently always return false...
  bool DoConsiderInliningStaticAlgCons() {
    Walk2Inline<LinConRange>(0);
    Walk2Inline<LinConLE>(1);
    Walk2Inline<LinConEQ>(2);
    Walk2Inline<LinConGE>(3);

    Walk2Inline<QuadConRange>(4);
    Walk2Inline<QuadConLE>(5);
    Walk2Inline<QuadConEQ>(6);
    Walk2Inline<QuadConGE>(7);

    WalkObjectives();

    return false;
  }

  /// Any new ones would have been linearized before this action.
  /// @todo Currently always return false...
  bool DoConsiderInliningInIndicators() {
    Walk2Inline<IndicatorConstraintLinLE>(8);
    Walk2Inline<IndicatorConstraintLinEQ>(9);
    Walk2Inline<IndicatorConstraintLinGE>(10);

    Walk2Inline<IndicatorConstraintQuadLE>(11);
    Walk2Inline<IndicatorConstraintQuadEQ>(12);
    Walk2Inline<IndicatorConstraintQuadGE>(13);

    return false;
  }

  /// @param iType: con type index
  template <class Con>
  void Walk2Inline(size_t iType) {
    auto& ck = MPD( GetConstraintKeeper( (Con*)nullptr) );
    ck.ForEachActive(
        [this](const Con& con, int i) {
          return MPD( InlineAlgExpr(con, i) );
        }, last_con_.at(iType), ck.Size());
    last_con_[iType] = ck.Size();  // @todo save across types
  }

  /// @todo just once?
  void WalkObjectives() {
    auto& objectives = MPD( get_objectives() );
    for (size_t iobj=0; iobj < objectives.size(); ++iobj) {
      auto& obj = objectives[iobj];
      if (HasAlgExpr(obj.GetLinTerms())) {
        auto qexpr
            = CollectAlgSubExpr(obj.GetLinTerms(), obj.GetQPTerms());
        auto obj_src = MPD( GetObjValueSourceNode() ).Select(iobj);
        pre::AutoLinkScope auto_link_scope
            { *(Impl*)this, obj_src };
        MPD( UncountArgRefs(obj) );
        obj.GetLinTerms() = std::move(qexpr.GetBody().GetLinTerms());
        obj.GetQPTerms() = std::move(qexpr.GetBody().GetQPTerms());
        if (qexpr.constant_term()) {
          obj.GetLinTerms().add_term(
                qexpr.constant_term(),
                int( MPD( MakeFixedVar(1.0) ) ) );
          obj.GetLinTerms().sort_terms();     // @todo merge would be faster?
        }
        MPD( CountArgRefs(obj) );
      }
    }
  }

  /// Inline algebraic subexpressions
  /// into an algebraic constraint.
  /// @return true iff made this one redundant
  template <class Body, class RangeOrRHS>
  bool InlineAlgExpr(
      const AlgebraicConstraint<Body, RangeOrRHS>& con, int i) {
    if (HasAlgExpr(con.GetBody().GetLinTerms())) {
      auto qexpr = CollectAlgSubExpr(con.GetBody());
      auto range_or_rhs = con.GetRhsOrRange();
      range_or_rhs.add_to_rhs( -qexpr.constant_term() );   // subtract
      auto auto_link_scope = MPD( MakeAutoLinker(con, i) );
      if (qexpr.GetQPTerms().size())
        MPD( AddConstraint(        // not _AS_ROOT
            AlgebraicConstraint<QuadAndLinTerms, RangeOrRHS>{
            { std::move(qexpr.GetLinTerms()),
             std::move(qexpr.GetQPTerms()) },
            range_or_rhs
            }) );
      else
        MPD( AddConstraint(
            AlgebraicConstraint<LinTerms, RangeOrRHS>{
            { std::move(qexpr.GetLinTerms()) },
            range_or_rhs
            }) );
      return true;
    }
    return false;
  }

  /// Inline algebraic subexpressions
  /// into an indicator constraint.
  /// @return true iff made this one redundant
  template <class Body, class RangeOrRHS>
  bool InlineAlgExpr(
      const IndicatorConstraint<
              AlgebraicConstraint<Body, RangeOrRHS> >& con, int i) {
    const auto& algcon = con.get_constraint();
    if (HasAlgExpr(algcon.GetBody().GetLinTerms())) {
      auto qexpr = CollectAlgSubExpr(algcon.GetBody());
      auto range_or_rhs = algcon.GetRhsOrRange();
      range_or_rhs.add_to_rhs( -qexpr.constant_term() );   // subtract
      auto auto_link_scope = MPD( MakeAutoLinker(con, i) );
      if (qexpr.GetQPTerms().size()) {
        using IndConQP = IndicatorConstraint<
            AlgebraicConstraint<QuadAndLinTerms, RangeOrRHS> >;
        // Only if accepted, redefine into IndQuad:
        if (MPCD( UserAcceptsAndRecommends((const IndConQP*)0) )) {
          MPD( AddConstraint(        // not _AS_ROOT
              IndConQP
              { con.get_binary_var(), con.get_binary_value(),
               { { std::move(qexpr.GetLinTerms()),
                std::move(qexpr.GetQPTerms()) },
               range_or_rhs
               } }) );
          return true;
        }
      } else {
        MPD( AddConstraint(
            IndicatorConstraint<
                AlgebraicConstraint<LinTerms, RangeOrRHS> >
            { con.get_binary_var(), con.get_binary_value(),
             { { std::move(qexpr.GetLinTerms()) },
             range_or_rhs
             } }) );
        return true;
      }
    }
    return false;
  }

  /// @return recursively collect linear and/or quadratic
  /// subexpressions
  /// @todo Consider \a fQuad_
  template <class Body>
  EExpr CollectAlgSubExpr(const Body& body) {
    if constexpr (std::is_same_v<Body, QuadAndLinTerms>) {
      return CollectAlgSubExpr(
            body.GetLinTerms(), body.GetQPTerms());
    }
    return CollectAlgSubExpr(body.GetLinTerms(), {});
  }

  /// Recursively collect linear and/or quadratic
  /// subexpressions
  /// @return The result
  /// @todo Consider \a fQuad_
  /// @todo Join with HasAlgExpr()
  EExpr CollectAlgSubExpr(
      const LinTerms& lt0, const QuadTerms& qt0) {
    BucketAccumulator<EExpr> buckets;
    AffineExpr ae_untouched;     // unmodified linear terms

    auto inline_alg_subexpr
        = [&](const auto& subexpr, double ci, int /*vi*/) {
      EExpr collected;
      if (HasAlgExpr(subexpr.GetBody().GetLinTerms())) {
        collected = CollectAlgSubExpr(
              subexpr.GetBody());
        collected.add_to_constant(subexpr.constant_term());
      } else {
        collected = EExpr{subexpr};
      }
      collected *= ci;
      buckets.Add(std::move(collected));
      // No: it is done when removing the old top-level constraint:
      // MPD( DecrementVarUsage(vi) );
    };

    for (auto i=lt0.size(); i--; ) {
      auto ci = lt0.coef(i);
      auto vi = lt0.var(i);
      const LinearFunctionalConstraint* pLFC;
      const QuadraticFunctionalConstraint* pQFC;
      if ((pLFC = MPCD(
               template GetActiveInitExpressionOfType<
                   LinearFunctionalConstraint>(vi) ))
          && CanBeEliminated(vi)) {
        inline_alg_subexpr(pLFC->GetAffineExpr(), ci, vi);
      } else if ((pQFC = MPCD(
                      template GetActiveInitExpressionOfType<
                          QuadraticFunctionalConstraint>(vi) ))
                 && CanBeEliminated(vi)) {
        inline_alg_subexpr(pQFC->GetArguments(), ci, vi);
      } else
        ae_untouched.add_term(ci, vi);
    }

    // untouched as well: reuse all QP terms
    EExpr ee_untouched
    { std::move(ae_untouched), qt0, 0.0};
    buckets.Add(ee_untouched);
    auto result = buckets.ExtractSum();
    assert(
          result.GetLinTerms()!=lt0
        || result.GetQPTerms()!=qt0);
    return result;
  }

  /// Whether the alg con/expr body has alg subexpressions.
  /// We only consider subexpressions whose result var
  /// has no stronger submitted bounds
  /// than those from the subexpression.
  /// @todo Possibly we should cache this, as otherwise
  /// we'd take as expression some linear terms later
  /// if the expression starts implying tighter result
  bool HasAlgExpr(const LinTerms& lt) {
    for (auto i=lt.size(); i--; ) {
      auto vi = lt.var(i);
      if (auto pLFC = MPCD(
              template GetActiveInitExpressionOfType<
                  LinearFunctionalConstraint>(vi) ))
        if (CanBeEliminated(vi))
          return true;
      if (auto pQFC = MPCD(
              template GetActiveInitExpressionOfType<
                  QuadraticFunctionalConstraint>(vi) )) {
        if (pQFC->GetArguments().GetLinTerms().size()
            || fQuad_)
          if (CanBeEliminated(vi))
            return true;
      }
    }
    return false;
  }


private:
  bool fLin_ {}, fQuad_ {};
  std::array<int, 14> last_con_ {};   // does this init to 0's?
};

}  // namespace mp

#endif // EXPR_ALG_INLINE_H
