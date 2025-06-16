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
      bool fLin, bool fQuad) {
    auto nlfc = MPCD( GetNumberOfAddable(
        (LinearFunctionalConstraint*)nullptr) );
    auto nqfc = MPCD( GetNumberOfAddable(
        (QuadraticFunctionalConstraint*)nullptr) );
    if ((nlfc && (fLin || fQuad)) || (nqfc && fQuad)) {
      fLin_ = fLin;
      fQuad_ = fQuad;
      return DoConsiderInlining();
    }
    return false;
  }

protected:
  /// @todo Conditionals, indicators?
  /// But any new ones would be linearized before this action.
  bool DoConsiderInlining() {
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
    for (auto& obj: MPD( get_objectives() )) {
      if (HasAlgExpr(obj.GetLinTerms())) {
        auto qexpr
            = CollectAlgSubExpr(obj.GetLinTerms(), obj.GetQPTerms());
        // @todo some linking for this...
        // but we modify in-place
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
                 { qexpr.GetLinTerms(), qexpr.GetQPTerms() },
                 range_or_rhs
               }) );
      else
        MPD( AddConstraint(
               AlgebraicConstraint<LinTerms, RangeOrRHS>{
                 { qexpr.GetLinTerms() },
                 range_or_rhs
               }) );
      return true;
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

  /// @return recursively collect linear and/or quadratic
  /// subexpressions
  /// @todo Consider \a fQuad_
  EExpr CollectAlgSubExpr(
      const LinTerms& lt0, const QuadTerms& qt0) {
    BucketAccumulator<EExpr> buckets;
    AffineExpr ae_untouched;     // unmodified linear terms

    auto inline_alg_subexpr
        = [&](const auto& subexpr, double ci, int vi) {
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
      if (auto pLFC = MPCD(
              template GetActiveInitExpressionOfType<
                  LinearFunctionalConstraint>(vi) ))
        inline_alg_subexpr(pLFC->GetAffineExpr(), ci, vi);
      else if (auto pQFC = MPCD(
              template GetActiveInitExpressionOfType<
                  QuadraticFunctionalConstraint>(vi) ))
        inline_alg_subexpr(pQFC->GetArguments(), ci, vi);
      else
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
  /// has no stronger bounds than those from the expression.
  /// @todo Possibly we should cache this, as otherwise
  /// we'd take as expression some linear terms later
  /// if the expression starts implying tighter result
  bool HasAlgExpr(const LinTerms& lt) {
    for (auto i=lt.size(); i--; ) {
      auto vi = lt.var(i);
      if (auto pLFC = MPCD(
              template GetActiveInitExpressionOfType<
                  LinearFunctionalConstraint>(vi) ))
        return !MPD( IfVarBoundsStrongerThanInitExpr(vi) );
      if (auto pQFC = MPCD(
              template GetActiveInitExpressionOfType<
                  QuadraticFunctionalConstraint>(vi) )) {
        if (pQFC->GetArguments().GetLinTerms().size()
            || fQuad_)
          return !MPD( IfVarBoundsStrongerThanInitExpr(vi) );
      }
    }
    return false;
  }


private:
  bool fLin_ {}, fQuad_ {};
  std::array<int, 8> last_con_ {};   // does this init to 0's?
};

}  // namespace mp

#endif // EXPR_ALG_INLINE_H
