#ifndef EXPR_ALG_INLINE_H
#define EXPR_ALG_INLINE_H

#include <array>
#include <type_traits>

#include "mp/flat/constr_std.h"
#include "mp/flat/bucketaccum.h"
#include "mp/flat/eexpr.h"

namespace mp {

/// A mix-in class inlining algebraic expressions
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
    if (nlfc && (fLin || fQuad) || nqfc && fQuad) {
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
  void WalkObjectives() {}

  /// @return true iff made this one redundant
  template <class Body, class RangeOrRHS>
  bool InlineAlgExpr(
      const AlgebraicConstraint<Body, RangeOrRHS>& con, int i) {
    if (HasAlgExpr(con)) {
      auto qexpr = CollectAlgSubExpr(con.GetBody().GetLinTerms());
      if constexpr (std::is_same_v<Body, LinTerms>) {  // a LinCon..
        assert(
            qexpr.GetLinTerms()!=con.GetBody().GetLinTerms()
            || qexpr.GetQPTerms().size());
      } else {                                         // a QuadCon..
        assert(qexpr.GetBody()!=con.GetBody());  // Do we receive the final body?
      }
      return true;
    }
    return false;
  }

  /// @return recursively collect linear and/or quadratic
  /// subexpressions
  QuadraticExpr CollectAlgSubExpr(const LinTerms& lt) {
    BucketAccumulator<EExpr> buckets;
    AffineExpr ae_untouched;

    buckets.Add(ae_untouched);
    return buckets.ExtractSum();
  }

  /// Whether the alg con has alg subexpressions
  template <class Body, class RangeOrRHS>
  bool HasAlgExpr(
      const AlgebraicConstraint<Body, RangeOrRHS>& con) {
    const auto& lt = con.GetBody().GetLinTerms();
    for (auto i=lt.size(); i--; ) {
      auto vi = lt.var(i);
      if (auto pLFC = MPCD(
              template GetActiveInitExpressionOfType<
                  LinearFunctionalConstraint>(vi) ))
        return true;
      if (auto pQFC = MPCD(
              template GetActiveInitExpressionOfType<
                  QuadraticFunctionalConstraint>(vi) )) {
        if (pQFC->GetArguments().GetLinTerms().size()
            || fQuad_)
          return true;
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
