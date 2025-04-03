#ifndef PREPRO_PROD_H
#define PREPRO_PROD_H

#include <tuple>
#include <algorithm>
#include <cassert>

#include "mp/common.h"
#include "mp/expr-visitor.h"
#include "mp/flat/constr_std.h"

#include "mp/utils-vec.h"


namespace mp {

/// Preprocess and flatten product of NL expressions.
///
/// Scheme:
/// 1. Collect an array of terms (flattened.)
/// 2. Group binary terms into a FORALL. If desired.
/// 3. Quadratize remaining terms.
template <class Flattener>
class PreproProd
    : public ExprConverter<PreproProd<Flattener>, void> {
public:
  /// Construct
  PreproProd(Flattener& flt) : flt_(flt) { }

  /// typedef flat expr
  using FlatExpr = typename Flattener::EExprType;

  /// Term comparator.
  /// .first <==> category (const, bin/negbin, int, real)
  /// .second - bounds range
  using TermCmp = std::pair<int, double>;

  /// Flatten product
  FlatExpr FlattenProduct(BinaryExpr nl_e) {
    CollectFactors(nl_e);
    return CombineFactors();
  }


protected:
  /// Collect product factors from a binary expression tree
  void CollectFactors(BinaryExpr nl_e) {
    assert(expr::MUL == nl_e.kind());
    this->Visit(nl_e);
  }

public:
  /// Visit a multiplication node
  void VisitMul(BinaryExpr e) {
    this->Visit(e.lhs());
    this->Visit(e.rhs());
  }

  /// Visit any other node
  void VisitUnsupported(Expr e) {
    auto t =
      std::make_tuple< TermCmp, FlatExpr, std::pair<double, double> >(
        {}, GetFlt().Visit(e), {}
      );
    terms_flt_.push_back( std::move(t) );
  }

protected:
  /// Quadratize / logicalize (sub)products
  FlatExpr CombineFactors() {
    SortTerms();
    return CombineOrderedFactors();
  }

  void SortTerms() {
    for (auto& tpl: terms_flt_) {
      auto bnds =
          GetFlt().GetFlatCvt().ComputeBoundsAndType(
          std::get<1>(tpl));
      bool is_int = var::INTEGER==bnds.type();
      assert(bnds.lb() <= bnds.ub());
      bool is_const = bnds.ub()<=bnds.lb();
      n_terms_const_ += is_const;
      bool is_bin_or_neg_bin =
          is_int &&
          ((!bnds.lb() && 1.0==bnds.ub()) ||
           (-1.0==bnds.lb() && !bnds.ub()) );
      n_terms_binary_ += is_bin_or_neg_bin;
      int category = is_const ? 0 :
          is_bin_or_neg_bin ? 1 :
          is_int ? 2 : 3;
      std::get<0>(tpl) = {category, bnds.ub()-bnds.lb() };
      std::get<2>(tpl) =
          {bnds.lb(), bnds.ub()};
    }
    std::stable_sort(terms_flt_.begin(), terms_flt_.end(),
              [](const auto& tpl1, const auto& tpl2) {
      return std::get<0>(tpl1) < std::get<0>(tpl2);
    });
  }

  FlatExpr CombineOrderedFactors() {
    FlatExpr result = typename FlatExpr::Constant {1.0};
    int i=0;
    double coef00 = 1.0;

    for ( ; i<n_terms_const_; ++i) {      // Collect constant factors
      // no: can have constant args
      // assert(std::get<1>(terms_flt_[i]).is_constant());
      auto bnds = std::get<2>(terms_flt_[i]);
      assert(bnds.first == bnds.second);
      coef00 *= bnds.first;
    }

    if ((n_terms_binary_==2 && GetFlt().prepro_products()&2) ||
        (n_terms_binary_>=3 && GetFlt().prepro_products()&4)) {
      AndConstraint::Arguments args_forall;
      args_forall.reserve(n_terms_binary_);   // Logicalize binary product
      for ( ; i<n_terms_const_+n_terms_binary_; ++i) {
        assert(1.0 ==                 // term bounds range 1
               std::get<2>(terms_flt_[i]).second -
                   std::get<2>(terms_flt_[i]).first );
        int binvar = GetFlt().Convert2Var(std::move(std::get<1>(terms_flt_[i])));
        auto is_lb_minus1 = (-1.0 == std::get<2>(terms_flt_[i]).first);
        if (is_lb_minus1) {    // negated binary
          coef00 *= -1;
          binvar = GetFlt().Convert2Var( { {-1.0}, {binvar} } );
        } else {
          assert(0.0 == std::get<2>(terms_flt_[i]).first); // normal binary
        }
        args_forall.push_back( binvar );
      }
      result = GetFlt().AssignResult2Args( AndConstraint{args_forall} );
      // Context should be set when adding the top contraint
    }

    result *= coef00;

    for ( ; i<(int)terms_flt_.size(); ++i) {
      result = GetFlt().QuadratizeOrLinearize(
          result, std::get<1>(terms_flt_[i]));
    }

    return result;
  }

  /// Obtain flattener, const
  const Flattener& GetFlt() const { return flt_; }
  /// Obtain flattener
  Flattener& GetFlt() { return flt_; }


private:
  Flattener& flt_;
  /// tuple: comparator, term, bounds
  SmallVec<
      std::tuple< TermCmp, FlatExpr, std::pair<double, double> >,
      32 >                 // 32 elements preallocated
      terms_flt_;
  int n_terms_const_ = 0;
  int n_terms_binary_ = 0;
};

 }  // namespace mp

#endif // PREPRO_PROD_H
