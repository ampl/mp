#ifndef PREPRO_PROD_H
#define PREPRO_PROD_H

#include <cmath>
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
    isRemoved_.resize(terms_flt_.size());
    SortTerms();
    return CombineOrderedFactors();
  }

  void SortTerms() {
    for (auto& tpl: terms_flt_) {
      auto bnds =
          GetFlt().GetFlatCvt().ComputeBoundsAndType(
          std::get<1>(tpl));
      auto pic = GetFlt().GetFlatCvt().ClassifyPreproInfo(bnds);
      n_terms_const_ += pic.is_const_;
      n_terms_binary_ += pic.is_bin_or_neg_bin_;
      int category = pic.is_const_ ? 0 :
          pic.is_bin_or_neg_bin_ ? 1 :
          pic.is_int_ ? 2 : 3;
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
    CombineConstants();
    CombineBinaryFactors();
    RecognizeSignpows();
    CombineRemainingFactors();

    return result;
  }

  void CombineConstants() {
    for ( ; iTerm<n_terms_const_; ++iTerm) {      // Collect constant factors
      // no: can have constant args
      // assert(std::get<1>(terms_flt_[i]).is_constant());
      auto bnds = std::get<2>(terms_flt_[iTerm]);
      assert(bnds.first == bnds.second);
      coef00 *= bnds.first;
    }
  }

  void CombineBinaryFactors() {
    if ((n_terms_binary_==2 && GetFlt().prepro_products()&2) ||
        (n_terms_binary_>=3 && GetFlt().prepro_products()&4)) {
      AndConstraint::Arguments args_forall;
      args_forall.reserve(n_terms_binary_);   // Logicalize binary product
      for ( ; iTerm<n_terms_const_+n_terms_binary_; ++iTerm) {
        assert(1.0 ==                 // term bounds range 1
               std::get<2>(terms_flt_[iTerm]).second -
                   std::get<2>(terms_flt_[iTerm]).first );
        int binvar = GetFlt().Convert2Var(std::move(std::get<1>(terms_flt_[iTerm])));
        auto is_lb_minus1 = (-1.0 == std::get<2>(terms_flt_[iTerm]).first);
        if (is_lb_minus1) {    // negated binary
          coef00 *= -1;
          binvar = GetFlt().Convert2Var( { {-1.0}, {binvar} } );
        } else {
          assert(0.0 == std::get<2>(terms_flt_[iTerm]).first); // normal binary
        }
        args_forall.push_back( binvar );
      }
      result = GetFlt().AssignResult2Args( AndConstraint{args_forall} );
      // Context should be set when adding the top contraint
    }
  }

  /// Information on a potential argument x
  /// of signpow(x, n)
  struct SignpowArgInfo {
    std::unordered_map<int, double> abspow_, justpow_;
    double totalabspow_ {}, totaljustpow_ {};
    bool invalid_ {};    // if cannot be
  public:
    bool IsValid() const { return !invalid_; }
    void MarkInvalid() { invalid_ = true; }
    void AddAbsPow(int iTerm, double p) {
      if (IsValid()) {
        auto res = abspow_.insert(std::pair{iTerm, p});
        assert(res.second);
        totalabspow_ += p;
      }
    }
    void AddJustPow(int iTerm, double p) {
      if (IsValid()) {
        auto res = justpow_.insert(std::pair{iTerm, p});
        assert(res.second);
        totaljustpow_ += p;
      }
    }
  };

  /// Tests: nonlinear/signpow_...
  void RecognizeSignpows() {
    if (GetFlt().recognize_signpow()) {
      // Pass 1: recognize abs(), sqrt(x^2k), possibly under pow()
      isSignpowAbs_.resize(terms_flt_.size());
      for (auto i1 = iTerm; i1<(int)terms_flt_.size(); ++i1) {
        if (CheckIfSignpowFactorAbsSqrtPossUnderPow(i1))
          isSignpowAbs_[i1] = true;
      }
      if (sparginfos_.size()) { // Pass 2: fitting simple factors
        for (auto i1 = iTerm; i1<(int)terms_flt_.size(); ++i1) {
          if (!isSignpowAbs_[i1])
            ConsiderSignpowFactorSimple(i1);
        }
        FinalizeSignpows();
      }
    }
  }

  bool CheckIfSignpowFactorAbsSqrtPossUnderPow(int iTerm) {
    const auto& term0 = std::get<1>(terms_flt_[iTerm]);
    if (term0.is_variable()) {   // Independent var, or a func con
      auto resvar0 = term0.get_representing_variable();
      if (const auto pConPow0 =      // pow(...)
          GetFlt().GetFlatCvt().template
          GetInitExpressionOfType<PowConstExpConstraint>(resvar0)) {
        auto pow0 = pConPow0->GetParameters()[0];
        if (pow0 > 0.0) {
          auto argvar0 = pConPow0->GetArguments()[0];
          return
              CheckIfSingpowArgAbsSqrt(iTerm, argvar0, pow0);
        }
      }
      return         // no pow()
          CheckIfSingpowArgAbsSqrt(iTerm, resvar0, 1.0);
    } else {                     // not a variable
      if (!term0.constant_term()   // simple case: abs(x)*abs(x)
          && term0.GetLinTerms().empty()  // can be out-mult'd from
          && 1==term0.GetQPTerms().size()) {     // abs(x)^2
        auto x = term0.GetQPTerms().var1(0);
        if (x == term0.GetQPTerms().var2(0)) {
          return
              CheckIfSingpowArgAbsSqrt(iTerm, x, 2.0);
        }
      }
    }
    return false;
  }

  void ConsiderSignpowFactorSimple(int iTerm) {
    const auto& term0 = std::get<1>(terms_flt_[iTerm]);
    if (term0.is_variable()) {   // Independent var, or a func con
      auto resvar0 = term0.get_representing_variable();
      if (const auto pConPow0 =      // pow(...)
          GetFlt().GetFlatCvt().template
          GetInitExpressionOfType<PowConstExpConstraint>(resvar0)) {
        auto argvar0 = pConPow0->GetArguments()[0];
        auto pow0 = pConPow0->GetParameters()[0];
        if (pow0<0.0
            || std::round(pow0)!=pow0) {
          MarkSignpowArgInvalid(iTerm, argvar0);
        } else {
          AddSignpowSimpleArg(iTerm, argvar0, pow0);
        }
      } else {                   // Indep var, or another func con
        AddSignpowSimpleArg(iTerm, resvar0, 1.0);
      }
    } else {                     // not a variable
      if (!term0.constant_term()   // simple case: x*x
          && term0.GetLinTerms().empty()  // can be out-mult'd from x^2
          && 1==term0.GetQPTerms().size()) {
        auto x = term0.GetQPTerms().var1(0);
        if (x == term0.GetQPTerms().var2(0)) {
          AddSignpowSimpleArg(iTerm, x, 2.0);
        }
      }
      // @todo handle (x+4), (x+4)^2 (would be x*x+8x+16),
      //   x*y+17 etc.
      // Find the expression among functional cons, otherwise skip.
      // Reason: to have a signpow(), we need at least an abs()
      // or sqrt(expr^2), and then expr would be a func con.
    }
  }

  bool CheckIfSingpowArgAbsSqrt(int iTerm, int resvar1, double pow1) {
    if (0.5 == pow1) { // sqrt = pow(..., 0.5) and we already in
      if (CheckIfSingpowSqrtArg(iTerm, resvar1, 1.0))
        return true;
    }
    if (const auto pConAbs1 =    // pow(abs(...), pow1)
        GetFlt().GetFlatCvt().template
        GetInitExpressionOfType<AbsConstraint>(resvar1)) {
      auto argvar1 = pConAbs1->GetArguments()[0];
      AddSignpowAbsPowArg(iTerm, argvar1, pow1);
      return true;
    }
    if (const auto pConSqrt1 =    // pow(sqrt(...), pow1)
        GetFlt().GetFlatCvt().template
        GetInitExpressionOfType<PowConstExpConstraint>(resvar1)) {
      auto resvar2 = pConSqrt1->GetArguments()[0];
      if (0.5 == pConSqrt1->GetParameters()[0]) {
        return CheckIfSingpowSqrtArg(iTerm, resvar2, pow1);
      }
    }
    return false;
  }

  /// This receives the argument of an sqrt()
  /// @param pow1 is the pow above sqrt()
  bool CheckIfSingpowSqrtArg(int iTerm, int resvar2, double pow1) {
    if (const auto pConQFC2 =          // pow(sqrt(x*x), pow1)
        GetFlt().GetFlatCvt().template // reason:
        GetInitExpressionOfType<QuadraticFunctionalConstraint>(
            resvar2)) {                // x^2 can be outmultiplied
      if (!pConQFC2->constant_term()   // @todo: (x+4)^2 etc.
          && pConQFC2->GetLinTerms().empty() // -recognize x^2+8x+16?
          && 1==pConQFC2->GetQPTerms().size()) {
        auto x = pConQFC2->GetQPTerms().var1(0);
        if (x == pConQFC2->GetQPTerms().var2(0)) {
          AddSignpowAbsPowArg(iTerm, x, pow1);
          return true;
        }
      }
    }
    if (const auto pConPow2 =    // pow(sqrt(pow(x, 2k)), pow1)
        GetFlt().GetFlatCvt().template
        GetInitExpressionOfType<PowConstExpConstraint>(resvar2)) {
      auto argvar2 = pConPow2->GetArguments()[0];
      auto pow3 = pConPow2->GetParameters()[0];
      auto pow3half = pow3 / 2.0;
      if (std::round(pow3half) == pow3half) {
        AddSignpowAbsPowArg(iTerm, argvar2, pow1*pow3half);
        return true;
      }
    }
    return false;
  }

  /// Register a [pow](abs/sqrt(pow^2k)) factor
  void AddSignpowAbsPowArg(int iTerm, int argvar, double powX) {
    sparginfos_[argvar].AddAbsPow(iTerm, powX);
    // @todo hash the argument LFC/QFC of argvar, if exists
  }

  /// Register a "simple" [pow](x) factor, x variable.
  /// Only if an abs(x) exists
  void AddSignpowSimpleArg(int iTerm, int argvar, double powN) {
    if (sparginfos_.end() != sparginfos_.find(argvar))
      sparginfos_[argvar].AddJustPow(iTerm, powN);
  }

  /// Mark invalid variable x for a signpow,
  /// due to x^2.7 etc.
  void MarkSignpowArgInvalid(int iTerm, int argvar) {
    if (sparginfos_.end() != sparginfos_.find(argvar))
      sparginfos_[argvar].MarkInvalid();
  }

  void FinalizeSignpows() {
    assert(sparginfos_.size());
    for (const auto& [var, info]: sparginfos_) {
      assert(info.totalabspow_>0);
      assert(info.abspow_.size());
      if (info.IsValid()
          && info.totaljustpow_>0     // justpow is odd positive
          && 1.0==std::fmod(info.totaljustpow_, 2.0)) {
        auto term1 = info.abspow_.begin();  // to be replaced
        assert(info.justpow_.size());
        for (const auto& [iTerm, pw]: info.abspow_)
          isRemoved_[iTerm] = true;
        for (const auto& [iTerm, pw]: info.justpow_)
          isRemoved_[iTerm] = true;
        isRemoved_[term1->first] = false;   // keep the 1st term
        if (std::fmod(info.totalabspow_, 2.0)) { // abspow not even
          std::get<1>(terms_flt_[term1->first])
              = typename FlatExpr::Variable {
                  GetFlt().GetFlatCvt().AssignResultVar2Args(
                      SignpowConstExpConstraint
                      { VarArray1{var},
                       DblParamArray1
                       {info.totalabspow_ + info.totaljustpow_}})
              };
        } else {                            // abspow even
          std::get<1>(terms_flt_[term1->first])
              = typename FlatExpr::Variable {
                  GetFlt().GetFlatCvt().AssignResultVar2Args(
                      PowConstExpConstraint
                      { VarArray1{var},
                       DblParamArray1
                       {info.totalabspow_ + info.totaljustpow_}})
              };
        }
      }
    }
  }

  /// @todo signpow factors - skip
  void CombineRemainingFactors() {
    result *= coef00;
    for ( ; iTerm<(int)terms_flt_.size(); ++iTerm) {
      if (!isRemoved_.at(iTerm))
        result = GetFlt().QuadratizeOrLinearize(
            result, std::get<1>(terms_flt_[iTerm]));
    }
  }

  /// Obtain flattener, const
  const Flattener& GetFlt() const { return flt_; }
  /// Obtain flattener
  Flattener& GetFlt() { return flt_; }


private:
  Flattener& flt_;
  /// tuple: comparator, term, bounds
  SmallVec<
      std::tuple<
          TermCmp, FlatExpr, std::pair<double, double> >,
      32 >                 // 32 elements preallocated
      terms_flt_;
  int n_terms_const_ = 0;
  int n_terms_binary_ = 0;

  FlatExpr result = typename FlatExpr::Constant {1.0};
  int iTerm=0;
  double coef00 = 1.0;

  /// Here the key is either an indep. var.,
  /// or the resvar of an LFC/QFC constraint
  std::unordered_map<int, SignpowArgInfo> sparginfos_;

  SmallVec<bool, 32> isSignpowAbs_;
  SmallVec<bool, 32> isRemoved_;
};

}  // namespace mp

#endif // PREPRO_PROD_H
