#ifndef QUAD_EXPR_H
#define QUAD_EXPR_H

#include <map>

#include "mp/convert/affine_expr.h"

namespace mp {

class QuadTerms {
private:
  std::vector<double> coefs_;
  std::vector<int> vars1_;
  std::vector<int> vars2_;

public:
  QuadTerms() { }
  QuadTerms(std::initializer_list<std::tuple<double, int, int>> quad_terms) {
    Reserve(quad_terms.size());
    for (const auto& term: quad_terms)
      AddTerm(std::get<0>(term), std::get<1>(term), std::get<2>(term));
    sort_terms();
  }

  bool empty() const { return coefs_.empty(); }
  int num_terms() const { return static_cast<int>(coefs_.size()); }
  int capacity() const { return static_cast<int>(coefs_.capacity()); }

  const double* coefs() const { return coefs_.data(); }
  const int* vars1() const { return vars1_.data(); }
  const int* vars2() const { return vars2_.data(); }

  double coef(int i) const { return coefs_[i]; }
  void set_coef(int i, double c) { coefs_[i] = c; }
  int var1(int i) const { return vars1_[i]; }
  int var2(int i) const { return vars2_[i]; }

  void AddTerm(double coef, int var1, int var2) {
    coefs_.push_back(coef);
    vars1_.push_back(var1);
    vars2_.push_back(var2);
  }

  void AddTerms(const QuadTerms& li) {
    coefs_.insert(coefs_.end(), li.coefs_.begin(), li.coefs_.end());
    vars1_.insert(vars1_.end(), li.vars1_.begin(), li.vars1_.end());
    vars2_.insert(vars2_.end(), li.vars2_.begin(), li.vars2_.end());
  }

  void Reserve(std::size_t num_terms) {
    coefs_.reserve(num_terms);
    vars1_.reserve(num_terms);
    vars2_.reserve(num_terms);
  }

  /// Arithmetic
  void Negate() {
    for (auto& cf: coefs_)
      cf = -cf;
  }

  void Add(const QuadTerms& ae) {
    this->Reserve(this->num_terms() + ae.num_terms());
    this->AddTerms(ae); // eliminate duplicates when?
  }

  void Subtract(QuadTerms&& ae) {
    ae.Negate();
    Add(ae);
  }


  void sort_terms() {
    auto sort_pair = [](int a, int b) {
      return a<b ? std::pair<int, int>(a, b) : std::pair<int, int>(b, a);
    };
    std::map<std::pair<int, int>, double> var_coef_map;
    for (int i=0; i<num_terms(); ++i)
      if (0.0!=std::fabs(coefs_[i]))
        var_coef_map[sort_pair(vars1_[i], vars2_[i])] += coefs_[i];
    if (true) {                                // would check size if using hash map
      coefs_.clear();
      vars1_.clear();
      vars2_.clear();
      for (const auto& vc: var_coef_map) {
        if (0.0!=std::fabs(vc.second))         // Need tolerance?
          AddTerm(vc.second, vc.first.first, vc.first.second);
      }
    }
  }

  /// Testing API
  bool operator==(const QuadTerms& qt) const {
    return coefs_==qt.coefs_ && vars1_==qt.vars1_ && vars2_==qt.vars2_;
  }
};

class QuadExpr {
public:
  QuadExpr() { }
  QuadExpr(AffineExpr&& ae) : ae_(std::move(ae)) { }

  using Constant = AffineExpr::Constant;
  using Variable = AffineExpr::Variable;

  /// Getters
  bool is_constant() const { return is_affine() && GetAE().is_constant(); }
  double constant_term() const { return GetAE().constant_term(); }
  bool is_variable() const { return is_affine() && GetAE().is_variable(); }
  int get_representing_variable() const {
    assert(is_affine());
    return GetAE().get_representing_variable();
  }

  bool is_affine() const { return GetQT().empty(); }
  const AffineExpr& GetAE() const { return ae_; }
  AffineExpr& GetAE() { return ae_; }
  bool is_quadratic() const { return !is_affine(); }
  const QuadTerms& GetQT() const { return qt_; }
  QuadTerms& GetQT() { return qt_; }


  /// Modifiers
  void constant_term(double v) { GetAE().constant_term(v); }
  void add_to_constant(double a) { GetAE().add_to_constant(a); }
  void AddLinearTerm(int var_index, double coef) {
    GetAE().AddTerm(var_index, coef);
  }
  void AddQuadraticTerm(int v1, int v2, double coef) {
    GetQT().AddTerm(coef, v1, v2);
  }

  void Add(const QuadExpr& qe) {
    GetAE().Add(qe.GetAE());
    GetQT().Add(qe.GetQT());
  }

  void Subtract(QuadExpr&& qe) {
    GetAE().Subtract(std::move(qe.GetAE()));
    GetQT().Subtract(std::move(qe.GetQT()));
  }

  void Negate() {
    GetAE().Negate();
    GetQT().Negate();
  }
private:
  AffineExpr ae_;
  QuadTerms qt_;
};

}

#endif // QUAD_EXPR_H
