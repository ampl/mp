#ifndef QUAD_EXPR_H
#define QUAD_EXPR_H

#include <tuple>

#include "mp/flat/expr_affine.h"

namespace mp {

/// Quadratic terms
class QuadTerms {
public:
  /// Name
  static constexpr const char* GetTypeName() { return "QuadTerms"; }

  /// Default constructor
  QuadTerms() = default;

  /// Construct from 3 vectors
  QuadTerms(std::vector<double> c,
           std::vector<int> v1, std::vector<int> v2) noexcept
    : coefs_(std::move(c)), vars1_(std::move(v1)), vars2_(std::move(v2))
  { assert(check()); }

  /// Validate
  bool check() const {
    return
        coefs_.size()==vars1_.size() &&
        coefs_.size()==vars2_.size() &&
        (!size() || 0<=*std::min_element(vars1_.begin(), vars1_.end())) &&
        (!size() || 0<=*std::min_element(vars2_.begin(), vars2_.end()));
  }

  /// Empty?
  bool empty() const { return coefs_.empty(); }

  /// Size
  size_t size() const { return coefs_.size(); }

  /// Capacity
  size_t capacity() const { return coefs_.capacity(); }

  const double* pcoefs() const { return coefs_.data(); }
  const int* pvars1() const { return vars1_.data(); }
  const int* pvars2() const { return vars2_.data(); }

  const std::vector<double>& coefs() const { return coefs_; }
  const std::vector<int>& vars1() const { return vars1_; }
  const std::vector<int>& vars2() const { return vars2_; }

  double coef(int i) const { return coefs_[i]; }
  void set_coef(int i, double c) { coefs_[i] = c; }
  int var1(int i) const { return vars1_[i]; }
  int var2(int i) const { return vars2_[i]; }

  /// Compute value given a dense vector of variable values
  double ComputeValue(ArrayRef<double> x) const {
    double s=0.0;
    for (size_t i=coefs().size(); i--; )
      s += coefs()[i] * x[vars1()[i]] * x[vars2()[i]];
    return s;
  }

  void add_term(double coef, int var1, int var2) {
    coefs_.push_back(coef);
    vars1_.push_back(var1);
    vars2_.push_back(var2);
  }

  void reserve(std::size_t num_terms) {
    coefs_.reserve(num_terms);
    vars1_.reserve(num_terms);
    vars2_.reserve(num_terms);
  }

  /// Arithmetic
  void negate() {
    for (auto& cf: coefs_)
      cf = -cf;
  }

  void add(const QuadTerms& li) {
    this->reserve(size() + li.size());
    /// eliminate duplicates when?
    coefs_.insert(coefs_.end(), li.coefs_.begin(), li.coefs_.end());
    vars1_.insert(vars1_.end(), li.vars1_.begin(), li.vars1_.end());
    vars2_.insert(vars2_.end(), li.vars2_.begin(), li.vars2_.end());
  }

  void subtract(QuadTerms&& ae) {
    ae.negate();
    add(ae);
  }

  void operator*=(double n) {
    for (auto& c: coefs_)
      c *= n;
  }

  /// Sort and eliminate duplicates
  void sort_terms();

  /// Test equality
  bool equals(const QuadTerms& qt) const {
    return *this == qt;
  }

  /// Testing API
  bool operator==(const QuadTerms& qt) const {
    return coefs_==qt.coefs_ && vars1_==qt.vars1_ && vars2_==qt.vars2_;
  }


private:
  std::vector<double> coefs_;
  std::vector<int> vars1_;
  std::vector<int> vars2_;
};


/// Quadratic expression
class QuadExp {
public:
  /// Name
  static constexpr const char* GetTypeName() { return "QuadExp"; }

  /// Default constructor
  QuadExp() = default;

  /// Construct from AE [+QT]
  QuadExp(AffExp ae, QuadTerms qt={}) noexcept :
    ae_(std::move(ae)), qt_(std::move(qt)) { }

  using Constant = AffExp::Constant;
  using Variable = AffExp::Variable;

  /// Getters
  bool is_constant() const { return is_affine() && GetAE().is_constant(); }
  double constant_term() const { return GetAE().constant_term(); }
  bool is_variable() const { return is_affine() && GetAE().is_variable(); }
  int get_representing_variable() const {
    assert(is_affine());
    return GetAE().get_representing_variable();
  }
  /// Check if affine
  bool is_affine() const { return GetQT().empty(); }
  /// Get affine expr, const
  const AffExp& GetAE() const { return ae_; }
  /// Get affine expr
  AffExp& GetAE() { return ae_; }
  /// Is quadratic?
  bool is_quadratic() const { return !is_affine(); }
  /// Get quadratic terms, const
  const QuadTerms& GetQT() const { return qt_; }
  /// Get quadratic terms
  QuadTerms& GetQT() { return qt_; }


  /// Modifiers
  void constant_term(double v) { GetAE().constant_term(v); }
  void add_to_constant(double a) { GetAE().add_to_constant(a); }
  void add_linear_term(double c, int v) {
    GetAE().add_term(c, v);
  }
  void add_qp_term(double coef, int v1, int v2) {
    GetQT().add_term(coef, v1, v2);
  }

  void add(const QuadExp& qe) {
    GetAE().add_aff_exp(qe.GetAE());
    GetQT().add(qe.GetQT());
  }

  void subtract(QuadExp&& qe) {
    GetAE().subtract(std::move(qe.GetAE()));
    GetQT().subtract(std::move(qe.GetQT()));
  }

  void negate() {
    GetAE().negate();
    GetQT().negate();
  }

  void sort_terms() {
    GetAE().sort_terms();
    GetQT().sort_terms();
  }


private:
  AffExp ae_;
  QuadTerms qt_;
};

}

#endif // QUAD_EXPR_H
