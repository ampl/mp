#ifndef QUAD_EXPR_H
#define QUAD_EXPR_H

#include <tuple>

#include "mp/flat/expr_affine.h"

namespace mp {

/// Quadratic terms x'Qx
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
  int size() const { return (int)coefs_.size(); }

  /// Capacity
  int capacity() const { return (int)coefs_.capacity(); }

  const double* pcoefs() const { return coefs_.data(); }
  const int* pvars1() const { return vars1_.data(); }
  const int* pvars2() const { return vars2_.data(); }

  const std::vector<double>& coefs() const { return coefs_; }
  const std::vector<int>& vars1() const { return vars1_; }
  const std::vector<int>& vars2() const { return vars2_; }

  double coef(int i) const { return coefs_[i]; }
  void set_coef(int i, double c) { coefs_[i] = c; }
  int var1(size_t i) const { return vars1_[i]; }
  int var2(size_t i) const { return vars2_[i]; }

  /// Compute value given a dense vector of variable values
  template <class VarInfo>
  long double ComputeValue(const VarInfo& x) const {
    long double s=0.0;
    for (size_t i=coefs().size(); i--; )
      s += (long double)(coefs()[i]) * x[vars1()[i]] * x[vars2()[i]];
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

  /// Is normalized? Assume sorted.
  bool is_normalized() const {
    assert(size());
    return coef(0) > 0.0;
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

  /// Clear
  void clear() {
    coefs_.clear();
    vars1_.clear();
    vars2_.clear();
  }

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


////////////////////////////////////////////////////////////////////////
/// Quadratic and linear terms.
/// Body of a quadratic constraint
class QuadAndLinTerms :
    protected LinTerms, protected QuadTerms {
public:
  /// Name
  static constexpr const char* GetTypeName() { return "QuadAndLinTerms"; }

  /// Default constructor
  QuadAndLinTerms() = default;

  /// Construct from linear + QP terms
  QuadAndLinTerms(LinTerms lt, QuadTerms qt) :
    LinTerms(std::move(lt)), QuadTerms(std::move(qt)) {
    sort_terms();
  }

  /// Get LinTerms, const
  const LinTerms& GetLinTerms() const { return *this; }
  /// Get LinTerms
  LinTerms& GetLinTerms() { return *this; }

  /// Get QuadTerms, const
  const QuadTerms& GetQPTerms() const { return *this; }
  /// Get QuadTerms
  QuadTerms& GetQPTerms() { return *this; }

  /// empty?
  bool empty() const { return GetLinTerms().empty() && GetQPTerms().empty(); }

  /// is linear?
  bool is_linear() const { return GetQPTerms().empty(); }

  /// Is quadratic?
  bool is_quadratic() const { return !is_linear(); }

  /// true when 1 linear variable with coef 1.0
  bool is_variable() const {
    return is_linear() && GetLinTerms().is_variable();
  }

  /// return the single variable assuming true==is_variable()
  int get_representing_variable() const {
    assert(is_variable());
    return GetLinTerms().get_representing_variable();
  }

  /// add_term(c, v)
  using LinTerms::add_term;

  /// add_term(c, v1, v2)
  using QuadTerms::add_term;

  /// Is normalized? Assume sorted.
  bool is_normalized() const {
    assert(QuadTerms::size());
    return
        LinTerms::size()
        ? LinTerms::is_normalized()
        : QuadTerms::is_normalized();
  }

  /// Negate
  void negate() {
    LinTerms::negate();
    QuadTerms::negate();
  }

  /// add body
  void add(const QuadAndLinTerms& qlt) {
    LinTerms::add(qlt.GetLinTerms());
    QuadTerms::add(qlt.GetQPTerms());
  }

  /// Value at given variable vector
  template <class VarInfo>
  long double ComputeValue(const VarInfo& x) const {
    return LinTerms::ComputeValue(x) + QuadTerms::ComputeValue(x);
  }

  /// Sort terms
  void sort_terms() {
    LinTerms::sort_terms();
    QuadTerms::sort_terms();
  }

  /// Test equality
  bool equals(const QuadAndLinTerms& qlc) const {
    return LinTerms::equals(qlc.GetLinTerms()) &&
        QuadTerms::equals(qlc.GetQPTerms());
  }

  /// Test equality
  bool operator==(const QuadAndLinTerms& qlc) const { return equals(qlc); }
};


/// Typedef QuadraticExpr
using QuadraticExpr = AlgebraicExpression<QuadAndLinTerms>;


/// Extract (move out) affine expr
inline AffineExpr MoveOutAffineExpr(QuadraticExpr&& qe) {
  return { std::move(qe.GetLinTerms()), qe.constant_term() };
}

} // namespace mp

#endif // QUAD_EXPR_H
