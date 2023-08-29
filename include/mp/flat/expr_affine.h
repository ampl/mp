#ifndef AFFINE_EXPR_H
#define AFFINE_EXPR_H

#include <array>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <cassert>

#include "mp/arrayref.h"
#include "mp/flat/expr_algebraic.h"

namespace mp {

/// Linear terms: c'x
class LinTerms {
public:
  /// Name
  static constexpr const char* GetTypeName() { return "LinTerms"; }

  /// Default constructor
  LinTerms() = default;

  /// Construct from 2 vectors
  LinTerms(std::vector<double> c, std::vector<int> v) noexcept
    : coefs_(std::move(c)), vars_(std::move(v))
  { assert(check()); }

  /// Construct from 2 std::array's
  template <size_t N>
  LinTerms(std::array<double, N>& c, std::array<int, N>& v)
    : coefs_(c.begin(), c.end()), vars_(v.begin(), v.end())
  { assert(check()); }

  /// Validate
  bool check() const {
    return coefs_.size()==vars_.size() &&
        (!size() || 0<=*std::min_element(vars_.begin(), vars_.end()));
  }

  /// empty()
  bool empty() const { return coefs_.empty(); }
  /// size()
  size_t size() const { return coefs_.size(); }
  /// coef[i]
  double coef(size_t i) const { assert(i<size()); return coefs_[i]; }
  /// var[i]
  int var(size_t i) const { assert(i<size()); return vars_[i]; }
  /// const vec& coefs()
  const std::vector<double>& coefs() const { return coefs_; }
  /// const vec& vars()
  const std::vector<int>& vars() const { return vars_; }
  /// Ptr to coefs
  const double* pcoefs() const { return coefs_.data(); }
  /// Ptr to vars
  const int* pvars() const { return vars_.data(); }

  /// true when 1 variable with coef 1.0
  bool is_variable() const {
    return 1==size() && 1.0==this->coef(0);
  }

  /// return the single variable assuming true==is_variable()
  int get_representing_variable() const {
    assert(is_variable());
    return this->var(0);
  }

  /// Always linear
  static constexpr bool is_linear() { return true; }

  /// Produce itself
  const LinTerms& GetLinTerms() const { return *this; }

  /// Compute value given a dense vector of variable values
  template <class VarInfo>
  long double ComputeValue(const VarInfo& x) const {
    long double s=0.0;
    for (size_t i=coefs().size(); i--; )
      s += (long double)(coefs()[i]) * x[vars()[i]];
    return s;
  }

  /// Set coef
  void set_coef(size_t i, double c)
  { assert(i<size()); coefs_[i]=c; }

  /// Reserve size
  void reserve(size_t s) {
    coefs_.reserve(s);
    vars_.reserve(s);
  }

  /// Add linear term
  void add_term(double c, int v)
  { coefs_.push_back(c); vars_.push_back(v); }

  /// Add another LinTerms
  void add(const LinTerms& le) {
    reserve(size() + le.size());
    for (size_t i=0; i<le.size(); ++i)
      add_term(le.coefs_[i], le.vars_[i]);
  }

  /// Add terms from a vector of pairs {c, v}
  template <class Vec>
  void add_terms(const Vec& v2) {
    reserve(size() + v2.size());
    for (const auto& term: v2)
      add_term(term.first, term.second);
  }

  /// Is normalized? Assume terms are sorted.
  bool is_normalized() const {
    assert(size());
    return coef(0) > 0.0;
  }

  /// Negate
  void negate() {
    for (auto& c: coefs_)
      c = -c;
  }

  /// Multiply by const
  void operator*=(double n) {
    for (auto& c: coefs_)
      c *= n;
  }

  /// preprocess / canonicalize
  void preprocess() { sort_terms(); }

  /// This a a NASTY one (when not used).
  /// Use it before adding constraints / objectives.
  /// Add same variables, eliminate 0's.
  /// Can be used by LinCon's etc
  /// Gurobi complains when 0's / repeated entries.
  void sort_terms(bool force_sort=false);

  /// Equality. Assumes being sorted
  bool equals(const LinTerms& lt) const {
    return coefs_==lt.coefs_ && vars_==lt.vars_;
  }

  /// operator== for hashing and testing
  bool operator==(const LinTerms& lt) const {
    return equals(lt);
  }


private:
  std::vector<double> coefs_;
  std::vector<int> vars_;
};


/// Typedef AffineExpr
using AffineExpr = AlgebraicExpression<LinTerms>;

} // namespace mp

#endif // AFFINE_EXPR_H
