#ifndef AFFINE_EXPR_H
#define AFFINE_EXPR_H

#include <array>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <functional>
#include <cassert>

#include "mp/arrayref.h"
#include "mp/utils-vec.h"
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
  LinTerms(const std::vector<double>& c, const std::vector<int>& v)
      : coefs_(c.begin(), c.end()), vars_(v.begin(), v.end())
  { assert(check()); }

  /// Construct from 2 std::array's
  template <size_t N>
  LinTerms(std::array<double, N>& c, std::array<int, N>& v)
    : coefs_(c.begin(), c.end()), vars_(v.begin(), v.end())
  { assert(check()); }

  /// Construct from 2 SmallVec's
  template <unsigned int N>
  LinTerms(const SmallVec<double, N>& c, const SmallVec<int, N>& v)
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
  ArrayRef<double> coefs() const { return {coefs_.data(), coefs_.size()}; }
  /// const vec& vars()
  ArrayRef<int> vars() const { return {vars_.data(), vars_.size()}; }
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

  /// Clear
  void clear() {
    coefs_.clear();
    vars_.clear();
  }

  /// Reserve size
  void reserve(size_t s) {
    coefs_.reserve(s);
    vars_.reserve(s);
  }

  /// Resize
  void resize(size_t s) {
    coefs_.resize(s);
    vars_.resize(s);
  }

  /// shrink_to_fit.
  /// Takes time, so use only when necessary.
  void shrink_to_fit() {
    coefs_.shrink_to_fit();
    vars_.shrink_to_fit();
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

  /// Fold the terms into a vector of (var, coef) pairs
  template <class Vec>
  void fold_into(Vec& vec);

  /// Unfold the terms from a vector of (var, coef) pairs
  template <class Vec>
  void unfold_from(const Vec& vec);

  /// preprocess / canonicalize
  void preprocess() { sort_terms(); }

  /// Is the expression sorted,
  /// all elements non-0 and unique?
  bool is_sorted() const;

  /// Is the expression sorted,
  /// all elements unique but maybe 0 coefs?
  bool is_sorted__maybe_0s() const;

  /// This a NASTY one (when not used).
  /// Use it before adding
  /// constraints / objectives / expressions.
  /// Unify same variables, eliminate 0's.
  /// Can be used by LinCon's etc
  /// Gurobi complains when 0's / repeated entries.
  /// @param force_sort: sort also when
  /// all elements unique and non-0.
  void sort_terms(bool force_sort=true);

  /// Sort and unify but leave 0's
  /// which can represent sparsity pattern
  void sort_terms__leave_0s();

  /// (var, coef)
  std::pair<int, double> IndexValue(size_t i) const
  { return {var(i), coef(i)}; }

  /// Add (var, coef)
  void add_index_value(std::pair<int, double> iv)
  { add_term(iv.second, iv.first); }

  /// Set [i] = (var, coef)
  void set_index_value(size_t i, std::pair<int, double> iv) {
    assert(i<size());
    vars_[i] = iv.first; coefs_[i] = iv.second;
  }

  /// Equality. Assumes being sorted
  bool equals(const LinTerms& lt) const {
    return coefs_==lt.coefs_ && vars_==lt.vars_;
  }

  /// operator== for hashing and testing
  bool operator==(const LinTerms& lt) const {
    return equals(lt);
  }

  /// operator!=
  bool operator!=(const LinTerms& lt) const
  { return !(*this==lt); }


private:
  SmallVec<double, 6> coefs_;
  SmallVec<int, 6> vars_;
};

/// Merge 2 sorted LinTerms
LinTerms Merge(const LinTerms& , const LinTerms& );

/// Specialize
template <>
void WriteJSON(JSONW jw, const LinTerms& qt);

/// Specialize
void VisitArguments(const LinTerms& lt, std::function<void (int) > argv);

/// Typedef AffineExpr
using AffineExpr = AlgebraicExpression<LinTerms>;

} // namespace mp

#endif // AFFINE_EXPR_H
