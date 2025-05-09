#ifndef QUAD_EXPR_H
#define QUAD_EXPR_H

#include <tuple>

#include "mp/flat/expr_affine.h"

namespace mp {

/// Quadratic terms x'Qx.
class QuadTerms {
public:
  /// Name
  static constexpr const char* GetTypeName() { return "QuadTerms"; }

  /// Folded entry type
  using TupleType = std::pair<std::pair<int, int>, double>;

  /// Default constructor
  QuadTerms() = default;

  /// Construct from 3 vectors
  QuadTerms(const std::vector<double>& c,
            const std::vector<int>& v1, const std::vector<int>& v2);

  /// Empty?
  bool empty() const { return folded_.empty(); }

  /// Size
  size_t size() const { return folded_.size(); }

  /// Capacity
  size_t capacity() const { return folded_.capacity(); }

  /// Folded vector reference
  ArrayRef<TupleType> get_folded() const
  { return {folded_.data(), folded_.size()}; }

  /// coef vector pointer
  const double* pcoefs() const
  { unfold_if_need(); return coefs_aux_.data(); }
  const int* pvars1() const
  { unfold_if_need(); return vars1_aux_.data(); }
  const int* pvars2() const
  { unfold_if_need(); return vars2_aux_.data(); }

  /// Coef vector as ArrayRef<>
  ArrayRef<double> coefs() const
  { unfold_if_need(); return {coefs_aux_.data(), coefs_aux_.size()}; }
  ArrayRef<int> vars1() const
  { unfold_if_need(); return {vars1_aux_.data(), vars1_aux_.size()}; }
  ArrayRef<int> vars2() const
  { unfold_if_need(); return {vars2_aux_.data(), vars2_aux_.size()}; }

  /// coef(i)
  double coef(int i) const { return folded_[i].second; }
  void set_coef(int i, double c) { folded_[i].second = c; }
  int var1(size_t i) const { return folded_[i].first.first; }
  int var2(size_t i) const { return folded_[i].first.second; }

  /// Compute value given a dense vector of variable values
  template <class VarInfo>
  long double ComputeValue(const VarInfo& x) const {
    long double s=0.0;
    for (size_t i=size(); i--; )
      s += ((long double)(coef(i))) * x[var1(i)] * x[var2(i)];
    return s;
  }

  void add_term(double coef, int var1, int var2)
	{ add_index_value({{var1, var2}, coef}); clear_unfolded(); }

  void reserve(std::size_t num_terms)
  { folded_.reserve(num_terms); }

  /// shrink_to_fit.
  /// Takes time, so use only when necessary.
  void shrink_to_fit() {
    folded_.shrink_to_fit();
    coefs_aux_.shrink_to_fit();
    vars1_aux_.shrink_to_fit();
    vars2_aux_.shrink_to_fit();
  }

  /// Is normalized? Assume sorted.
  bool is_normalized() const {
    assert(size());
    return coef(0) > 0.0;
  }

  /// Arithmetic
  void negate() {
    for (auto& iv: folded_)
      iv.second = -iv.second;
		clear_unfolded();
  }

  void add(const QuadTerms& li) {
    this->reserve(size() + li.size());
    /// eliminate duplicates when?
    auto fld = li.get_folded();
    folded_.insert(folded_.end(), fld.begin(), fld.end());
		clear_unfolded();
  }

  void subtract(QuadTerms&& ae) {
    ae.negate();
    add(ae);
  }

  void operator*=(double n) {
    for (auto& c: folded_)
      c.second *= n;
		clear_unfolded();
  }

  /// Is the expression sorted,
  /// all elements non-0 and unique?
  bool is_sorted() const;

  /// Sort and eliminate duplicates
  void sort_terms();

  /// ({var1, var2}, coef)
  TupleType IndexValue(size_t i) const
  { assert(i<=size()); return folded_[i]; }

  /// Add ({var1, var2}, coef)
  void add_index_value(TupleType iv)
	{ folded_.push_back(iv); clear_unfolded(); }

  /// Clear
  void clear() {
    folded_.clear();
    clear_unfolded();
  }

  /// Test equality
  bool equals(const QuadTerms& qt) const {
    return *this == qt;
  }

  /// Testing API
  bool operator==(const QuadTerms& qt) const {
    return folded_ == qt.folded_;
  }

protected:
  /// Unfold if alternative empty
  void unfold_if_need() const;

  /// Unfold into the alternative
  void unfold() const;

  /// Clear alternative
  void clear_unfolded() {
    coefs_aux_.clear();
    vars1_aux_.clear();
    vars2_aux_.clear();
  }

  /// Sort index pairs
  template <class Vec>
  static void sort_index_pairs(Vec& );

  /// Fold the alternative terms into a vector of ((var1, var2), coef) tuples
  template <class Vec>
  void fold_into(Vec& vec);

  /// Unfold the terms from a vector of ((var1, var2), coef) tuples
  template <class Vec>
  void unfold_from(const Vec& vec) const;

  /// Fold the terms into a vector of ((var1, var2), coef) tuples
  template <class Vec, class ArrayC, class ArrayV>
  static void fold_into(Vec& vec,
                 const ArrayC& coefs,
                 const ArrayV& vars1, const ArrayV& vars2);

  /// Unfold the terms from a vector of ((var1, var2), coef) tuples
  template <class Vec, class ArrayC, class ArrayV>
  static void unfold_from(const Vec& vec,
                          ArrayC& coefs,
                          ArrayV& vars1, ArrayV& vars2);

private:
  /// Currently always there.
  /// @todo switch between alternatives but not too often.
  SmallVec<TupleType, 8> folded_;
  /// Alternative
  mutable SmallVec<double, 6> coefs_aux_;
  mutable SmallVec<int, 6> vars1_aux_;
  mutable SmallVec<int, 6> vars2_aux_;
};

/// Merge 2 sorted QuadTerms
QuadTerms Merge(const QuadTerms& , const QuadTerms& );

/// Multiply out two LinTerms.
/// Sort both factors - only in Debug,
/// then produce sorted result.
QuadTerms MultiplyOut(const LinTerms& e1, const LinTerms& e2);

/// Specialize
template <>
void WriteJSON(JSONW jw, const QuadTerms& qt);

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

  /// Multiply by const
  void operator*=(double n) {
    GetLinTerms() *= n;
    GetQPTerms() *= n;
  }

  /// Clear
  void clear() {
    GetLinTerms().clear();
    GetQPTerms().clear();
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

/// Specialize
template <>
void WriteJSON(JSONW jw, const QuadAndLinTerms& qt);

/// Specialize
void VisitArguments(const QuadTerms& lt, std::function<void (int) > argv);

/// Specialize
void VisitArguments(const QuadAndLinTerms& lt, std::function<void (int) > argv);


/// Typedef QuadraticExpr
using QuadraticExpr = AlgebraicExpression<QuadAndLinTerms>;


/// Extract (move out) affine expr
inline AffineExpr MoveOutAffineExpr(QuadraticExpr&& qe) {
  return { std::move(qe.GetLinTerms()), qe.constant_term() };
}

} // namespace mp

#endif // QUAD_EXPR_H
