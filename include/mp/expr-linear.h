#ifndef EXPR_LINEAR_H
#define EXPR_LINEAR_H

#include "mp/utils-vec.h"

namespace mp {

/// Linear expression (not affine: no constant term)
/// used in `mp::BasicProblem<>`
class LinearExpr {
private:
  class Term {
  private:
    int var_index_;
    double coef_;

    friend class LinearExpr;

    Term(int var_index, double coef) : var_index_(var_index), coef_(coef) {}

  public:
    int var_index() const { return var_index_; }
    double coef() const { return coef_; }
    void set_coef(double c) { coef_=c; }
    void operator*=(double n) { coef_*=n; }
  };
  /// Typedef term vector
  using TermVec = SmallVec<Term, 6>;
  TermVec terms_;

public:
  LinearExpr() { }

  template <class CoefVec=std::vector<double>, class VarVec=std::vector<int> >
  LinearExpr(CoefVec&& c, VarVec&& v) {
    ConstructFrom(std::forward<CoefVec>(c), std::forward<VarVec>(v));
  }
  template <int N>
  LinearExpr(const std::array<double, N>& c, const std::array<int, N>& v) {
    ConstructFrom(c, v);
  }

  int num_terms() const { return static_cast<int>(terms_.size()); }
  int capacity() const { return static_cast<int>(terms_.capacity()); }

  int var_index(int i) const { return terms_[i].var_index(); }
  double coef(int i) const { return terms_[i].coef(); }
  void set_coef(int i, double c) { terms_[i].set_coef(c); }

  typedef TermVec::const_iterator const_iterator;

  const_iterator begin() const { return terms_.begin(); }
  const_iterator end() const { return terms_.end(); }

  typedef TermVec::iterator iterator;

  iterator begin() { return terms_.begin(); }
  iterator end() { return terms_.end(); }

  void AddTerm(int var_index, double coef) {
    terms_.push_back(Term(var_index, coef));
  }

  void AddTerms(const LinearExpr& li) {
    terms_.insert(end(), li.begin(), li.end());
  }

  template <class CoefVec=std::vector<double>, class VarVec=std::vector<int> >
  void ConstructFrom(CoefVec&& c, VarVec&& v) {
    assert(c.size()==v.size());
    assert(0==num_terms());
    Reserve(c.size());
    for (size_t i=0; i<c.size(); ++i)
      AddTerm(v[i], c[i]);
  }

  void Reserve(std::size_t num_terms) {
    terms_.reserve(num_terms);
  }

  void SortTerms();
};

}  // namespace mp

#endif  // EXPR_LINEAR_H
