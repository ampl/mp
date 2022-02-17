#ifndef AFFINE_EXPR_H
#define AFFINE_EXPR_H

#include <array>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cassert>

namespace mp {

/// Linear terms: c'x
/// TODO use SSO
class LinTerms {
public:
  /// Default constructor
  LinTerms() = default;
  /// Construct from 2 vectors
  template <class CV=std::vector<double>, class VV=std::vector<int> >
  LinTerms(CV&& c, VV&& v) noexcept
    : coefs_(std::forward<CV>(c)), vars_(std::forward<VV>(v))
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
  void add_lin_exp(const LinTerms& le) {
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

  /// Add same variables, eliminate 0's.
  /// Can be used by LinCon's etc
  /// Gurobi complains when 0's / repeated entries.
  /// TODO use hash map when sorting not needed?
  void sort_terms(bool force_sort=false);


private:
  std::vector<double> coefs_;
  std::vector<int> vars_;
};

class AffExp : public LinTerms {
public:
  /// Default constructor
  AffExp() = default;
  /// From LinTerms&&
  AffExp(LinTerms ae, double ct=0.0) noexcept :
    LinTerms(std::move(ae)), constant_term_(ct) { }
  /// From const AE&
  AffExp(const AffExp& ae) = default;
  /// From AE&&
  AffExp(AffExp&& ae) = default;
  /// Assign AE&&
  AffExp& operator= (AffExp&& ae) = default;

  /// Helper to construct AffExp from something special
  template <class Value>
  struct ConstructorHelper {
    Value v;
  };
  using Constant = ConstructorHelper<double>;
  using Variable = ConstructorHelper<int>;
  AffExp(Constant c) : constant_term_(c.v) { }
  AffExp(Variable i) { add_term(1.0, i.v); }
  AffExp& operator+=(std::pair< std::pair<double, int>, double > a) {
    add_term(a.first.first, a.first.second);
    add_to_constant(a.second);
    return *this;
  }

  /// Whether AffExp represents a constant
  bool is_constant() const { return 0==size(); }
  /// true when constant=0 and 1 variable with coef 1.0
  bool is_variable() const {
    return 0.0==std::fabs(constant_term()) && 1==size() &&
        1.0==this->coefs()[0];
  }
  /// return the single variable assuming true==is_variable()
  int get_representing_variable() const {
    assert(is_variable());
    return this->vars()[0];
  }

  /// Get the lin exp, const
  const LinTerms& get_lin_exp() const { return (const LinTerms&)(*this); }

  /// The constant term
  double constant_term() const { return constant_term_; }
  /// Set constant term
  void constant_term(double v) { constant_term_ = v; }
  /// Add to constant
  void add_to_constant(double a) { constant_term_ += a; }

  /// Negate the ae
  void negate() {
    LinTerms::negate();
    constant_term(-constant_term());
  }
  /// Add another ae
  void add_aff_exp(const AffExp& ae) {
    add_lin_exp(ae);                // eliminate duplicates when?
    this->add_to_constant(ae.constant_term());
  }
  /// Subtract another ae
  void subtract(AffExp ae) {
    ae.negate();
    add_aff_exp(ae);
  }
  /// Multiply by const
  void operator*=(double n) {
    ((LinTerms&)(*this)) *= n;
    constant_term_ *= n;
  }
private:
  double constant_term_ = 0.0;
};

} // namepsace mp

#endif // AFFINE_EXPR_H
