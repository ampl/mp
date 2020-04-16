#ifndef AFFINE_EXPR_H
#define AFFINE_EXPR_H

#include <cassert>

#include "mp/problem.h"

namespace mp {

class AffineExpr : public LinearExpr {
  double constant_term_ = 0.0;
public:
  AffineExpr() {}
  AffineExpr(AffineExpr&& ae) = default;
  AffineExpr& operator = (AffineExpr&& ae) = default;
  /// Helper struct to construct AffineExpr from something special
  template <class Value>
  struct ConstructorHelper {
    Value v;
  };
  using Constant = ConstructorHelper<double>;
  using Variable = ConstructorHelper<int>;
  AffineExpr(Constant c) : constant_term_(c.v) {}
  AffineExpr(Variable i) { AddTerm(i.v, 1.0); }

  bool is_constant() const { return 0==num_terms(); }
  /// true when constant=0 and 1 variable with coef 1.0
  bool is_variable() const {
    return 0.0==std::fabs(constant_term()) && 1==num_terms() &&
        1.0==this->begin()->coef();
  }
  /// return the single variable assuming true==is_variable()
  int get_representing_variable() const {
    assert(is_variable());
    return this->begin()->var_index();
  }

  double constant_term() const { return constant_term_; }
  void constant_term(double v) { constant_term_ = v; }
  void add_to_constant(double a) { constant_term_ += a; }

  void Negate() {
    for (auto& term: *this)
      term.set_coef(-term.coef());
    constant_term(-constant_term());
  }

  void Add(const AffineExpr& ae) {
    this->Reserve(this->num_terms() + ae.num_terms());
    this->AddTerms(ae); // eliminate duplicates when?
    this->add_to_constant(ae.constant_term());
  }
};

} // namepsace mp

#endif // AFFINE_EXPR_H
