#ifndef EXPR_ALGEBRAIC_H
#define EXPR_ALGEBRAIC_H

/**
 * A template algebraic expression:
 * Body (variable_terms) + constant_term.
 */

#include <string>
#include <cmath>

#include "mp/arrayref.h"

namespace mp {

/// A template algebraic expression:
/// Body (variable_terms) + constant_term.
template <class Body>
class AlgebraicExpression : public Body {
public:
  /// Name
  static std::string GetTypeName(){
    return std::string("AlgebraicExpression::") +
        Body::GetTypeName();
  }

  /// Typedef BodyType
  using BodyType = Body;

  /// Default constructor
  AlgebraicExpression() = default;

  /// From Body and const_term
  AlgebraicExpression(Body bt, double ct) noexcept :
    Body(std::move(bt)), constant_term_(ct) { }

  /// Helper to construct AlgebraicExpression
  /// to represent something special.
  /// @param Value: stored value type
  template <class Value>
  struct ConstructorHelper {
    Value v;
  };
  /// Helper "constant type"
  using Constant = ConstructorHelper<double>;
  /// Helper "variable type"
  using Variable = ConstructorHelper<int>;
  /// Constructor from helper "constant"
  AlgebraicExpression(Constant c) : constant_term_(c.v) { }
  /// Constructor from helper "variable"
  AlgebraicExpression(Variable i) { Body::add_term(1.0, i.v); }

  /// Whether AlgebraicExpression represents a constant
  bool is_constant() const { return Body::empty(); }
  /// true when constant=0 and 1 variable with coef 1.0
  bool is_variable() const {
    return 0.0==std::fabs(constant_term()) && GetBody().is_variable();
  }
  /// if affine
  bool is_affine() const { return Body::is_linear(); }

  /// Get the body, const
  const Body& GetBody() const { return (const Body&)(*this); }
  /// Get the body
  Body& GetBody() { return (Body&)(*this); }

  /// Get the body (variable terms)
  /// of a corresponding algebraic constraint
  const Body& GetAlgConBody() const { return GetBody(); }

  /// The constant term
  double constant_term() const { return constant_term_; }
  /// Set constant term
  void constant_term(double v) { constant_term_ = v; }
  /// Add to constant
  void add_to_constant(double a) { constant_term_ += a; }

  /// Compute value given a dense vector of variable values
  template <class VarInfo>
  double ComputeValue(const VarInfo& x) const
  { return Body::ComputeValue(x) + constant_term(); }

  /// Negate the ae
  void negate() {
    Body::negate();
    constant_term(-constant_term());
  }
  /// Add another ae
  void add(const AlgebraicExpression& ae) {
    Body::add(ae.GetBody());  // eliminate duplicates when? User code?
    this->add_to_constant(ae.constant_term());
  }
  /// Subtract another ae
  void subtract(AlgebraicExpression ae) {
    ae.negate();
    add(ae);
  }
  /// Multiply by const
  void operator*=(double n) {
    GetBody() *= n;
    constant_term_ *= n;
  }

  /// operator==
  bool operator==(const AlgebraicExpression& ae) const {
    return GetBody()==ae.GetBody() &&
        constant_term()==ae.constant_term();
  }

private:
  double constant_term_ = 0.0;
};

} // namespace mp

#endif // EXPR_ALGEBRAIC_H
