#ifndef EXPR2CONSTRAINT_H
#define EXPR2CONSTRAINT_H

#include "mp/problem.h"
#include "mp/convert/affine_expr.h"
#include "mp/convert/constraint_keeper.h"

namespace mp {

class EExpr : public AffineExpr {   // Result expression for expression conversions
public:
  EExpr() = default;
  EExpr(Constant c) : AffineExpr(c) {}
  EExpr(Variable v) : AffineExpr(v) {}
  EExpr(int i, double c) { AddTerm(i, c); }
};

/// Helper class providing a default framework for converting an expression
/// to a variable plus a constraint equating that variable to that expression
template <class Impl, class Converter, class Constraint>
class BasicExprToConstraintConverter {
  Converter& converter_;
public:
  using ArgArray = typename Constraint::ArgArray;
  ArgArray args_;                       // variables
  int result_var_;
protected:
  Converter& GetConverter() { return converter_; }
  ArgArray&& MoveOutArguments() {        // this returns rvalue and invalidates args_
    return std::move(args_);
  }
  void AddArgument(int v) { args_.push_back(v); }
  int GetResultVar() { return result_var_; }
  void SetResultVar(int v) { result_var_ = v; }
public:
  BasicExprToConstraintConverter(Converter& cvt) : converter_(cvt) { }
  template <class ExprArray>
  EExpr ConvertArray(ExprArray e) {
    MP_DISPATCH( ConvertArguments(e) );
    MP_DISPATCH( AddConstraint() );
    return EExpr::Variable{ MP_DISPATCH( GetResultVar() ) };
  }
};

template <class Impl, class Converter, class Constraint>
class Expr2Constr : public BasicExprToConstraintConverter<Impl, Converter, Constraint> {
  using Base = BasicExprToConstraintConverter<Impl, Converter, Constraint>;
public:
  Expr2Constr(Converter& cvt) : Base(cvt) { }
  template <class ExprArray>
  void ConvertArguments(ExprArray ea) {
    for (const auto& e: ea)
      MP_DISPATCH( AddArgument(MP_DISPATCH( GetConverter() ).Convert2Var(e)) );
  }
  void AddConstraint() {
    MP_DISPATCH( SetResultVar( MP_DISPATCH( GetConverter() ).AddVar() ) );
    /// TODO propagate bounds from arguments
    MP_DISPATCH( GetConverter() ).AddConstraint(
          makeConstraint<Converter, Constraint>
            (MP_DISPATCH( MoveOutArguments() ), MP_DISPATCH( GetResultVar() )));
  }
};

/// This is a helper to produce a 'final' E2C converter avoiding using Impl
template <template <typename, typename, typename> class E2C,
          class Converter, class Constraint>
class E2CImpl : public E2C<E2CImpl<E2C, Converter, Constraint>, Converter, Constraint > {
  using Base = E2C<E2CImpl<E2C, Converter, Constraint>, Converter, Constraint >;
public:
  E2CImpl(Converter& cvt) : Base(cvt) { }
};

template <template <typename, typename, typename> class E2C,
          class Converter, class Constraint, class Converter2>
E2CImpl<E2C, Converter, Constraint >
makeE2CConverter(Converter2& cvt) {
  return E2CImpl<E2C, Converter, Constraint>( static_cast<Converter&>(cvt) );
}

} // namespace mp

#endif // EXPR2CONSTRAINT_H
