#ifndef EXPR2CONSTRAINT_H
#define EXPR2CONSTRAINT_H

#include "mp/problem.h"


namespace mp {

class EExpr : public LinearExpr {   // Result expression for expression conversions
public:
  EExpr(int i) { AddTerm(i, 1.0); }
};

/// Helper class providing a default framework for converting an expression
/// to a variable plus a constraint equating that variable to that expression
template <class Impl, class Converter, class Constraint>
class BasicExprToConstraintConverter {
  Converter& converter_;
  std::vector<EExpr> args_;
  int result_var_;
protected:
  Converter& GetConverter() { return converter_; }
  std::vector<EExpr>& GetArguments() { return args_; }
  void AddArgument(EExpr&& ee) { args_.push_back(ee); }
  int GetResultVar() { return result_var_; }
  void SetResultVar(int v) { result_var_ = v; }
public:
  BasicExprToConstraintConverter(Converter& cvt) : converter_(cvt) { }
  template <class ExprArray>
  EExpr ConvertArray(ExprArray e) {
    MP_DISPATCH( ConvertArguments(e) );
    /// TODO Preprocessing, Common Subexpression Elimination, Explanation
    MP_DISPATCH( AddConstraint() );
    return MP_DISPATCH( GetResultVar() );
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
      MP_DISPATCH( AddArgument(MP_DISPATCH( GetConverter() ).Visit(e)) );
  }
  void AddConstraint() {
    MP_DISPATCH( SetResultVar( MP_DISPATCH( GetConverter() ).AddVar().index() ) );
    /// TODO propagate bounds from arguments
    MP_DISPATCH( GetConverter() ).AddConstraint();
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
          class Converter, template <class, class> class Constraint, class Converter2>
E2CImpl<E2C, Converter, Constraint<Converter, typename Converter::BackendType> >
makeE2CConverter(Converter2& cvt) {
  return E2CImpl<E2C, Converter,
      Constraint<Converter, typename Converter::BackendType> >( static_cast<Converter&>(cvt) );
}

} // namespace mp

#endif // EXPR2CONSTRAINT_H
