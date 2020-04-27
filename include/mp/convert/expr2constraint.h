#ifndef EXPR2CONSTRAINT_H
#define EXPR2CONSTRAINT_H

#include "mp/convert/affine_expr.h"

namespace mp {

/// Helper class providing a default framework for converting an expression
/// to a variable plus a constraint equating that variable to that expression
template <class Impl, class Converter, class Constraint>
class BasicFunctionalConstraintConverter {
  Converter& converter_;
public:
  using ArgArray = typename Converter::VarArray;
  ArgArray args_;
  bool result_is_known_ = false;
  int result_var_ = -1;
  double lb_ = Converter::MinusInfinity();
  double ub_ = Converter::PlusInfinity();
  var::Type type_ = var::CONTINUOUS;
protected:
  Converter& GetConverter() { return converter_; }
  double lb() const { return lb_; }
  double ub() const { return ub_; }
  var::Type type() const { return type_; }
public:
  ArgArray& Arguments() const { return args_; }
  bool ResultIsKnown() const { return result_is_known_; }
  int ResultVar() const { assert(result_var_>=-1); return result_var_; }
protected:
  void SetResultVar(int v) { result_var_ = v; }
public:
  BasicFunctionalConstraintConverter(Converter& cvt, ArgArray&& aa) :
    converter_(cvt), args_(aa) { }
  /// Convert array of arguments into a result variable
  /// plus constraint(s)
  int Convert() {
    MP_DISPATCH( AnalyzeArguments() );
    if (ResultIsKnown())
      return ResultVar();
    MP_DISPATCH( AddResultVariable() );
    MP_DISPATCH( AddConstraint() );
    return ResultVar();
  }
  void AnalyzeArguments() {

  }
  void AddResultVariable() {
    SetResultVar( GetConverter().AddVar(lb(), ub(), type()) );
  }
  void AddConstraint() {
    GetConverter().AddConstraint( Constraint (
                                    std::move(args_), ResultVar() ) );
  }
};


/// This is a helper to produce a 'final' FC converter avoiding using Impl
/// Then it could be specialized for individual constraint types
template <class Converter, class Constraint>
class FunctionalConstraintConverter :
    public BasicFunctionalConstraintConverter<
              FunctionalConstraintConverter<Converter, Constraint>, Converter, Constraint > {
  using Base = BasicFunctionalConstraintConverter<
      FunctionalConstraintConverter<Converter, Constraint>, Converter, Constraint >;
public:
  FunctionalConstraintConverter(
      Converter& cvt, typename Converter::VarArray&& aa) : Base(cvt, std::move(aa)) { }
};

template <class Converter, class Constraint, class Converter2>
FunctionalConstraintConverter<Converter, Constraint>
makeFuncConstrConverter(Converter2& cvt, typename Converter::VarArray&& aa) {
  return FunctionalConstraintConverter<Converter, Constraint>(
        static_cast<Converter&>(cvt), std::move(aa) );
}

} // namespace mp

#endif // EXPR2CONSTRAINT_H
