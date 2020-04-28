#ifndef EXPR2CONSTRAINT_H
#define EXPR2CONSTRAINT_H

#include "mp/convert/prepro_args.h"

namespace mp {

/// Helper class providing a default framework for assigning result
/// to a functional constraint,
/// possibly adding a constraint on the result variable
template <class Impl, class Converter, class Constraint>
class BasicFCC {
  Converter& converter_;
  Constraint constr_;
public:
  using EExprType = typename Converter::EExprType;
  BasicPreprocessInfo<Constraint> prepro_;
protected:
  Converter& GetConverter() { return converter_; }
  Constraint& GetConstraint() { return constr_; }
  double lb() const { return prepro_.lb_; }
  double ub() const { return prepro_.ub_; }
  var::Type type() const { return prepro_.type_; }
protected:
  bool ResultIsConstant() const { return prepro_.is_constant(); }
  bool ResultVarIsKnown() const { return prepro_.is_result_var_known(); }
  int GetResultVar() const { return prepro_.get_result_var(); }
protected:
  void SetResultVar(int r) { prepro_.set_result_var(r); }
public:
  BasicFCC(Converter& cvt, Constraint&& fc) :
    converter_(cvt), constr_(std::move(fc)) { }
  /// Convert array of arguments into a result expression
  /// possible adding extra constraint(s)
  EExprType Convert() {
    MP_DISPATCH( PreprocessArguments() );
    if (ResultIsConstant())
      return typename EExprType::Constant{ lb() };
    if (ResultVarIsKnown())
      return typename EExprType::Variable{ GetResultVar() };
    MP_DISPATCH( AddResultVariable() );
    MP_DISPATCH( AddConstraint() );
    return typename EExprType::Variable{ GetResultVar() };
  }
  void PreprocessArguments() {
    PreprocessConstraint(GetConverter().GetModel(), GetConstraint(), prepro_);
  }
  void AddResultVariable() {
    auto r = GetConverter().AddVar(lb(), ub(), type());
    SetResultVar( r );
    GetConstraint().SetResultVar( r );
  }
  void AddConstraint() {
    GetConverter().AddConstraint( std::move(GetConstraint()) );
  }
};

/// This is a helper to produce a 'final' FC converter avoiding using Impl
/// Then it could be specialized for individual constraint types
template <class Converter, class Constraint>
class FCC : public BasicFCC< FCC<Converter, Constraint>, Converter, Constraint > {
  using Base = BasicFCC< FCC<Converter, Constraint>, Converter, Constraint >;
public:
  FCC(Converter& cvt, Constraint&& fc) : Base(cvt, std::move(fc)) { }
};

template <class Converter, class Constraint, class Converter2>
FCC<Converter, Constraint>
MakeFuncConstrConverter(Converter2& cvt, Constraint&& fc) {
  return FCC<Converter, Constraint>(
        static_cast<Converter&>(cvt), std::move(fc) );
}


//////////////////////////// SPECIALIZED FCCs and/or their components /////////////////////////////
///
///////////////////////////////////////////////////////////////////////////////////////////////////


} // namespace mp

#endif // EXPR2CONSTRAINT_H
