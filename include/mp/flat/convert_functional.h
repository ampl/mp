#ifndef EXPR2CONSTRAINT_H
#define EXPR2CONSTRAINT_H

#include <utility>
#include <cassert>

#include "mp/flat/preprocess.h"

namespace mp {

/// Helper class providing a default framework for assigning result
/// to a functional expression,
/// possibly adding a corresponding constraint on the result variable
template <class Impl, class Converter, class Constraint>
class BasicFCC {
  Converter& converter_;
  Constraint constr_;
public:
  PreprocessInfo<Constraint> prepro_;
  using Var = typename Converter::Var;
protected:
  Converter& GetConverter() { return converter_; }
  Constraint& GetConstraint() { return constr_; }
  double lb() const { return prepro_.lb_; }
  double ub() const { return prepro_.ub_; }
  var::Type type() const { return prepro_.type_; }
protected:
  bool ResultIsConstant() const { return prepro_.is_constant(); }
  bool ResultVarIsKnown() const { return prepro_.is_result_var_known(); }
  bool MapFind() {
    const auto pck = GetConverter().MapFind(GetConstraint());
    if (pck) {
      SetResultVar(pck.GetResultVar());
      return true;
    }
    return false;
  }
  int GetResultVar() const { return prepro_.get_result_var(); }
protected:
  void SetResultVar(int r) { prepro_.set_result_var(r); }
public:
  BasicFCC(Converter& cvt, Constraint&& fc) noexcept :
    converter_(cvt), constr_(std::move(fc)) { }
  /// Holder for the conversion result of an expression
  class VarOrConst {
    const bool is_v_;
    union {
      double c_;
      Var var_;
    };
    VarOrConst(bool isv, Var v) :
      is_v_(isv), var_(v) { assert(is_var()); }
    VarOrConst(bool isv, double c) :
      is_v_(isv), c_(c) { assert(is_const()); }
  public:
    bool is_var() const { return is_v_; }
    bool is_const() const { return !is_var(); }
    double get_const() const { assert(is_const()); return c_; }
    Var get_var() const { assert(is_var()); return var_; }
    static VarOrConst MakeVar(Var v)
    { assert(Converter::VoidVar()!=v); return VarOrConst(true, v); }
    static VarOrConst MakeConst(double c)
    { return VarOrConst(false, c); }
  };
  /// Convert array of arguments into a result (var or const),
  /// possibly adding extra constraint(s).
  /// @return either a constant or a variable
  VarOrConst Convert() {
    MP_DISPATCH( PreprocessArguments() );
    if (ResultIsConstant())
      return VarOrConst::MakeConst( lb() );
    if (ResultVarIsKnown())
      return VarOrConst::MakeVar( GetResultVar() );
    if (MapFind())
      return VarOrConst::MakeVar( GetResultVar() );
    MP_DISPATCH( AddResultVariable() );
    MP_DISPATCH( AddConstraint() );
    return VarOrConst::MakeVar( GetResultVar() );
  }
  void PreprocessArguments() {
    GetConverter().PreprocessConstraint(GetConstraint(), prepro_);
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
  FCC(Converter& cvt, Constraint&& fc) noexcept : Base(cvt, std::move(fc)) { }
};

template <class Converter, class Constraint, class Converter2>
FCC<Converter, Constraint>
MakeFuncConstrConverter(Converter2& cvt, Constraint&& fc) {
  return FCC<Converter, Constraint>(
        static_cast<Converter&>(cvt), std::forward<Constraint>(fc) );
}


//////////////////////////// SPECIALIZED FCCs and/or their components /////////////////////////////
///
///////////////////////////////////////////////////////////////////////////////////////////////////


} // namespace mp

#endif // EXPR2CONSTRAINT_H
