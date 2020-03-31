#ifndef MP2MIP_H
#define MP2MIP_H

#include "mp/convert/expr2constraint.h"
#include "mp/convert/constraints/maximum.h"

namespace mp {

/// MPToMIPConverter: one of the converters requiring a "minimal" output interface
template <class Impl, class Backend,
          class Model = BasicProblem<std::allocator<char> > >
class MPToMIPConverter
    : public MPFlatConverter<Impl, Backend, Model>
{
  using ClassName = MPToMIPConverter<Impl, Backend, Model>;
  using Base = MPFlatConverter<Impl, Backend, Model>;
public:

  void Convert(typename Model::MutCommonExpr e) {
    throw std::runtime_error("MPToMIPConverter: No common exprs convertible yet TODO");
  }

  void Convert(typename Model::MutObjective obj) {
    if (obj.nonlinear_expr())
      throw std::runtime_error("MPToMIPConverter: Only linear objectives allowed TODO");
  }

  void Convert(typename Model::MutAlgebraicCon con) {
    LinearExpr &linear = con.linear_expr();
    if (NumericExpr e = con.nonlinear_expr()) {
      linear.AddTerms(this->Visit(e));
      con.unset_nonlinear_expr();                  // delete the non-linear expr
      auto ne = NumericExpr();
      assert(!ne);
    } // Modifying the original constraint by replacing the expr
  }

  void Convert(typename Model::MutLogicalCon e) {
    throw std::runtime_error("MPToMIPConverter: Only algebraic constraints implemented TODO");
  }

  EExpr VisitMinus(UnaryExpr e) {
    auto ee = this->Visit(e.arg());
    for (auto& term: ee)
      term.set_coef(-term.coef());
    return ee;
  }

  EExpr VisitMax(typename Base::VarArgExpr e) {       // TODO why need Base:: here in g++ 9.2.1?
    auto e2c = makeE2CConverter<Expr2Constr, ClassName, MaximumConstraint<ClassName, Backend> >(*this);
    return e2c.ConvertArray(e);
  }

  EExpr VisitVariable(Reference r) {
    return r.index();
  }


};

} // namespace mp

#endif // MP2MIP_H
