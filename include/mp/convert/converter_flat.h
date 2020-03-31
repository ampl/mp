#ifndef CONVERTER_FLAT_H
#define CONVERTER_FLAT_H

#include "mp/converter.h"
#include "mp/expr-visitor.h"
#include "mp/convert/expr2constraint.h"
#include "mp/convert/model.h"

#include "mp/convert/constraints/maximum.h"

namespace mp {


/// BasicMPFlatConverter: it "flattens" most expressions by replacing them by a result variable and constraints
/// Such constraints might need to be decomposed, which is handled by redefined virtual methods in derived classes
template <class Impl, class Backend,
          class Model = BasicModel<std::allocator<char> > >
class BasicMPFlatConverter
    : public BasicMPConverter<Impl, Backend, Model>,
      public ExprVisitor<Impl, EExpr>,
      public BasicConstraintConverter
{
public:

  using ClassName = BasicMPFlatConverter<Impl, Backend, Model>;
  using BaseExprVisitor = ExprVisitor<Impl, EExpr>;
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

  EExpr VisitMax(typename BaseExprVisitor::VarArgExpr e) {       // TODO why need Base:: here in g++ 9.2.1?
    auto e2c = makeE2CConverter<Expr2Constr, Impl, MaximumConstraint<Impl, Backend> >(*this);
    return e2c.ConvertArray(e);
  }

  EExpr VisitVariable(Reference r) {
    return r.index();
  }

};


} // namespace mp

#endif // CONVERTER_FLAT_H
