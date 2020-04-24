#ifndef CONVERTER_FLAT_H
#define CONVERTER_FLAT_H

#include <unordered_map>

#include "mp/convert/basic_converters.h"
#include "mp/expr-visitor.h"
#include "mp/convert/expr2constraint.h"
#include "mp/convert/model.h"

#include "mp/convert/std_constr.h"

namespace mp {


/// BasicMPFlatConverter: it "flattens" most expressions by replacing them by a result variable and constraints
/// Such constraints might need to be decomposed, which is handled by overloaded methods in derived classes
template <class Impl, class Backend,
          class Model = BasicModel<std::allocator<char> > >
class BasicMPFlatConverter
    : public BasicMPConverter<Impl, Backend, Model>,
      public ExprVisitor<Impl, EExpr>,
      public BasicConstraintConverter
{
  std::unordered_map<double, int> map_fixed_vars_;

protected:

  using ClassName = BasicMPFlatConverter<Impl, Backend, Model>;
  using BaseExprVisitor = ExprVisitor<Impl, EExpr>;

  int MakeFixedVar(double value) {
    auto it = map_fixed_vars_.find(value);
    if (map_fixed_vars_.end()!=it)
      return it->second;
    auto v = this->AddVar(value, value);
    map_fixed_vars_[value] = v;
    return v;
  }

public:

  //////////////////////////// COVERTERS OF STANDRAD MP ITEMS ///////////////////////////////
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////
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
    } // Modifying the original constraint by replacing the expr
  }

  void Convert(typename Model::MutLogicalCon e) {
    throw std::runtime_error("MPToMIPConverter: Only algebraic constraints implemented TODO");
  }


  //////////////////////////// CUSTOM CONSTRAINTS CONVERSION ////////////////////////////
  ///
  //////////////////////////// THE CONVERSION LOOP: BREADTH-FIRST ///////////////////////
  void ConvertExtraItems() {
    for (int endConstraintsThisLoop = 0, endPrevious = 0;
         (endConstraintsThisLoop = this->GetModel().num_custom_cons()) > endPrevious;
         endPrevious = endConstraintsThisLoop
         ) {
      PreprocessIntermediate();                        // preprocess before each level
      ConvertExtraItemsInRange(endPrevious, endConstraintsThisLoop);
    }
    PreprocessFinal();                                 // final prepro
  }

  void ConvertExtraItemsInRange(int first, int after_last) {
    for (; first<after_last; ++first) {
      auto* pConstraint = this->GetModel().custom_con(first);
      if (!pConstraint->IsRemoved()) {
        if (BasicConstraintAdder::Recommended !=
            pConstraint->BackendAcceptance(this->GetBackend())) {
          pConstraint->ConvertWith(*this);
          pConstraint->Remove();
        }
      }
    }
  }

  //////////////////////////// CUSTOM CONSTRAINTS CONVERSION ////////////////////////////
  ///
  //////////////////////////// SPECIFIC CONSTRAINT CONVERTERS ///////////////////////////

  USE_BASE_CONSTRAINT_CONVERTERS(BasicConstraintConverter)      // reuse default converters

  /// If backend does not like LDC, we can redefine it
  void Convert(const LinearDefiningConstraint& ldc) {
    this->AddConstraint(ldc.to_linear_constraint());
  }

  //////////////////////// PREPROCESSING /////////////////////////
  void PreprocessIntermediate() { }
  void PreprocessFinal() { }

public:
  /// Add custom constraint. Takes ownership
  void AddConstraint(BasicConstraintKeeper* pbc) {
    MP_DISPATCH( GetModel() ).AddConstraint(pbc);
  }
  template <class Constraint>
  void AddConstraint(Constraint&& con) {
    AddConstraint(makeConstraint<Impl, Constraint>(std::move(con)));
  }

public:
  //////////////////////////////////// Visitor Adapters /////////////////////////////////////////

  EExpr Convert2EExpr(Expr e) {
    return this->Visit(e);
  }

  /// Adds a result variable r and constraint r = expr
  int Convert2Var(Expr e) {
    return Convert2Var( Convert2EExpr(e) );
  }
  int Convert2Var(EExpr ee) {
    if (ee.is_variable())
      return ee.get_representing_variable();
    if (ee.is_constant())
      return MakeFixedVar(ee.constant_term());
    auto r = this->AddVar();    // TODO use a helper class for propagations etc
    auto lck = makeConstraint<Impl, LinearDefiningConstraint>(std::move(ee), r);
    AddConstraint(lck);
    return r;
  }

  //////////////////////////////////// Specialized Visitors /////////////////////////////////////////

  EExpr VisitNumericConstant(NumericConstant n) {
    return EExpr::Constant{ n.value() };
  }

  EExpr VisitVariable(Reference r) {
    return EExpr::Variable{ r.index() };
  }

  EExpr VisitMinus(UnaryExpr e) {
    auto ee = Convert2EExpr(e.arg());
    ee.Negate();
    return ee;
  }

  EExpr VisitAdd(BinaryExpr e) {
    auto ee = Convert2EExpr(e.lhs());
    ee.Add( Convert2EExpr(e.rhs()) );
    return ee;
  }

  EExpr VisitSub(BinaryExpr e) {
    auto el = Convert2EExpr(e.lhs());
    auto er = Convert2EExpr(e.rhs());
    er.Negate();
    el.Add(er);
    return el;
  }

  EExpr VisitMax(typename BaseExprVisitor::VarArgExpr e) {       // TODO why need Base:: here in g++ 9.2.1?
    auto e2c = makeE2CConverter<Expr2Constr, Impl, MaximumConstraint>(*this);
    return e2c.ConvertArray(e);
  }

  EExpr VisitMin(typename BaseExprVisitor::VarArgExpr e) {
    auto e2c = makeE2CConverter<Expr2Constr, Impl, MinimumConstraint>(*this);
    return e2c.ConvertArray(e);
  }

};


} // namespace mp

#endif // CONVERTER_FLAT_H
