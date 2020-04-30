#ifndef CONVERTER_FLAT_H
#define CONVERTER_FLAT_H

#include <unordered_map>
#include <cmath>

#include "mp/convert/basic_converters.h"
#include "mp/expr-visitor.h"
#include "mp/convert/convert_functional.h"
#include "mp/convert/model.h"
#include "mp/convert/std_constr.h"

namespace mp {

/// Result expression type for expression conversions
class EExpr : public AffineExpr {
public:
  EExpr() = default;
  EExpr(AffineExpr&& ae) : AffineExpr(std::move(ae)) { }
  EExpr(Constant c) : AffineExpr(c) {}
  EExpr(Variable v) : AffineExpr(v) {}
  EExpr(int i, double c) { AddTerm(i, c); }
};

/// BasicMPFlatConverter: it "flattens" most expressions by replacing them by a result variable and constraints
/// Such constraints might need to be decomposed, which is handled by overloaded methods in derived classes
template <class Impl, class Backend,
          class Model = BasicModel<std::allocator<char> > >
class BasicMPFlatConverter
    : public BasicMPConverter<Impl, Backend, Model>,
      public ExprVisitor<Impl, EExpr>,
      public BasicConstraintConverter
{
public:
  using EExprType = EExpr;
  using VarArray = std::vector<int>;

protected:
  using ClassName = BasicMPFlatConverter<Impl, Backend, Model>;
  using BaseConverter = BasicMPConverter<Impl, Backend, Model>;
  using BaseExprVisitor = ExprVisitor<Impl, EExpr>;

  using EExprArray = std::vector<EExpr>;

private:
  std::unordered_map<double, int> map_fixed_vars_;

  //////////////////////////// UTILITIES /////////////////////////////////
public:

  //////////////////////////// CREATE OR FIND A FIXED VARIABLE //////////////////////////////
  int MakeFixedVar(double value) {
    auto it = map_fixed_vars_.find(value);
    if (map_fixed_vars_.end()!=it)
      return it->second;
    auto v = this->AddVar(value, value);
    map_fixed_vars_[value] = v;
    return v;
  }

  /// Create or find a fixed variable
  int AddVar(double lb, double ub, var::Type type = var::CONTINUOUS) {
    if (lb!=ub)
      return BaseConverter::AddVar(lb, ub, type);
    return MakeFixedVar(lb);
  }

  double lb(int var) const { return this->GetModel().var(var).lb(); }
  double ub(int var) const { return this->GetModel().var(var).ub(); }

  int MakeComplementVar(int bvar) {
    if (! (lb(bvar)==0.0 && ub(bvar)==1.0) )
      throw std::logic_error("Asked to complement variable with bounds "
                             + std::to_string(lb(bvar)) + ".." + std::to_string(ub(bvar)));
    AffineExpr ae({-1.0}, {bvar}, 1.0);
    return MP_DISPATCH( Convert2Var(std::move(ae)) );
  }

  //////////////////////////// CONVERTERS OF STANDRAD MP ITEMS //////////////////////////////
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

  //////////////////////// MODEL PREPROCESSING /////////////////////////
  void PreprocessIntermediate() { }
  void PreprocessFinal() { }

public:
  //////////////////////// ADD CUSTOM CONSTRAINT ///////////////////////
  //////////////////////// Takes ownership /////////////////////////////
  void AddConstraint(BasicConstraintKeeper* pbc) {
    MP_DISPATCH( GetModel() ).AddConstraint(pbc);
  }
  template <class Constraint>
  void AddConstraint(Constraint&& con) {
    AddConstraint(makeConstraint<Impl, Constraint>(std::forward<Constraint>(con)));
  }

public:
  //////////////////////////////////// VISITOR ADAPTERS /////////////////////////////////////////

  /// Convert an expression to an EExpr
  EExpr Convert2EExpr(Expr e) {
    return this->Visit(e);
  }

  /// From an expression:
  /// Adds a result variable r and constraint r == expr
  int Convert2Var(Expr e) {
    return Convert2Var( Convert2EExpr(e) );
  }
  int Convert2Var(EExpr ee) {
    if (ee.is_variable())
      return ee.get_representing_variable();
    if (ee.is_constant())
      return MakeFixedVar(ee.constant_term());
    PreprocessInfoStd bnt;
    ComputeBoundsAndType(this->GetModel(), ee, bnt);
    auto r = this->AddVar(bnt.lb_, bnt.ub_, bnt.type_);
    auto lck = makeConstraint<Impl, LinearDefiningConstraint>(std::move(ee), r);
    AddConstraint(lck);
    return r;
  }

  /// Generic functional expression array visitor
  /// Can produce a new variable/expression and specified constraints on it
  template <class FuncConstraint, class ExprArray>
  EExpr VisitFunctionalExpression(ExprArray ea) {
    FuncConstraint fc;
    Exprs2Vars(ea, fc.GetArguments());
    return AssignResultToArguments( std::move(fc) );
  }

  template <class ExprArray, class Vars>
  void Exprs2Vars(const ExprArray& ea, Vars& result) {
    assert(result.empty());
    result.reserve(ea.num_args());
    for (const auto& e: ea)
      result.push_back( MP_DISPATCH( Convert2Var(e) ) );
  }

  template <class FuncConstraint>
  EExpr AssignResultToArguments(FuncConstraint&& fc) {
    auto fcc = MakeFuncConstrConverter<Impl, FuncConstraint>(
          *this, std::forward<FuncConstraint>(fc));
    return fcc.Convert();
  }


  ///////////////////////////////// EXPRESSION VISITORS ////////////////////////////////////
  ///
  //////////////////////////////////////////////////////////////////////////////////////////

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
    return VisitFunctionalExpression<MaximumConstraint>(e);
  }

  EExpr VisitMin(typename BaseExprVisitor::VarArgExpr e) {
    return VisitFunctionalExpression<MinimumConstraint>(e);
  }

};


} // namespace mp

#endif // CONVERTER_FLAT_H
