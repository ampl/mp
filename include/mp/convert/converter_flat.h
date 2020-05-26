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

/// BasicMPFlatConverter: it "flattens" most expressions
/// by replacing them by a result variable and constraints.
/// Such constraints might need to be decomposed, which is
/// handled by overloaded methods in derived classes
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


public:
  //////////////////////////// CONVERTERS OF STANDRAD MP ITEMS //////////////////////////////
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////

  void Convert(typename Model::MutCommonExpr e) {
    throw std::runtime_error("BasicMPFlatConverter: No common exprs convertible yet TODO");
  }

  void Convert(typename Model::MutObjective obj) {
    if (NumericExpr e = obj.nonlinear_expr()) {
      LinearExpr &linear = obj.linear_expr();
      const auto affine_expr=MP_DISPATCH( Visit(e) );
      linear.AddTerms(affine_expr);
      if (std::fabs(affine_expr.constant_term())!=0.0) {
        linear.AddTerm(MakeFixedVar(affine_expr.constant_term()), 1.0);
      }
      obj.unset_nonlinear_expr();
    } // Modifying the original objective by replacing the expr
  }

  void Convert(typename Model::MutAlgebraicCon con) {
    if (NumericExpr e = con.nonlinear_expr()) {
      LinearExpr &linear = con.linear_expr();
      const auto affine_expr=MP_DISPATCH( Visit(e) );
      linear.AddTerms(affine_expr);
      con.set_lb(con.lb() + affine_expr.constant_term());
      con.set_ub(con.ub() + affine_expr.constant_term());
      con.unset_nonlinear_expr();                  // delete the non-linear expr
    } // Modifying the original constraint by replacing the expr
  }

  void Convert(typename Model::MutLogicalCon e) {
    const auto resvar = MP_DISPATCH( Convert2Var(e.expr()) );
    PropagateResult(resvar, 1.0, 1.0, +Context());
  }

  void PropagateResult(int var, double lb, double ub, Context ctx) {
    this->GetModel().narrow_var_bounds(var, lb, ub);
    if (HasInitExpression(var))
      GetInitExpression(var)->PropagateResult(*this, lb, ub, ctx);
  }


public:
  //////////////////////////////////// VISITOR ADAPTERS /////////////////////////////////////////

  /// Convert an expression to an EExpr
  EExpr Convert2EExpr(Expr e) {
    return MP_DISPATCH(Visit(e));
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
    auto r = MP_DISPATCH( AddVar(bnt.lb_, bnt.ub_, bnt.type_) );
    auto lck = makeConstraint<Impl, LinearDefiningConstraint>(r, std::move(ee));
    AddConstraint(lck);
    return r;
  }

  /// Generic functional expression array visitor
  /// Can produce a new variable/expression and specified constraints on it
  template <class FuncConstraint, class ExprArray=std::initializer_list<Expr> >
  EExpr VisitFunctionalExpression(ExprArray ea) {
    FuncConstraint fc;
    Exprs2Vars(ea, fc.GetArguments());
    return AssignResultToArguments( std::move(fc) );
  }

  template <class ExprArray>
  void Exprs2Vars(const ExprArray& ea, std::vector<int>& result) {
    assert(result.empty());
    result.reserve(ea.num_args());
    for (const auto& e: ea)
      result.push_back( MP_DISPATCH( Convert2Var(e) ) );
  }

  template <class Expr>
  void Exprs2Vars(const std::initializer_list<Expr>& ea, std::vector<int>& result) {
    assert(result.empty());
    result.reserve(ea.size());
    for (const auto& e: ea)
      result.push_back( MP_DISPATCH( Convert2Var(e) ) );
  }

  template <class ExprArray, size_t N>
  void Exprs2Vars(const ExprArray& ea, std::array<int, N>& result) {
    assert(ea.size() == result.size());
    auto itea = ea.begin();
    for (int i=0; i<N; ++i, ++itea)
      result[i] = MP_DISPATCH( Convert2Var(*itea) );
  }

  template <class FuncConstraint>
  EExpr AssignResultToArguments(FuncConstraint&& fc) {
    auto fcc = MakeFuncConstrConverter<Impl, FuncConstraint>(
          *this, std::forward<FuncConstraint>(fc));
    return fcc.Convert();
  }

  /// Generic relational expression visitor
  /// Can produce a new variable/expression and specified constraints on it
  template <class FuncConstraint, class ExprArray=std::initializer_list<Expr> >
  EExpr VisitRelationalExpression(ExprArray ea) {
    return VisitFunctionalExpression<FuncConstraint>(ea);
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

  EExpr VisitEQ(RelationalExpr e) {
    return VisitRelationalExpression<EQConstraint>({ e.lhs(), e.rhs() });
  }

  EExpr VisitNE(RelationalExpr e) {
    auto EQ = this->GetModel().MakeRelational(expr::EQ, e.lhs(), e.rhs());
    return VisitFunctionalExpression<NotConstraint>({ EQ });
  }

  EExpr VisitLE(RelationalExpr e) {
    return VisitRelationalExpression<LEConstraint>({ e.lhs(), e.rhs() });
  }

  EExpr VisitNot(NotExpr e) {
    return VisitFunctionalExpression<NotConstraint>({ e.arg() });
  }

  EExpr VisitOr(BinaryLogicalExpr e) {
    return VisitFunctionalExpression<DisjunctionConstraint>({ e.lhs(), e.rhs() });
  }


public:

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
        if (Recommended !=
            pConstraint->BackendAcceptance(this->GetBackend())) {
          pConstraint->ConvertWith(*this);
          pConstraint->Remove();
        }
      }
    }
  }

  //////////////////////// WHOLE-MODEL PREPROCESSING /////////////////////////
  void PreprocessIntermediate() { }
  void PreprocessFinal() { }



  //////////////////////////// CUSTOM CONSTRAINTS ////////////////////////////
  ///
  //////////////////////////// CONSTRAINT PROPAGATORS ///////////////////////////////////


  /// Preprocess minimum
  void PreprocessConstraint(
      MinimumConstraint& c, PreprocessInfo<MinimumConstraint>& prepro) {
    auto& m = MP_DISPATCH( GetModel() );
    auto& args = c.GetArguments();
    prepro.narrow_result_bounds( m.lb_array(args),
                          m.ub_min_array(args) );
    prepro.set_result_type( m.common_type(args) );
  }

  /// Preprocess maximum
  void PreprocessConstraint(
      MaximumConstraint& c, PreprocessInfo<MaximumConstraint>& prepro) {
    auto& m = MP_DISPATCH( GetModel() );
    auto& args = c.GetArguments();
    prepro.narrow_result_bounds( m.lb_max_array(args),
                          m.ub_array(args) );
    prepro.set_result_type( m.common_type(args) );
  }

  /// Preprocess EQ
  void PreprocessConstraint(
      EQConstraint& c, PreprocessInfo<EQConstraint>& prepro) {
    auto& m = MP_DISPATCH( GetModel() );
    auto& args = c.GetArguments();
    prepro.narrow_result_bounds(0.0, 1.0);
    prepro.set_result_type( var::INTEGER );
    if (m.is_fixed(args[0]) && m.is_fixed(args[1])) {
      auto res = (double)int(m.fixed_value(args[0])==m.fixed_value(args[1]));
      prepro.narrow_result_bounds(res, res);
      return;
    }
    if (m.is_fixed(args[0])) {                 // Constant on the right
      std::swap(args[0], args[1]);
    }
    if (m.is_fixed(args[1])) {                 // See if this is binary var==const
      if (m.is_binary_var(args[0])) {
        if (1.0==std::fabs(m.fixed_value(args[1])))
          prepro.set_result_var( args[0] );
        else if (0.0==m.fixed_value(args[1]))
          prepro.set_result_var( MakeComplementVar(args[0]) );
        else
          prepro.narrow_result_bounds(0.0, 0.0);    // not 0/1 value, result false
        return;
      }
    }
  }

  /// Preprocess NE
  void PreprocessConstraint(
      LEConstraint& c, PreprocessInfo<LEConstraint>& prepro) {
    prepro.narrow_result_bounds(0.0, 1.0);
    prepro.set_result_type( var::INTEGER );
  }

  /// Preprocess Disjunction
  void PreprocessConstraint(
      DisjunctionConstraint& c, PreprocessInfo<DisjunctionConstraint>& prepro) {
    prepro.narrow_result_bounds(0.0, 1.0);
    prepro.set_result_type( var::INTEGER );
  }

  /// Preprocess Not
  template <class Converter>
  void PreprocessConstraint(
      NotConstraint& c, PreprocessInfo<NotConstraint>& prepro) {
    prepro.narrow_result_bounds(0.0, 1.0);
    prepro.set_result_type( var::INTEGER );
  }


  //////////////////////////// CUSTOM CONSTRAINTS //////////////////////
  ///
  //////////////////////////// SPECIFIC CONSTRAINT RESULT-TO-ARGUMENTS PROPAGATORS //////

  void PropagateResult(LinearDefiningConstraint& con, double lb, double ub, Context ctx) {
  }

  void PropagateResult(NotConstraint& con, double lb, double ub, Context ctx) {
    con.AddContext(ctx);
    for (const auto a: con.GetArguments())
      PropagateResult(a, 1.0-ub, 1.0-lb, -ctx);
  }

  void PropagateResult(DisjunctionConstraint& con, double lb, double ub, Context ctx) {
    con.AddContext(ctx);
    for (const auto a: con.GetArguments())
      PropagateResult(a, 0.0, ub, +ctx);
  }

  void PropagateResult(LEConstraint& con, double lb, double ub, Context ctx) {
    con.AddContext(ctx);
    const auto& args = con.GetArguments();
    PropagateResult(args[0], this->MinusInfty(), this->Infty(), -ctx);
    PropagateResult(args[1], this->MinusInfty(), this->Infty(), +ctx);
  }


  //////////////////////////// CUSTOM CONSTRAINTS CONVERSION ////////////////////////////
  ///
  //////////////////////////// SPECIFIC CONSTRAINT CONVERTERS ///////////////////////////

  USE_BASE_CONSTRAINT_CONVERTERS(BasicConstraintConverter)      // reuse default converters

  /// If backend does not like LDC, we can redefine it
  void Convert(const LinearDefiningConstraint& ldc) {
    this->AddConstraint(ldc.to_linear_constraint());
  }

public:
  //////////////////////// ADD CUSTOM CONSTRAINT ///////////////////////
  //////////////////////// Takes ownership /////////////////////////////
  void AddConstraint(BasicConstraintKeeper* pbc) {
    MP_DISPATCH( GetModel() ).AddConstraint(pbc);
    const auto resvar = pbc->GetResultVar();
    if (resvar>=0)
      AddInitExpression(resvar, pbc);
  }
  template <class Constraint>
  void AddConstraint(Constraint&& con) {
    AddConstraint(makeConstraint<Impl, Constraint>(std::forward<Constraint>(con)));
  }


  //////////////////////////// UTILITIES /////////////////////////////////
  ///

private:
  std::unordered_map<double, int> map_fixed_vars_;

public:

  //////////////////////////// CREATE OR FIND A FIXED VARIABLE //////////////////////////////
  int MakeFixedVar(double value) {
    auto it = map_fixed_vars_.find(value);
    if (map_fixed_vars_.end()!=it)
      return it->second;
    auto v = BaseConverter::AddVar(value, value);
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
  bool is_fixed(int var) const { return this->GetModel().is_fixed(var); }
  double fixed_value(int var) const { return this->GetModel().fixed_value(var); }

  int MakeComplementVar(int bvar) {
    if (! (lb(bvar)==0.0 && ub(bvar)==1.0) )
      throw std::logic_error("Asked to complement variable with bounds "
                             + std::to_string(lb(bvar)) + ".." + std::to_string(ub(bvar)));
    AffineExpr ae({-1.0}, {bvar}, 1.0);
    return MP_DISPATCH( Convert2Var(std::move(ae)) );
  }

  struct VarInfo {
    BasicConstraintKeeper *pInitExpr=nullptr;
  };

private:
  std::vector<VarInfo> var_info_;

public:
  void AddInitExpression(int var, BasicConstraintKeeper* pie) {
    var_info_.resize(std::max(var_info_.size(), (size_t)var+1));
    var_info_[var].pInitExpr = pie;
  }

  bool HasInitExpression(int var) const {
    return var_info_.size()>var && nullptr!=var_info_[var].pInitExpr;
  }

  BasicConstraintKeeper* GetInitExpression(int var) {
    assert(HasInitExpression(var));
    return var_info_[var].pInitExpr;
  }


};


} // namespace mp

#endif // CONVERTER_FLAT_H
