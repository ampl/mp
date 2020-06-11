#ifndef CONVERTER_FLAT_H
#define CONVERTER_FLAT_H

#include <unordered_map>
#include <cmath>

#include "mp/convert/preprocess.h"
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
  BasicMPFlatConverter() {
    InitOptions();
  }

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

  void Convert(typename Model::MutCommonExpr ) {
    /// Convert on demand, not here
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
  int Convert2Var(EExpr&& ee) {
    if (ee.is_variable())
      return ee.get_representing_variable();
    if (ee.is_constant())
      return MakeFixedVar(ee.constant_term());
    PreprocessInfoStd bnt = ComputeBoundsAndType(ee);
    auto r = MP_DISPATCH( AddVar(bnt.lb_, bnt.ub_, bnt.type_) );
    AddConstraint(LinearDefiningConstraint(r, std::move(ee)));
    return r;
  }

  PreprocessInfoStd ComputeBoundsAndType(const AffineExpr& ae) {
    PreprocessInfoStd result;
    result.lb_ = result.ub_ = ae.constant_term();    // TODO reuse bounds if supplied
    result.type_ = is_integer(result.lb_) ? var::INTEGER : var::CONTINUOUS;
    result.linexp_type_ = var::INTEGER;
    auto& model = MP_DISPATCH( GetModel() );
    for (const auto& term: ae) {
      auto v = model.var(term.var_index());
      if (term.coef() >= 0.0) {
        result.lb_ += term.coef() * v.lb();
        result.ub_ += term.coef() * v.ub();
      } else {
        result.lb_ += term.coef() * v.ub();
        result.ub_ += term.coef() * v.lb();
      }
      if (var::INTEGER!=v.type() || !is_integer(term.coef())) {
        result.type_=var::CONTINUOUS;
        result.linexp_type_=var::CONTINUOUS;
      }
    }
    return result;
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
    std::array<EExpr, 2> ee;
    Exprs2EExprs(ea, ee);
    ee[0].Subtract(ee[1]);
    return AssignResultToArguments( FuncConstraint(ee[0]) );
  }

  template <class ExprArray, size_t N>
  void Exprs2EExprs(const ExprArray& ea, std::array<EExpr, N>& result) {
    assert(ea.size() == result.size());
    auto itea = ea.begin();
    for (int i=0; i<N; ++i, ++itea)
      result[i] = MP_DISPATCH( Convert2EExpr(*itea) );
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

  EExpr VisitCommonExpr(Reference r) {
    const auto index = r.index();
    if (index >= common_exprs_.size())
      common_exprs_.resize(index+1, -1);          // init by -1, "no variable"
    if (common_exprs_[index]<0) {                 // not yet converted
      auto ce = MP_DISPATCH( GetModel() ).common_expr(index);
      EExpr eexpr(ce.linear_expr());
      if (ce.nonlinear_expr())
        eexpr.Add( Convert2EExpr(ce.nonlinear_expr()) );
      common_exprs_[index] = Convert2Var(std::move(eexpr));
    }
    return EExpr::Variable{ common_exprs_[index] };
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

  EExpr VisitSum(typename BaseExprVisitor::SumExpr expr) {
    EExpr sum;
    for (typename BaseExprVisitor::SumExpr::iterator i =
         expr.begin(), end = expr.end(); i != end; ++i)
      sum.Add( MP_DISPATCH( Convert2EExpr(*i) ) );
    return sum;
  }

  EExpr VisitMax(typename BaseExprVisitor::VarArgExpr e) {       // TODO why need Base:: here in g++ 9.2.1?
    return VisitFunctionalExpression<MaximumConstraint>(e);
  }

  EExpr VisitMin(typename BaseExprVisitor::VarArgExpr e) {
    return VisitFunctionalExpression<MinimumConstraint>(e);
  }

  EExpr VisitAbs(UnaryExpr e) {
    return VisitFunctionalExpression<AbsConstraint>({ e.arg() });
  }

  EExpr VisitEQ(RelationalExpr e) {
    return VisitRelationalExpression<EQ0Constraint>({ e.lhs(), e.rhs() });
  }

  EExpr VisitNE(RelationalExpr e) {
    auto EQ = this->GetModel().MakeRelational(expr::EQ, e.lhs(), e.rhs());
    return VisitFunctionalExpression<NotConstraint>({ EQ });
  }

  EExpr VisitLE(RelationalExpr e) {
    return VisitRelationalExpression<LE0Constraint>({ e.lhs(), e.rhs() });
  }

  EExpr VisitGE(RelationalExpr e) {
    return VisitRelationalExpression<LE0Constraint>({ e.rhs(), e.lhs() });
  }

  EExpr VisitNot(NotExpr e) {
    return VisitFunctionalExpression<NotConstraint>({ e.arg() });
  }

  EExpr VisitAnd(BinaryLogicalExpr e) {
    return VisitFunctionalExpression<ConjunctionConstraint>({ e.lhs(), e.rhs() });
  }

  EExpr VisitOr(BinaryLogicalExpr e) {
    return VisitFunctionalExpression<DisjunctionConstraint>({ e.lhs(), e.rhs() });
  }

  EExpr VisitIf(IfExpr e) {
    return VisitFunctionalExpression<IfThenConstraint>({
                e.condition(), e.then_expr(), e.else_expr() });
  }

  /////////////// NONLINEAR FUNCTIONS ////////////////
  EExpr VisitPowConstExp(BinaryExpr e) {
    return AssignResultToArguments( PowConstraint(
      PowConstraint::Arguments{ Convert2Var(e.lhs()) },
      PowConstraint::Parameters{ Cast<NumericConstant>(e.rhs()).value() } ) );
  }

  EExpr VisitPow2(UnaryExpr e) {     // MIP could have better conversion for pow2
    return AssignResultToArguments( PowConstraint(
      PowConstraint::Arguments{ Convert2Var(e.arg()) },
      PowConstraint::Parameters{ 2.0 } ) );
  }

  EExpr VisitSqrt(UnaryExpr e) {
    return AssignResultToArguments( PowConstraint(
      PowConstraint::Arguments{ Convert2Var(e.arg()) },
      PowConstraint::Parameters{ 0.5 } ) );
  }

  EExpr VisitExp(UnaryExpr e) {
    return VisitFunctionalExpression<ExpConstraint>({ e.arg() });
  }

  EExpr VisitPowConstBase(BinaryExpr e) {
    return AssignResultToArguments( ExpAConstraint(
      ExpAConstraint::Arguments{ Convert2Var(e.rhs()) },
      ExpAConstraint::Parameters{ Cast<NumericConstant>(e.lhs()).value() } ) );
  }

  EExpr VisitLog(UnaryExpr e) {
    return VisitFunctionalExpression<LogConstraint>({ e.arg() });
  }

  EExpr VisitLog10(UnaryExpr e) {
    return AssignResultToArguments( LogAConstraint(
      LogAConstraint::Arguments{ Convert2Var(e.arg()) },
      LogAConstraint::Parameters{ 10.0 } ) );
  }

  EExpr VisitSin(UnaryExpr e) {
    return VisitFunctionalExpression<SinConstraint>({ e.arg() });
  }

  EExpr VisitCos(UnaryExpr e) {
    return VisitFunctionalExpression<CosConstraint>({ e.arg() });
  }

  EExpr VisitTan(UnaryExpr e) {
    return VisitFunctionalExpression<TanConstraint>({ e.arg() });
  }



public:

  //////////////////////////// CUSTOM CONSTRAINTS CONVERSION ////////////////////////////
  ///
  //////////////////////////// THE CONVERSION LOOP: BREADTH-FIRST ///////////////////////
  void ConvertExtraItems() {
    try {
      for (int endConstraintsThisLoop = 0, endPrevious = 0;
           (endConstraintsThisLoop = this->GetModel().num_custom_cons()) > endPrevious;
           endPrevious = endConstraintsThisLoop
           ) {
        MP_DISPATCH( PreprocessIntermediate() );                        // preprocess before each level
        ConvertExtraItemsInRange(endPrevious, endConstraintsThisLoop);
      }
      MP_DISPATCH( ConvertMaps(); );
      MP_DISPATCH( PreprocessFinal() );                                 // final prepro
    } catch (const ConstraintConversionFailure& cff) {
      throw std::logic_error(cff.message());
    }
  }

  void ConvertExtraItemsInRange(int first, int after_last) {
    for (; first<after_last; ++first) {
      auto* pConstraint = this->GetModel().custom_con(first);
      if (!pConstraint->IsRemoved()) {
        const auto acceptanceLevel =
            pConstraint->BackendAcceptance(this->GetBackend());
        if (NotAccepted == acceptanceLevel) {
          pConstraint->ConvertWith(*this);
          pConstraint->Remove();
        }
        else if (AcceptedButNotRecommended == acceptanceLevel) {
          try {
            pConstraint->ConvertWith(*this);
            pConstraint->Remove();
          } catch (const ConstraintConversionFailure& ccf) {
            MP_DISPATCH( Print(
                           "WARNING: {}. Will pass the constraint "
                           "to the backend {}. Continuing\n",
                           ccf.message(), typeid(Backend).name() ) );
          }
        }
      }
    }
  }

  //////////////////////// WHOLE-MODEL PREPROCESSING /////////////////////////
  void PreprocessIntermediate() { }
  void ConvertMaps() { }
  void PreprocessFinal() { }



  //////////////////////////// CUSTOM CONSTRAINTS ////////////////////////////
  ///
  //////////////////////////// CONSTRAINT PROPAGATORS ///////////////////////////////////


  template <class PreprocessInfo>
  void PreprocessConstraint(
      MinimumConstraint& c, PreprocessInfo& prepro) {
    auto& m = MP_DISPATCH( GetModel() );
    auto& args = c.GetArguments();
    prepro.narrow_result_bounds( m.lb_array(args),
                          m.ub_min_array(args) );
    prepro.set_result_type( m.common_type(args) );
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      MaximumConstraint& c, PreprocessInfo& prepro) {
    auto& m = MP_DISPATCH( GetModel() );
    auto& args = c.GetArguments();
    prepro.narrow_result_bounds( m.lb_max_array(args),
                          m.ub_array(args) );
    prepro.set_result_type( m.common_type(args) );
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      AbsConstraint& c, PreprocessInfo& prepro) {
    prepro.narrow_result_bounds(0.0, this->Infty());
    prepro.set_result_type( var_type(c.GetArguments()[0]) );
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      EQ0Constraint& c, PreprocessInfo& prepro) {
    prepro.narrow_result_bounds(0.0, 1.0);
    prepro.set_result_type( var::INTEGER );
    if (0!=CanPreprocess( options_.preprocessEqualityResultBounds_ ))
      if (FixEqualityResult(c, prepro))
        return;
    PreprocessEqVarConst__unifyCoef(c);
    if (0!=CanPreprocess( options_.preprocessEqualityBvar_ ))
      if (ReuseEqualityBinaryVar(c, prepro))
        return;
  }

  template <class PreprocessInfo>
  bool FixEqualityResult(
      EQ0Constraint& c, PreprocessInfo& prepro) {
    AffineExpr& ae = c.GetArguments();
    if (ae.is_constant()) {                  // const==0
      auto res = (double)int(0.0==ae.constant_term());
      prepro.narrow_result_bounds(res, res);
      return true;
    }
    auto bndsNType = ComputeBoundsAndType(ae);
    if (bndsNType.lb() > 0.0 || bndsNType.ub() < 0.0) {
      prepro.narrow_result_bounds(0.0, 0.0);
      return true;
    }
    if (bndsNType.lb()==0.0 && bndsNType.ub()==0.0) {
      prepro.narrow_result_bounds(1.0, 1.0);
      return true;
    }
    if (var::INTEGER==bndsNType.linexp_type_ &&
        !is_integer(ae.constant_term())) {
      prepro.narrow_result_bounds(0.0, 0.0);
      return true;
    }
    return false;
  }

  static void PreprocessEqVarConst__unifyCoef(EQ0Constraint& c) {
    AffineExpr& ae = c.GetArguments();
    if (1==ae.num_terms()) {
      const double coef = ae.coef(0);
      if (1.0!=coef) {
        assert(0.0!=std::fabs(coef));
        ae.constant_term(ae.constant_term() / coef);
        ae.set_coef(0, 1.0);
      }
    }
  }

  template <class PreprocessInfo>
  bool ReuseEqualityBinaryVar(
      EQ0Constraint& c, PreprocessInfo& prepro) {
    auto& m = MP_DISPATCH( GetModel() );
    AffineExpr& ae = c.GetArguments();
    if (1==ae.num_terms()) {                           // var==const
      assert( 1.0==ae.coef(0) );
      int var = ae.var_index(0);
      if (m.is_binary_var(var)) {            // See if this is binary var==const
        const double rhs = -ae.constant_term();
        if (1.0==rhs)
          prepro.set_result_var( var );
        else if (0.0==rhs)
          prepro.set_result_var( MakeComplementVar(var) );
        else
          prepro.narrow_result_bounds(0.0, 0.0);    // not 0/1 value, result false
        return true;
      }
    }
    return false;
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      LE0Constraint& , PreprocessInfo& prepro) {
    prepro.narrow_result_bounds(0.0, 1.0);
    prepro.set_result_type( var::INTEGER );
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      ConjunctionConstraint& , PreprocessInfo& prepro) {
    prepro.narrow_result_bounds(0.0, 1.0);
    prepro.set_result_type( var::INTEGER );
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      DisjunctionConstraint& , PreprocessInfo& prepro) {
    prepro.narrow_result_bounds(0.0, 1.0);
    prepro.set_result_type( var::INTEGER );
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      NotConstraint& , PreprocessInfo& prepro) {
    prepro.narrow_result_bounds(0.0, 1.0);
    prepro.set_result_type( var::INTEGER );
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      IfThenConstraint& c, PreprocessInfo& prepro) {
    const auto& args = c.GetArguments();
    prepro.narrow_result_bounds(std::min(lb(args[1]), lb(args[2])),
        std::max(ub(args[1]), ub(args[2])));
    prepro.set_result_type( MP_DISPATCH(GetModel()).
                            common_type( { args[1], args[2] } ) );
  }

  ////////////////////// NONLINEAR FUNCTIONS //////////////////////
  template <class PreprocessInfo>
  void PreprocessConstraint(
      ExpConstraint& , PreprocessInfo& prepro) {
    prepro.narrow_result_bounds(0.0, this->Infty());
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      ExpAConstraint& , PreprocessInfo& prepro) {
    prepro.narrow_result_bounds(0.0, this->Infty());
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      LogConstraint& c, PreprocessInfo& ) {
    MP_DISPATCH( GetModel() ).narrow_var_bounds(
          c.GetArguments()[0], 0.0, this->Infty());
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      LogAConstraint& c, PreprocessInfo& ) {
    MP_DISPATCH( GetModel() ).narrow_var_bounds(
          c.GetArguments()[0], 0.0, this->Infty());
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      SinConstraint& , PreprocessInfo& prepro) {
    prepro.narrow_result_bounds(-1.0, 1.0);
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      CosConstraint& , PreprocessInfo& prepro) {
    prepro.narrow_result_bounds(-1.0, 1.0);
  }



  //////////////////////////// CUSTOM CONSTRAINTS //////////////////////
  ///
  //////////////////////////// SPECIFIC CONSTRAINT RESULT-TO-ARGUMENTS PROPAGATORS //////

  void PropagateResult(LinearDefiningConstraint& con, double lb, double ub, Context ctx) {
    con.AddContext(ctx);
  }

  void PropagateResult(NotConstraint& con, double lb, double ub, Context ctx) {
    con.AddContext(ctx);
    for (const auto a: con.GetArguments())
      PropagateResult(a, 1.0-ub, 1.0-lb, -ctx);
  }

  void PropagateResult(ConjunctionConstraint& con, double lb, double ub, Context ctx) {
    con.AddContext(ctx);
    for (const auto a: con.GetArguments())
      PropagateResult(a, lb, 1.0, +ctx);
  }

  void PropagateResult(DisjunctionConstraint& con, double lb, double ub, Context ctx) {
    con.AddContext(ctx);
    for (const auto a: con.GetArguments())
      PropagateResult(a, 0.0, ub, +ctx);
  }

  void PropagateResult(LE0Constraint& con, double lb, double ub, Context ctx) {
    con.AddContext(ctx);
  }

  void PropagateResult(EQ0Constraint& con, double lb, double ub, Context ctx) {
    con.AddContext(ctx);
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
  template <class Constraint>
  void AddConstraint(Constraint&& con) {
    const auto pck = makeConstraintKeeper<Impl, Constraint>(std::forward<Constraint>(con));
    AddConstraintAndTryNoteResultVariable(pck);
  }
  template <class ConstraintKeeper>
  void AddConstraintAndTryNoteResultVariable(ConstraintKeeper* pbc) {
    MP_DISPATCH( GetModel() ).AddConstraint(pbc);
    const auto resvar = pbc->GetResultVar();
    if (resvar>=0)
      AddInitExpression(resvar, pbc);
    if (! MP_DISPATCH( MapInsert(pbc) ))
      throw std::logic_error("Trying to map_insert() duplicated constraint: " +
                             pbc->GetDescription());
  }


  //////////////////////////// UTILITIES /////////////////////////////////
  ///

private:
  std::unordered_map<double, int> map_fixed_vars_;

  std::vector<int> common_exprs_;               // variables equal to the result

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
  var::Type var_type(int var) const { return this->GetModel().var(var).type(); }
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


  ///////////////////////////////////////////////////////////////////////
  /////////////////////// OPTIONS /////////////////////////
  ///
private:
  struct Options {
    int preprocessAnything_ = 1;
    int preprocessEqualityResultBounds_ = 1;
    int preprocessEqualityBvar_ = 1;
  };
  Options options_;

  void InitOptions() {
    this->AddOption("cvt:prepro:all",
        "0/1*: Set to 0 to disable all presolve in the converter",
        options_.preprocessAnything_);
    this->AddOption("cvt:prepro:eqresult",
        "0/1*: Preprocess reified equality comparison's boolean result bounds",
        options_.preprocessEqualityResultBounds_);
    this->AddOption("cvt:prepro:eqbinary",
        "0/1*: Preprocess reified equality comparison with a binary variable",
        options_.preprocessEqualityBvar_);
  }

protected:
  bool CanPreprocess(int f) const {
    return 0!=options_.preprocessAnything_ && 0!=f;
  }

};


} // namespace mp

#endif // CONVERTER_FLAT_H
