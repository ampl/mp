#ifndef CONVERTER_FLAT_H
#define CONVERTER_FLAT_H

#include <unordered_map>
#include <map>
#include <cmath>

#include "mp/convert/preprocess.h"
#include "mp/convert/basic_converters.h"
#include "mp/convert/converter_flat_query.h"
#include "mp/expr-visitor.h"
#include "mp/convert/eexpr.h"
#include "mp/convert/convert_functional.h"
#include "mp/convert/model.h"
#include "mp/convert/std_constr.h"

namespace mp {

/// BasicMPFlatConverter: it "flattens" most expressions
/// by replacing them by a result variable and constraints.
/// Such constraints might need to be converted to others, which is
/// handled by overloaded methods in derived classes
template <class Impl, class Backend,
          class Model = BasicModel< > >
class BasicMPFlatConverter
    : public BasicMPConverter<Impl, Backend, Model>,
      public ExprVisitor<Impl, EExpr>
{
public:
  using BackendType = Backend;
  using ModelType = Model;

public:
  using EExprType = EExpr;
  using VarArray = std::vector<int>;
  template <class Constraint>
    using ConstraintKeeperType = ConstraintKeeper<Impl, Backend, Constraint>;

protected:
  using ClassName = BasicMPFlatConverter<Impl, Backend, Model>;
  using BaseConverter = BasicMPConverter<Impl, Backend, Model>;
  using BaseExprVisitor = ExprVisitor<Impl, EExpr>;

  using EExprArray = std::vector<EExpr>;


public:
  static const char* GetConverterName() { return "BasicMPFlatConverter"; }

  BasicMPFlatConverter() {
    InitOptions();
  }

  std::unique_ptr<ConverterQuery> MakeConverterQuery() {
      return std::unique_ptr<FlatConverterQuery<Impl>>(
          new FlatConverterQuery<Impl>(*(Impl*)this));
  }


  //////////////////////////// CONVERTERS OF STANDRAD MP ITEMS //////////////////////////////
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////
public:
  void Convert(typename Model::MutCommonExpr ) {
    /// Converting on demand, see VisitCommonExpr
  }

  void Convert(typename Model::MutObjective obj) {
    if (NumericExpr e = obj.nonlinear_expr()) {
      LinearExpr &linear = obj.linear_expr();
      const auto eexpr=MP_DISPATCH( Visit(e) );
      linear.AddTerms(eexpr.GetAE());
      if (std::fabs(eexpr.constant_term())!=0.0) {   // TODO use constant (in the extra info)
        linear.AddTerm(MakeFixedVar(eexpr.constant_term()), 1.0);
      }
      if (!eexpr.is_affine())                     // higher-order terms
        obj.set_extra_info(
              typename Model::Params::ExtraItemInfo::ObjExtraInfo{
                0.0, std::move(eexpr.GetQT()) } );
      obj.unset_nonlinear_expr();
    } // Modifying the original objective by replacing the expr
  }

  void Convert(typename Model::MutAlgebraicCon con) {
    if (NumericExpr e = con.nonlinear_expr()) {
      LinearExpr &linear = con.linear_expr();
      const auto ee=MP_DISPATCH( Visit(e) );
      linear.AddTerms(ee.GetAE());
      if (!ee.is_affine())                       // higher-order terms
        con.set_extra_info( std::move(ee.GetQT()) );
      con.set_lb(con.lb() - ee.constant_term());
      con.set_ub(con.ub() - ee.constant_term());
      con.unset_nonlinear_expr();                  // delete the non-linear expr
    } // Modifying the original constraint by replacing the expr
  }

  void Convert(typename Model::MutLogicalCon e) {
    const auto resvar = MP_DISPATCH( Convert2Var(e.expr()) );
    PropagateResultOfInitExpr(resvar, 1.0, 1.0, +Context());
  }

  void PropagateResultOfInitExpr(int var, double lb, double ub, Context ctx) {
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
    if (ee.is_affine())
      AddConstraint(LinearDefiningConstraint(r, std::move(ee.GetAE())));
    else
      AddConstraint(QuadraticDefiningConstraint(r, std::move(ee)));
    return r;
  }
  /// Makes an affine expr representing just one variable
  AffineExpr Convert2VarAsAffineExpr(EExpr&& ee) {
    return AffineExpr::Variable{Convert2Var(std::move(ee))};
  }

  AffineExpr Convert2AffineExpr(EExpr&& ee) {
    if (ee.is_affine())
      return std::move(ee.GetAE());
    return Convert2VarAsAffineExpr(std::move(ee));           // just simple, whole QuadExpr
  }

  PreprocessInfoStd ComputeBoundsAndType(const QuadExpr& ee) {
    auto bntAE = ComputeBoundsAndType(ee.GetAE());
    auto bntQT = ComputeBoundsAndType(ee.GetQT());
    return AddBoundsAndType(bntAE, bntQT);
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

  PreprocessInfoStd ComputeBoundsAndType(const QuadTerms& qt) {
    PreprocessInfoStd result;
    result.lb_ = result.ub_ = 0.0;
    result.type_ = var::INTEGER;
    result.linexp_type_ = var::INTEGER;
    auto& model = MP_DISPATCH( GetModel() );
    for (int i=0; i<qt.num_terms(); ++i) {
      auto coef = qt.coef(i);
      auto v1 = model.var(qt.var1(i));
      auto v2 = model.var(qt.var2(i));
      auto prodBnd = ProductBounds(v1, v2);
      if (coef >= 0.0) {
        result.lb_ += coef * prodBnd.first;
        result.ub_ += coef * prodBnd.second;
      } else {
        result.lb_ += coef * prodBnd.second;
        result.ub_ += coef * prodBnd.first;
      }
      if (var::INTEGER!=v1.type() || var::INTEGER!=v2.type() || !is_integer(coef)) {
        result.type_=var::CONTINUOUS;
        result.linexp_type_=var::CONTINUOUS;
      }
    }
    return result;
  }

  template <class Var>
  std::pair<double, double> ProductBounds(Var x, Var y) const {
    auto lx=x.lb(), ly=y.lb(), ux=x.ub(), uy=y.ub();
    std::array<double, 4> pb{lx*ly, lx*uy, ux*ly, ux*uy};
    return {*std::min_element(pb.begin(), pb.end()), *std::max_element(pb.begin(), pb.end())};
  }

  PreprocessInfoStd AddBoundsAndType(const PreprocessInfoStd& bnt1,
                                     const PreprocessInfoStd& bnt2) {
    return {bnt1.lb()+bnt2.lb(), bnt1.ub()+bnt2.ub(),
      var::INTEGER==bnt1.type() && var::INTEGER==bnt2.type() ?
            var::INTEGER : var::CONTINUOUS};
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
    for (unsigned i=0; i<N; ++i, ++itea)
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
    ee[0].Subtract(std::move(ee[1]));
    return AssignResultToArguments(
          FuncConstraint(                  // comparison with linear expr only
            Convert2AffineExpr( std::move(ee[0]) ) ) );
  }

  template <class ExprArray, size_t N>
  void Exprs2EExprs(const ExprArray& ea, std::array<EExpr, N>& result) {
    assert(ea.size() == result.size());
    auto itea = ea.begin();
    for (size_t i=0; i<N; ++i, ++itea)
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
    if (index >= (int)common_exprs_.size())
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

  EExpr VisitMul(BinaryExpr e) {
    auto el = Convert2EExpr(e.lhs());
    auto er = Convert2EExpr(e.rhs());
    return QuadratizeOrLinearize( el, er );
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

  EExpr VisitAllDiff(PairwiseExpr e) {
    if (expr::ALLDIFF != e.kind())
      throw std::logic_error("NOT_ALLDIFF NOT IMPLEMENTED");
    return VisitFunctionalExpression<AllDiffConstraint>(e);
  }

  EExpr VisitPLTerm(PLTerm e) {
    int num_breakpoints = e.num_breakpoints();
    std::vector<double> slopes(num_breakpoints+1), breakpoints(num_breakpoints);
    for (int i = 0; i < num_breakpoints; ++i) {
      slopes[i] = e.slope(i);
      breakpoints[i] = e.breakpoint(i);
    }
    slopes.back() = e.slope(num_breakpoints);
    return AssignResultToArguments( PLConstraint(
          PLConstraint::Arguments{ Convert2Var(e.arg()) },
          PLConstraint::Parameters{ breakpoints, slopes, 0.0, 0.0 } ) );
  }


  ////////////////////////////////////////////////////
  /////////////// NONLINEAR FUNCTIONS ////////////////
  ////////////////////////////////////////////////////
  EExpr VisitPowConstExp(BinaryExpr e) {
    auto c = Cast<NumericConstant>(e.rhs()).value();
    if (2.0==c) {                            // Quadratic
      auto el = Convert2EExpr(e.lhs());
      return QuadratizeOrLinearize(el, el);
    }
    return AssignResultToArguments( PowConstraint(
      PowConstraint::Arguments{ Convert2Var(e.lhs()) },
      PowConstraint::Parameters{ c } ) );
  }

  EExpr VisitPow2(UnaryExpr e) {     // MIP could have better conversion for pow2
    auto el = Convert2EExpr(e.arg());
    return QuadratizeOrLinearize(el, el);
    /* TODO Can do better for integer variables if we redefine pow:
    return AssignResultToArguments( PowConstraint(
      PowConstraint::Arguments{ Convert2Var(e.arg()) },
      PowConstraint::Parameters{ 2.0 } ) ); */
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


  /// Depending on the target backend
  /// Currently only quadratize higher-order products
  /// Can change arguments. They could point to the same
  /// TODO multiply-out optional
  /// CAUTION: allows &el==&er, needed from Pow2
  EExpr QuadratizeOrLinearize(EExpr& el, EExpr& er) {
    if (!el.is_affine() && !er.is_constant())
      el = Convert2AffineExpr(std::move(el));      // will convert to a new var now
    if (!er.is_affine() && !el.is_constant())
      er = Convert2AffineExpr(std::move(er));
    return MultiplyOut(el, er);
  }

  EExpr MultiplyOut(const EExpr& el, const EExpr& er) {
    assert(el.is_affine() && er.is_affine());
    EExpr result;
    result.constant_term(el.constant_term() * er.constant_term());
    if (0.0!=std::fabs(er.constant_term()))
      for (const auto& term: el.GetAE()) {
        result.AddLinearTerm(term.var_index(), term.coef() * er.constant_term());
      }
    if (0.0!=std::fabs(el.constant_term()))
      for (const auto& term: er.GetAE()) {
        result.AddLinearTerm(term.var_index(), term.coef() * el.constant_term());
      }
    for (const auto& termL: el.GetAE()) {
      for (const auto& termR: er.GetAE()) {
        result.AddQuadraticTerm(termL.var_index(), termR.var_index(),
                                termL.coef() * termR.coef());
      }
    }
    return result;
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
                           ccf.message(),
                           MP_DISPATCH( GetBackend() ).GetBackendName() ) );
          }
        }
      }
    }
  }

  /// fAllSOS2: if false, only groups with sosno<0
  void ConvertSOSCollection(ArrayRef<int> sosno, ArrayRef<double> ref,
                            bool fAllSOS2) {
    assert(sosno.size() == ref.size());
    std::map< int, std::map< double, int > > sos_map;
    for (auto i=ref.size(); i--; )
      if (sosno[i]) {
        auto& sos_group = sos_map[sosno[i]];
        if (sos_group.end() != sos_group.find(ref[i]))
          MP_RAISE(fmt::format(
                     "In SOS group {}, repeated weight {}",
                     sosno[i], ref[i]));
        sos_group[ref[i]] = i;
      }
    for (const auto& group: sos_map) {
      std::vector<int> vars;
      vars.reserve(group.second.size());
      std::vector<double> weights;
      weights.reserve(group.second.size());
      for (const auto& wv: group.second) {
        weights.push_back(wv.first);
        vars.push_back(wv.second);
      }
      if (group.first<0 || fAllSOS2)
        AddConstraint(SOS2Constraint(vars, weights));
      else
        AddConstraint(SOS1Constraint(vars, weights));
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
    const auto argvar = c.GetArguments()[0];
    const auto lb = this->lb(argvar),
        ub = this->ub(argvar);
    if (lb>=0.0) {
      prepro.set_result_var(argvar);
      return;
    } else if (ub<=0.0) {
      prepro.set_result_var(                   // create newvar = -argvar
            AssignResultToArguments(
                  LinearDefiningConstraint({ {-1.0}, {argvar}, 0.0 })).
                              get_representing_variable());
      return;
    }
    prepro.narrow_result_bounds(0.0, std::max(-lb, ub));
    prepro.set_result_type( var_type(argvar) );
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
  /// Currently we should propagate to all arguments, be it always the CTX_MIX.

  /// By default, declare mixed context
  template <class Constraint>
  void PropagateResult(Constraint& con, double lb, double ub, Context ctx) {
    internal::Unused(con, lb, ub, ctx);
    con.SetContext(Context::CTX_MIX);
    for (const auto a: con.GetArguments())
      PropagateResultOfInitExpr(a, this->MinusInfty(), this->Infty(), Context::CTX_MIX);
  }

  void PropagateResult(LinearDefiningConstraint& con, double lb, double ub, Context ctx) {
    internal::Unused(con, lb, ub, ctx);
    con.AddContext(ctx);
    for (const auto& term: con.GetAffineExpr())
      PropagateResultOfInitExpr(term.var_index(), this->MinusInfty(), this->Infty(), Context::CTX_MIX);
  }

  void PropagateResult(QuadraticDefiningConstraint& con, double lb, double ub, Context ctx) {
    internal::Unused(con, lb, ub, ctx);
    con.AddContext(ctx);
    const auto& args = con.GetArguments();
    for (const auto& term: args.GetAE())
      PropagateResultOfInitExpr(term.var_index(), this->MinusInfty(), this->Infty(), Context::CTX_MIX);
    const auto& qt = args.GetQT();
    for (int i=0; i<qt.num_terms(); ++i) {
      PropagateResultOfInitExpr(qt.var1(i), this->MinusInfty(), this->Infty(), Context::CTX_MIX);
      PropagateResultOfInitExpr(qt.var2(i), this->MinusInfty(), this->Infty(), Context::CTX_MIX);
    }
  }

  void PropagateResult(LinearConstraint& con, double lb, double ub, Context ctx) {
    internal::Unused(con, lb, ub, ctx);
    for (const auto& v: con.vars())
      PropagateResultOfInitExpr(v, this->MinusInfty(), this->Infty(), Context::CTX_MIX);
  }

  void PropagateResult(IndicatorConstraintLinLE& con, double lb, double ub, Context ctx) {
    internal::Unused(con, lb, ub, ctx);
    PropagateResultOfInitExpr(con.get_binary_var(),
                              this->MinusInfty(), this->Infty(), Context::CTX_MIX);
    for (const auto& v: con.get_lin_vars())
      PropagateResultOfInitExpr(v, this->MinusInfty(), this->Infty(), Context::CTX_MIX);
  }

  template <int type>
  void PropagateResult(SOS_1or2_Constraint<type>& con, double lb, double ub, Context ctx) {
    internal::Unused(con, lb, ub, ctx);
    for (const auto& v: con.get_vars())
      PropagateResultOfInitExpr(v, this->MinusInfty(), this->Infty(), Context::CTX_MIX);
  }

  void PropagateResult(NotConstraint& con, double lb, double ub, Context ctx) {
    internal::Unused(con, lb, ub, ctx);
    con.AddContext(ctx);
    for (const auto a: con.GetArguments())
      PropagateResultOfInitExpr(a, 1.0-ub, 1.0-lb, -ctx);
  }

  void PropagateResult(ConjunctionConstraint& con, double lb, double ub, Context ctx) {
    internal::Unused(con, lb, ub, ctx);
    con.AddContext(ctx);
    for (const auto a: con.GetArguments())
      PropagateResultOfInitExpr(a, lb, 1.0, +ctx);
  }

  void PropagateResult(DisjunctionConstraint& con, double lb, double ub, Context ctx) {
    internal::Unused(con, lb, ub, ctx);
    con.AddContext(ctx);
    for (const auto a: con.GetArguments())
      PropagateResultOfInitExpr(a, 0.0, ub, +ctx);
  }

  void PropagateResult(IfThenConstraint& con, double lb, double ub, Context ctx) {
    internal::Unused(con, lb, ub, ctx);
    con.AddContext(ctx);
    auto& args = con.GetArguments();
    PropagateResultOfInitExpr(args[0], 0.0, 1.0, Context::CTX_MIX);
    PropagateResultOfInitExpr(args[1], this->MinusInfty(), this->Infty(), +ctx);
    PropagateResultOfInitExpr(args[2], this->MinusInfty(), this->Infty(), -ctx);
  }

  void PropagateResult(AllDiffConstraint& con, double lb, double ub, Context ctx) {
    internal::Unused(con, lb, ub, ctx);
    con.AddContext(ctx);
    // TODO go into arguments
  }

  void PropagateResult(LE0Constraint& con, double lb, double ub, Context ctx) {
    internal::Unused(con, lb, ub, ctx);
    con.AddContext(ctx);
  }

  void PropagateResult(EQ0Constraint& con, double lb, double ub, Context ctx) {
    internal::Unused(con, lb, ub, ctx);
    con.AddContext(ctx);
  }


  //////////////////////////// CUSTOM CONSTRAINTS CONVERSION ////////////////////////////
  ///
  //////////////////////////// SPECIFIC CONSTRAINT CONVERTERS ///////////////////////////

  USE_BASE_CONSTRAINT_CONVERTERS(BasicConstraintConverter)      // reuse default converters

  /// Assume mixed context if not set in the constraint
  /// TODO Make sure context is always propagated for all constraints and objectives
  template <class Constraint>
  void RunConversion(const Constraint& con) {
    if (con.HasContext())
      if (con.GetContext().IsNone())
        con.SetContext(Context::CTX_MIX);
    MP_DISPATCH(Convert(con););
  }


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
  template <class VarArray>
  double lb_array(const VarArray& va) const { return this->GetModel().lb_array(va); }
  template <class VarArray>
  double ub_array(const VarArray& va) const { return this->GetModel().ub_array(va); }
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
    return int(var_info_.size())>var &&
        nullptr!=var_info_[var].pInitExpr;
  }

  BasicConstraintKeeper* GetInitExpression(int var) {
    assert(HasInitExpression(var));
    return var_info_[var].pInitExpr;
  }


  ///////////////////////////////////////////////////////////////////////
  //////////////////// SOLUTION REPORTING FROM BACKEND //////////////////
  ///////////////////////////////////////////////////////////////////////
public:
  void HandleSolution(int status, fmt::CStringRef msg,
      const double *x, const double *y, double obj) {
    MP_DISPATCH( GetSolH() ).HandleSolution(status, msg, x, y, obj);
  }

  void HandleFeasibleSolution(fmt::CStringRef msg,
      const double *x, const double *y, double obj) {
    MP_DISPATCH( GetSolH() ).HandleFeasibleSolution(msg, x, y, obj);
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
    this->AddOption("cvt:pre:all",
        "0/1*: Set to 0 to disable all presolve in the converter.",
        options_.preprocessAnything_);
    this->AddOption("cvt:pre:eqresult",
        "0/1*: Preprocess reified equality comparison's boolean result bounds.",
        options_.preprocessEqualityResultBounds_);
    this->AddOption("cvt:pre:eqbinary",
        "0/1*: Preprocess reified equality comparison with a binary variable.",
        options_.preprocessEqualityBvar_);
  }

protected:
  bool CanPreprocess(int f) const {
    return 0!=options_.preprocessAnything_ && 0!=f;
  }

};

} // namespace mp

#endif // CONVERTER_FLAT_H
