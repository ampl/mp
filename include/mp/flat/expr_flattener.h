#ifndef EXPR_FLATTENER_H
#define EXPR_FLATTENER_H

#include <utility>
#include <unordered_map>
#include <map>
#include <cmath>

#include "mp/flat/preprocess.h"
#include "mp/flat/basic_converters.h"
#include "mp/flat/MIP/mp2mip.h"
#include "mp/expr-visitor.h"
#include "mp/flat/eexpr.h"
#include "mp/flat/convert_functional.h"
#include "mp/flat/model.h"
#include "mp/flat/std_constr.h"

namespace mp {

/// ExprFlattener: it walks and "flattens" most expressions
/// by replacing them by a result variable and constraints.
/// Such replacement is performed by a FlatConverter object.
/// Such constraints might need to be converted to others, which is
/// handled by overloaded methods in FlatConverter.
template <class Impl, class FlatConverter>
class ExprFlattener
    : public BasicMPConverter<Impl,  // TODO become owned by NLSolver
        typename FlatConverter::BackendType,
        typename FlatConverter::ModelType>,
      public ExprVisitor<Impl, EExpr>
{
  FlatConverter flat_cvt_;
public:
  const FlatConverter& FlatCvt() const { return flat_cvt_; }
  FlatConverter& FlatCvt() { return flat_cvt_; }
  using BackendType = typename FlatConverter::BackendType;
  using Model = typename FlatConverter::ModelType;

public:
  using Var = typename FlatConverter::Var;
  using EExprType = EExpr;
  using VarArray = std::vector<Var>;

protected:
  using ClassName = ExprFlattener;
  using BaseConverter = BasicMPConverter<Impl,  // TODO become owned by
    typename FlatConverter::BackendType,
    typename FlatConverter::ModelType>;
  using BaseExprVisitor = ExprVisitor<Impl, EExpr>;

  using EExprArray = std::vector<EExpr>;


public:
  static const char* GetConverterName() { return "BasicMPFlatConverter"; }

  ExprFlattener() { }

  std::unique_ptr<ConverterQuery> MakeConverterQuery() {
      return std::unique_ptr<FlatConverterQuery<Impl>>(
          new FlatConverterQuery<Impl>(*(Impl*)this));
  }


  //////////////////////////// CONVERTERS OF STANDRAD MP ITEMS //////////////////////////////
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////
  /// TOD might create new obj/constraints in the FlatConverter
  /// (in its internal model - to be added),
  /// instead of modifying Problem
public:
  void Convert(typename Model::MutCommonExpr ) {
    /// Converting on demand, see VisitCommonExpr
  }

  void Convert(typename Model::MutObjective obj) {
    if (NumericExpr e = obj.nonlinear_expr()) {
      LinearExpr &linear = obj.linear_expr();
      const auto eexpr=MP_DISPATCH( Visit(e) );
      linear.AddTerms(eexpr.GetAE());
      if (std::fabs(eexpr.constant_term())!=0.0) {
        // TODO use constant (in the extra info)
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
    FlatCvt().FixAsTrue(resvar);  // TODO avoid creting the variable
  }

  void ConvertExtraItems() { FlatCvt().ConvertExtraItems(); }

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
    return FlatCvt().Convert2Var(std::move(ee));
  }
  /// Makes an affine expr representing just one variable
  AffineExpr Convert2VarAsAffineExpr(EExpr&& ee) {
    return AffineExpr::Variable{Convert2Var(std::move(ee))};
  }
  AffineExpr Convert2AffineExpr(EExpr&& ee) {
    if (ee.is_affine())
      return std::move(ee.GetAE());
    return Convert2VarAsAffineExpr(std::move(ee)); // just simple, whole QuadExpr
  }

  /// Generic functional expression array visitor
  /// Can produce a new variable/expression and specified constraints on it
  template <class FuncConstraint, class ExprArray=std::initializer_list<Expr> >
  EExpr VisitFunctionalExpression(ExprArray ea) {
    FuncConstraint fc;
    Exprs2Vars(ea, fc.GetArguments());
    return AssignResult2Args( std::move(fc) );
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
  EExpr AssignResult2Args(FuncConstraint&& fc) {
    auto vc = FlatCvt().AssignResult2Args(
          std::forward<FuncConstraint>(fc));
    return vc.is_var() ? EExpr{EExpr::Variable{ vc.get_var() }} :
                         EExpr{EExpr::Constant{ vc.get_const() }};
  }

  /// Generic relational expression visitor
  /// Can produce a new variable/expression and specified constraints on it
  template <class FuncConstraint, class ExprArray=std::initializer_list<Expr> >
  EExpr VisitRelationalExpression(ExprArray ea) {
    std::array<EExpr, 2> ee;
    Exprs2EExprs(ea, ee);
    ee[0].Subtract(std::move(ee[1]));
    return AssignResult2Args(
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
    for (auto i =
         expr.begin(), end = expr.end(); i != end; ++i)
      sum.Add( MP_DISPATCH( Convert2EExpr(*i) ) );
    return sum;
  }

  EExpr VisitMax(typename BaseExprVisitor::VarArgExpr e) {
    // Why need BaseExprVisitor:: here in g++ 9.2.1?
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
    return AssignResult2Args( PLConstraint(
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
    return AssignResult2Args( PowConstraint(
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

  EExpr VisitPow(BinaryExpr e) {
    auto el = Convert2EExpr(e.lhs());
    auto er = Convert2EExpr(e.rhs());
    if (er.is_constant() && 2.0==er.constant_term())
      return QuadratizeOrLinearize(el, el);
    else
      MP_RAISE("Unsupported: general ^");
    /* TODO Can do better for integer variables if we redefine pow:
    return AssignResultToArguments( PowConstraint(
      PowConstraint::Arguments{ Convert2Var(e.arg()) },
      PowConstraint::Parameters{ 2.0 } ) ); */
  }

  EExpr VisitSqrt(UnaryExpr e) {
    return AssignResult2Args( PowConstraint(
      PowConstraint::Arguments{ Convert2Var(e.arg()) },
      PowConstraint::Parameters{ 0.5 } ) );
  }

  EExpr VisitExp(UnaryExpr e) {
    return VisitFunctionalExpression<ExpConstraint>({ e.arg() });
  }

  EExpr VisitPowConstBase(BinaryExpr e) {
    return AssignResult2Args( ExpAConstraint(
      ExpAConstraint::Arguments{ Convert2Var(e.rhs()) },
      ExpAConstraint::Parameters{ Cast<NumericConstant>(e.lhs()).value() } ) );
  }

  EExpr VisitLog(UnaryExpr e) {
    return VisitFunctionalExpression<LogConstraint>({ e.arg() });
  }

  EExpr VisitLog10(UnaryExpr e) {
    return AssignResult2Args( LogAConstraint(
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
    assert((el.is_affine() && er.is_affine()) ||
           (el.is_constant() || er.is_constant()));
    EExpr result;
    if (0.0!=std::fabs(er.constant_term())) {
      result.GetAE().Add(el.GetAE());
      result.GetAE() *= er.constant_term();
      result.GetQT().AddTerms(el.GetQT());
      result.GetQT() *= er.constant_term();
    }
    if (0.0!=std::fabs(el.constant_term())) {
      {
        AffineExpr ae2 = er.GetAE();
        ae2 *= el.constant_term();
        result.GetAE().AddTerms(ae2);
      }
      result.GetQT().Add(er.GetQT());
      result.GetQT() *= el.constant_term();
    }
    for (const auto& termL: el.GetAE()) {
      for (const auto& termR: er.GetAE()) {
        result.AddQuadraticTerm(termL.var_index(), termR.var_index(),
                                termL.coef() * termR.coef());
      }
    }
    result.SortTerms();
    return result;
  }

public:

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


public:
  //////////////////////// ADD CUSTOM CONSTRAINT ///////////////////////
  //////////////////////// Takes ownership /////////////////////////////
  template <class Constraint>
  void AddConstraint(Constraint&& con) {
    FlatCvt().AddConstraint(std::move(con));
  }

  //////////////////////////// UTILITIES /////////////////////////////////
  ///

private:
  std::unordered_map<double, int> map_fixed_vars_;

  std::vector<int> common_exprs_;               // variables equal to the result

protected:

  //////////////////////////// CREATE OR FIND A FIXED VARIABLE //////////////////////////////
  int MakeFixedVar(double value) // TODO use proper const term in obj
  { return FlatCvt().MakeFixedVar(value); }


  ///////////////////////////////////////////////////////////////////////
  //////////////////// SOLUTION REPORTING FROM BACKEND //////////////////
  ///////////////////////////////////////////////////////////////////////
public:
  /// TODO use pre/postsolve
  void HandleSolution(int status, fmt::CStringRef msg,
      const double *x, const double *y, double obj) {
    if ( MPD( HaveSolH() ) )
      MP_DISPATCH( GetSolH() ).HandleSolution(status, msg, x, y, obj);
    else
      MP_RAISE_WITH_CODE(0, msg);
  }

  void HandleFeasibleSolution(fmt::CStringRef msg,
      const double *x, const double *y, double obj) {
    MP_DISPATCH( GetSolH() ).HandleFeasibleSolution(msg, x, y, obj);
  }

  ///////////////////////////////////////////////////////////////////////
  /////////////////////// OPTIONS /////////////////////////
  ///
public:
  void InitOptions() {
    BaseConverter::InitOptions();
    FlatCvt().InitOptions( MPD(GetBackend()) );
  }
};

/// A 'final' ExprFlattener in a hierarchy
template <template <typename, typename> class ExprFlattener, class FlatCvt>
class ExprFlattenerImpl :
    public ExprFlattener<ExprFlattenerImpl<ExprFlattener, FlatCvt>, FlatCvt> { };

/// TODO Isolate NLSolver from ExprFlattener
template <class Backend,
          template <typename, typename, typename> class Converter,
          class Model = BasicModel< > >
using NLSolverWithFlatBackend = ExprFlattenerImpl<ExprFlattener,
                                  Interface<Converter, Backend, Model> >;


} // namespace mp

#endif // EXPR_FLATTENER_H
