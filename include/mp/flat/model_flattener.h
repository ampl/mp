#ifndef MODEL_FLATTENER_H
#define MODEL_FLATTENER_H

#include <utility>
#include <unordered_map>
#include <map>
#include <cmath>

#include "mp/problem.h"    // for ToLinTerms()
#include "mp/converter-base.h"
#include "mp/expr-visitor.h"
#include "mp/flat/eexpr.h"
#include "mp/flat/constraints_std.h"
#include "mp/flat/obj_std.h"
#include "mp/presolve.h"


namespace mp {

/// Convert mp::LinearExpr to LinTerms
LinTerms ToLinTerms(const LinearExpr& e) {
  LinTerms le;
  le.reserve(e.num_terms());
  for (auto it=e.begin(); it!=e.end(); ++it) {
    le.add_term(it->coef(), it->var_index());
  }
  return le;
}


/// ModelFlattener: it walks and "flattens" most expressions
/// by replacing them by a result variable and constraints.
/// Such replacement is performed by a FlatConverter object.
/// Such constraints might need to be converted to others, which is
/// handled by overloaded methods in descendants of FlatConverter.
/// @param Impl: final CRTP class
/// @param Model: the model class representing the input instance
/// @param FlatConverter: the FlatConverter type
template <class Impl, class Model, class FlatConverter>
class ModelFlattener :
    public ExprConverter<Impl, EExpr>,
    public BasicConverter<Model>
{
public:
  using ModelType = Model;
  using FlatConverterType = FlatConverter;

public:
  using Var = typename FlatConverter::Var;
  using EExprType = EExpr;
  using VarArray = std::vector<Var>;

protected:
  using ClassName = ModelFlattener;
  using BaseExprVisitor = ExprVisitor<Impl, EExpr>;
  using BaseConverter = BasicConverter<Model>;

  using EExprArray = std::vector<EExpr>;

  using BaseConverter::GetEnv;

public:
  static const char* GetName() { return "ExprFlattener"; }

  ModelFlattener(Env& e) : BaseConverter(e), flat_cvt_(e) { }

public:
  /// INCREMENTAL INTERFACE
  /// These guys used from outside to feed a model to be converted
  /// and forwarded to a backend
  /// Currently only used for testing
  void InputVariables(int n, const double* lb, const double* ub, const var::Type* ty) {
    GetModel().AddVars(n, lb, ub, ty);
  }
  /// Add vector of variables. Type: var::CONTINUOUS by default
  std::vector<int> AddVars(std::size_t nvars,
                           double lb=-INFINITY, double ub=INFINITY,
                           var::Type type = var::CONTINUOUS) {
    return GetModel().AddVars(nvars, lb, ub, type);
  }
  void InputObjective(obj::Type t,
                      int nnz, const double* c, const int* v, NumericExpr e=NumericExpr()) {
    typename Model::LinearObjBuilder lob = GetModel().AddObj(t, e);
    for (int i=0; i!=nnz; ++i) {
      lob.AddTerm(v[i], c[i]);
    }
  }
  void InputAlgebraicCon(int nnz, const double* c, const int* v,
                         double lb, double ub, NumericExpr e=NumericExpr()) {
    typename Model::MutAlgebraicCon mac = GetModel().AddCon(lb, ub);
    typename Model::LinearConBuilder lcb = mac.set_linear_expr(nnz);
    for (int i=0; i!=nnz; ++i)
      lcb.AddTerm(v[i], c[i]);
    mac.set_nonlinear_expr(e);
  }


  //////////////////////////// CONVERTERS OF STANDRAD MP ITEMS //////////////////////////////
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////
public:
  /// Convert the whole model, e.g., after reading from NL
  void ConvertModel() {
    GetFlatCvt().StartModelInput();
    MP_DISPATCH( ConvertStandardItems() );
    GetFlatCvt().FinishModelInput();      // Chance to flush to the Backend
  }

protected:
  void ConvertStandardItems() {
    ConvertVars();
    int num_common_exprs = GetModel().num_common_exprs();
    for (int i = 0; i < num_common_exprs; ++i)
      MP_DISPATCH( Convert( GetModel().common_expr(i) ) );
    if (int num_objs = GetModel().num_objs())
      for (int i = 0; i < num_objs; ++i)
        MP_DISPATCH( Convert( GetModel().obj(i) ) );
    if (int n_cons = GetModel().num_algebraic_cons())
      for (int i = 0; i < n_cons; ++i)
        MP_DISPATCH( ConvertAlgCon( i ) );
    if (int n_lcons = GetModel().num_logical_cons())
      for (int i = 0; i < n_lcons; ++i)
        MP_DISPATCH( ConvertLogicalCon( i ) );

    ////////////////////////////////////////////////
    MP_DISPATCH( ConvertSOSConstraints() );
  }

  /// TODO this can be slow
  void ConvertVars() {
    std::vector<double> lbs(GetModel().num_vars());
    std::vector<double> ubs(GetModel().num_vars());
    std::vector<var::Type> types(GetModel().num_vars());
    for (int i=GetModel().num_vars(); i--; ) {
      const auto mpvar = GetModel().var(i);
      lbs[i] = mpvar.lb();
      ubs[i] = mpvar.ub();
      types[i] = mpvar.type();
    }
    auto vnr = GetFlatCvt().AddVars(lbs, ubs, types);
    GetCopyBridge().AddEntry({
          GetPresolver().GetSourceNodes().GetVarValues().MakeSingleKey().
                               Add(lbs.size()),
          vnr });
  }

  void Convert(typename ModelType::MutCommonExpr ) {
    /// Converting on demand, see VisitCommonExpr
  }

  void Convert(typename ModelType::MutObjective obj) {
      auto le = ToLinTerms(obj.linear_expr());
      NumericExpr e = obj.nonlinear_expr();
      EExpr eexpr;
      if (e) {
        eexpr=MP_DISPATCH( Visit(e) );
        le.add_lin_exp(eexpr.GetAE());
        if (std::fabs(eexpr.constant_term())!=0.0) {
          /// TODO use constant (in the extra info)
          le.add_term(1.0, MakeFixedVar(eexpr.constant_term()));
        }
      }
      LinearObjective lo { obj.type(),
            std::move(le.coefs()), std::move(le.vars()) };
      GetFlatCvt().AddObjective( // TODO save & convert different types
            QuadraticObjective{std::move(lo),
                               std::move(eexpr.GetQT())});
  }

  void ConvertAlgCon(int i) {
    auto con = GetModel().algebraic_con(i);
    auto le = ToLinTerms(con.linear_expr());
    EExpr ee;
    if (NumericExpr e = con.nonlinear_expr()) {
      ee=MP_DISPATCH( Visit(e) );
      le.add_lin_exp(ee.GetAE());
    }
    auto lc = RangeLinCon{
    { std::move(le.coefs()), std::move(le.vars()) },
    { con.lb() - ee.constant_term(), con.ub() - ee.constant_term() } };
    lc.sort_terms();
    pre::NodeRange nr;
    if (ee.is_affine())
      nr = AddConstraint( std::move(lc) );
    else                                   // higher-order terms
      nr = AddConstraint(
          QuadraticConstraint{std::move(lc),
                              std::move(ee.GetQT())} );
    GetCopyBridge().AddEntry( {
            GetPresolver().GetSourceNodes().GetConValues()(0).Add(),
            nr
          } );
  }

  void ConvertLogicalCon(int i) {
    auto e = GetModel().logical_con(i);
    const auto resvar = MP_DISPATCH( Convert2Var(e.expr()) );
    GetFlatCvt().FixAsTrue(resvar);            // TODO avoid creating the variable
    assert(GetFlatCvt().HasInitExpression(resvar));
    const auto& ie = GetFlatCvt().GetInitExpression(resvar);
    /// TODO check that logical cons' values come after algebraic ones
    GetCopyBridge().AddEntry( {
            GetPresolver().GetSourceNodes().GetConValues()(0).Add(),
            ie.GetCK()->GetValueNode().Select( ie.GetIndex() )
          } );
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
    return GetFlatCvt().Convert2Var(std::move(ee));
  }
  /// Makes an affine expr representing just one variable
  AffExp Convert2VarAsAffineExpr(EExpr&& ee) {
    return AffExp::Variable{Convert2Var(std::move(ee))};
  }
  AffExp Convert2AffineExpr(EExpr&& ee) {
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
    auto vc = GetFlatCvt().AssignResult2Args(
          std::forward<FuncConstraint>(fc));
    return vc.is_var() ? EExpr{ EExpr::Variable{ vc.get_var() } } :
                         EExpr{ EExpr::Constant{ vc.get_const() } };
  }

  /// Generic relational expression visitor
  /// Can produce a new variable/expression and specified constraints on it
  template <class FuncConstraint, class ExprArray=std::initializer_list<Expr> >
  EExpr VisitRelationalExpression(ExprArray ea) {
    std::array<EExpr, 2> ee;
    Exprs2EExprs(ea, ee);
    ee[0].subtract(std::move(ee[1]));
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
      EExpr eexpr(ToLinTerms(ce.linear_expr()));
      if (ce.nonlinear_expr())
        eexpr.add( Convert2EExpr(ce.nonlinear_expr()) );
      common_exprs_[index] = Convert2Var(std::move(eexpr));
    }
    return EExpr::Variable{ common_exprs_[index] };
  }

  EExpr VisitMinus(UnaryExpr e) {
    auto ee = Convert2EExpr(e.arg());
    ee.negate();
    return ee;
  }

  EExpr VisitAdd(BinaryExpr e) {
    auto ee = Convert2EExpr(e.lhs());
    ee.add( Convert2EExpr(e.rhs()) );
    return ee;
  }

  EExpr VisitSub(BinaryExpr e) {
    auto el = Convert2EExpr(e.lhs());
    auto er = Convert2EExpr(e.rhs());
    er.negate();
    el.add(er);
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
      sum.add( MP_DISPATCH( Convert2EExpr(*i) ) );
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
  EExpr VisitForAll(IteratedLogicalExpr e) {
    return VisitFunctionalExpression<ConjunctionConstraint>(e);
  }

  EExpr VisitOr(BinaryLogicalExpr e) {
    return VisitFunctionalExpression<DisjunctionConstraint>({ e.lhs(), e.rhs() });
  }
  EExpr VisitExists(IteratedLogicalExpr e) {
    return VisitFunctionalExpression<DisjunctionConstraint>(e);
  }

  EExpr VisitIf(IfExpr e) {
    return VisitFunctionalExpression<IfThenConstraint>({
                e.condition(), e.then_expr(), e.else_expr() });
  }

  EExpr VisitIff(BinaryLogicalExpr e) {
    return VisitRelationalExpression<EQ0Constraint>({ e.lhs(), e.rhs() });
  }

  EExpr VisitAllDiff(PairwiseExpr e) {
    if (expr::ALLDIFF != e.kind())
      throw std::logic_error("NOT_ALLDIFF NOT IMPLEMENTED");
    return VisitFunctionalExpression<AllDiffConstraint>(e);
  }

  /// Numberof: const, var
  EExpr VisitNumberOf(typename BaseExprVisitor::NumberOfExpr e) {
    VarArray va;
    va.reserve(e.num_args());
    EExpr valexpr = Convert2EExpr(e.arg(0));
    if (valexpr.is_constant()) {
      for (int i=1; i<e.num_args(); ++i)
        va.push_back(Convert2Var(e.arg(i)));
      return AssignResult2Args( NumberofConstConstraint{
                                  va, { valexpr.constant_term() }
                                } );
    } else {
      va.push_back(Convert2Var(std::move(valexpr)));
      for (int i=1; i<e.num_args(); ++i)
        va.push_back(Convert2Var(e.arg(i)));
      return AssignResult2Args( NumberofVarConstraint{ va } );
    }
  }

  EExpr VisitCount(CountExpr ce) {
    return VisitFunctionalExpression<CountConstraint>(ce);
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

  /// TODO how to bridge them if they are no real items in NL?
  void ConvertSOSConstraints() {
    if (sos()) {
      auto sosno = GetModel().
          ReadIntSuffix( { "sosno", suf::Kind::VAR } );
      auto ref = GetModel().
          ReadDblSuffix( { "ref", suf::Kind::VAR } );
      if (sosno && ref)
        MP_DISPATCH( ConvertSOSCollection(sosno, ref, false) );
    }
    if (sos2_ampl_pl()) {
      auto sos = GetModel().
          ReadIntSuffix( { "sos", suf::Kind::VAR } );
      auto sosref = GetModel().
          ReadDblSuffix( { "sosref", suf::Kind::VAR } );
      if (sos && sosref)
        MP_DISPATCH( ConvertSOSCollection(sos, sosref, true) );
    }
  }

  /// fAllSOS2: if false, only groups with sosno<0 are treated as SOS2
  void ConvertSOSCollection(ArrayRef<int> sosno, ArrayRef<double> ref,
                            bool fAllSOS2) {
    assert(sosno.size() == ref.size());
    std::map< int, std::multimap< double, int > > sos_map;
    for (auto i=ref.size(); i--; )
      if (sosno[i]) {
        auto& sos_group = sos_map[sosno[i]];
        if (sos_group.end() != sos_group.find(ref[i])) {
          GetFlatCvt().AddWarning( "SOS_repeated_weight",
                                   "An SOS constraint has repeated weights, "
                                   "solver might reject it" );
        }
        sos_group.insert( {ref[i], i} );
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

  /// Depending on the target backend
  /// Currently only quadratize higher-order products
  /// Can change arguments. They could point to the same
  /// TODO multiply-out optional
  /// TODO Move to FlatCvt?
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
      result.GetAE().add_lin_exp(el.GetAE());  // no const here
      result.GetAE() *= er.constant_term();
      result.GetQT().add(el.GetQT());
      result.GetQT() *= er.constant_term();
    }
    if (0.0!=std::fabs(el.constant_term())) {
      {
        AffExp ae2 = er.GetAE();
        ae2 *= el.constant_term();
        result.GetAE().add_aff_exp(ae2);
      }
      result.GetQT().add(er.GetQT());
      result.GetQT() *= el.constant_term();
    }
    const auto& ae1 = el.GetAE();
    const auto& ae2 = er.GetAE();
    for (auto i1 = ae1.size(); i1--; ) {
      for (auto i2 = ae2.size(); i2--; ) {
        result.add_qp_term(ae1.coef(i1) * ae2.coef(i2),
                           ae1.var(i1), ae2.var(i2) );
      }
    }
    result.sort_terms();      // eliminate 0's and duplicates
    return result;
  }


protected:
  //////////////////////// ADD CUSTOM CONSTRAINT ///////////////////////
  //////////////////////// Takes ownership /////////////////////////////
  template <class Constraint>
  pre::NodeRange AddConstraint(Constraint&& con) {
    return GetFlatCvt().AddConstraint(std::move(con));
  }

  //////////////////////////// UTILITIES /////////////////////////////////
  ///
protected:
  pre::Presolver& GetPresolver() { return GetFlatCvt().GetPresolver(); }

private:
  std::unordered_map<double, int> map_fixed_vars_;

  std::vector<int> common_exprs_;               // variables equal to the result

protected:

  //////////////////////////// CREATE OR FIND A FIXED VARIABLE //////////////////////////////
  int MakeFixedVar(double value) // TODO use proper const term in obj
  { return GetFlatCvt().MakeFixedVar(value); }

  /// Presolve bridge copying values between model items
  pre::CopyBridge& GetCopyBridge() { return GetFlatCvt().GetCopyBridge(); }

  ///////////////////////////////////////////////////////////////////////
  /////////////////////// OPTIONS /////////////////////////
  ///
private:
  struct Options {
    int sos_ = 1;
    int sos2_ = 1;
  };
  Options options_;

protected:
  int sos() const { return options_.sos_; }
  int sos2_ampl_pl() const { return options_.sos2_; }

public:
  void InitOptions() {
    InitOwnOptions( );
    GetFlatCvt().InitOptions( );
  }

private:
  void InitOwnOptions() {
    GetEnv().AddOption("cvt:sos sos",
        "0/1*: Whether to honor declared suffixes .sosno and .ref describing "
        "SOS sets. Each distinct nonzero .sosno "
        "value designates an SOS set, of type 1 for "
        "positive .sosno values and of type 2 for "
        "negative values.  The .ref suffix contains "
        "corresponding reference values used to order the variables.",
        options_.sos_, 0, 1);
    GetEnv().AddOption("cvt:sos2 sos2",
        "0/1*: Whether to honor SOS2 constraints for nonconvex "
        "piecewise-linear terms, using suffixes .sos and .sosref "
        "provided by AMPL.",
        options_.sos2_, 0, 1);
  }


public:
  /// The model as input from NL
  const ModelType& GetModel() const { return model_; }
  /// The model as input from NL
  ModelType& GetModel() { return model_; }
  /// The model as input from NL
  const ModelType& GetInputModel() const { return GetModel(); }
  /// The model as input from NL
  ModelType& GetInputModel() { return GetModel(); }

  const FlatConverter& GetFlatCvt() const { return flat_cvt_; }
  FlatConverter& GetFlatCvt() { return flat_cvt_; }


private:
  ModelType model_;
  FlatConverter flat_cvt_;

};

/// A 'final' ModelFlattener in a hierarchy
template <template <typename, typename, typename> class ModelFlt,
          class Model, class FlatCvt>
class ModelFltImpl :
    public ModelFlt<
      ModelFltImpl<ModelFlt, Model, FlatCvt>, Model, FlatCvt> {
  using Base = ModelFlt<
    ModelFltImpl<ModelFlt, Model, FlatCvt>, Model, FlatCvt>;
public:
  ModelFltImpl(Env& e) : Base(e) { }
};


} // namespace mp

#endif // MODEL_FLATTENER_H
