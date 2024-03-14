#ifndef PROBLEM_FLATTENER_H
#define PROBLEM_FLATTENER_H

#include <utility>
#include <unordered_map>
#include <map>
#include <cmath>

#include "mp/problem.h"    // for ToLinTerms()
#include "mp/converter-base.h"
#include "mp/expr-visitor.h"
#include "mp/flat/eexpr.h"
#include "mp/flat/constr_std.h"
#include "mp/flat/obj_std.h"
#include "mp/valcvt.h"


namespace mp {

/// Convert mp::LinearExpr to LinTerms
inline
LinTerms ToLinTerms(const LinearExpr& e) {
  LinTerms le;
  le.reserve(e.num_terms());
  for (auto it=e.begin(); it!=e.end(); ++it) {
    le.add_term(it->coef(), it->var_index());
  }
  return le;
}

/// Write algebraic expression (linear + non-linear.)
template <typename ExprTypes,
          typename LinearExpr, typename NumericExpr,
          typename Namer>
void WriteExpr(fmt::Writer &w, const LinearExpr &linear,
               NumericExpr nonlinear, Namer);

/// Write logical expression
template <typename ExprTypes, typename LogicalExpr, class VN>
void WriteExpr(fmt::Writer &w, LogicalExpr expr, VN);

/// Write algebraic constraint.
template <class ExprTypes, class AlgCon, class VN>
void WriteAlgCon(fmt::Writer &w, const AlgCon &con, VN);


/// ProblemFlattener: it walks and "flattens" most expressions
/// by replacing them by a result variable and constraints.
/// Such replacement is performed by a FlatConverter object.
/// Such constraints might need to be converted to others, which is
/// handled by overloaded methods in descendants of FlatConverter.
/// @param Impl: final CRTP class
/// @param Problem: the model class representing the input instance.
/// Should implement the `mp::BasicProblem` interface.
/// @param FlatConverter: the FlatConverter type
template <class Impl, class Problem, class FlatConverter>
class ProblemFlattener :
    public ExprConverter<Impl, EExpr>,
    public BasicConverter<Problem>
{
public:
  using ProblemType = Problem;
  using FlatConverterType = FlatConverter;

public:
  using Var = typename FlatConverter::Var;
  using EExprType = EExpr;
  using VarArray = std::vector<Var>;

protected:
  using ClassName = ProblemFlattener;
  using BaseExprVisitor = ExprVisitor<Impl, EExpr>;
  using BaseConverter = BasicConverter<Problem>;

  using EExprArray = std::vector<EExpr>;

  using BaseConverter::GetEnv;

public:
  static const char* GetTypeName() { return "ProblemFlattener"; }

  ProblemFlattener(Env& e) : BaseConverter(e), flat_cvt_(e) { }

public:
  /// INCREMENTAL INTERFACE.
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
    typename Problem::LinearObjBuilder lob = GetModel().AddObj(t, e);
    for (int i=0; i!=nnz; ++i) {
      lob.AddTerm(v[i], c[i]);
    }
  }
  void InputAlgebraicCon(int nnz, const double* c, const int* v,
                         double lb, double ub, NumericExpr e=NumericExpr()) {
    typename Problem::MutAlgebraicCon mac = GetModel().AddCon(lb, ub);
    typename Problem::LinearConBuilder lcb = mac.set_linear_expr(nnz);
    for (int i=0; i!=nnz; ++i)
      lcb.AddTerm(v[i], c[i]);
    mac.set_nonlinear_expr(e);
  }


  //////////////////////////// CONVERTERS OF STANDARD MP ITEMS //////////////////////////////
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////
public:
  /// Convert the whole model, e.g., after reading from NL
  void ConvertModel() override {
    GetFlatCvt().StartModelInput();
    MP_DISPATCH( ConvertStandardItems() );
    GetFlatCvt().FinishModelInput();      // Chance to flush to the Backend
  }

  /// Fill model traits
  void FillModelTraits(AMPLS_ModelTraits& mt) override {
    GetFlatCvt().FillModelTraits(mt);
  }

  /// Solver-facing model info:
  /// has unfixed int vars?
  bool HasUnfixedIntVars() const override
  { return GetFlatCvt().GetModelInfo()->NumUnfixedIntVars(); }


protected:
  /// Convert problem items
  void ConvertStandardItems() {
    ////////////////////////// Variables
    ConvertVars();

    ////////////////////// SOS constraints //////////////////////////
    /// SOS2 come before algebraic constraints, so we can omit
    /// AMPL's linearization of SOS2 (in case we use them).
    MP_DISPATCH( ConvertSOSConstraints() );

    ////////////////////////// Common exprs
    int num_common_exprs = GetModel().num_common_exprs();
    for (int i = 0; i < num_common_exprs; ++i) {
      MPD( ExportCommonExpr(i) );
      MP_DISPATCH( Convert( GetModel().common_expr(i) ) );
    }

    ////////////////////////// Objectives
    ifFltCon_ = 0;
    if (int num_objs = GetModel().num_objs())
      for (int i = 0; i < num_objs; ++i) {
        MPD( ExportObj(i) );
        MP_DISPATCH( Convert( GetModel().obj(i) ) );
      }

    ////////////////////////// Algebraic constraints
    ifFltCon_ = 1;
    if (int n_cons = GetModel().num_algebraic_cons())
      for (int i = 0; i < n_cons; ++i) {
        MPD( ExportAlgCon(i) );
        MP_DISPATCH( ConvertAlgCon( i ) );
      }

    ////////////////////////// Logical constraints
    ifFltCon_ = 1;
    if (int n_lcons = GetModel().num_logical_cons())
      for (int i = 0; i < n_lcons; ++i) {
        MPD( ExportLogCon(i) );
        MP_DISPATCH( ConvertLogicalCon( i ) );
      }

    /// We could have produced variable names
    /// when exporting NL model info
    CopyItemNames();

    /// Signal we are not flattening anything
    ifFltCon_ = -1;
  }

  /// Export common expression \a i.
  void ExportCommonExpr(int i) {
    if (GetFlatCvt().GetFileAppender().IsOpen()) {
      fmt::MemoryWriter wrt;
      {
        MiniJSONWriter jw(wrt);
        jw["NL_COMMON_EXPR_index"] = i;
        // We don't receive defvar names from AMPL
        auto vn = GetModel().GetVarNamer();
        jw["name"] = vn.dvname(i);
        auto ce = GetModel().common_expr(i);
        fmt::MemoryWriter w2;
        w2 << "var " << vn.dvname(i) << " = ";
        WriteExpr<typename ProblemType::ExprTypes>(
              w2, ce.linear_expr(), ce.nonlinear_expr(), vn);
        jw["printed"] = w2.c_str();
      }
      wrt.write("\n");                     // EOL
      GetFlatCvt().GetFileAppender().Append(wrt);
    }
  }

  /// Export objective \a i.
  void ExportObj(int i) {
    if (GetFlatCvt().GetFileAppender().IsOpen()) {
      fmt::MemoryWriter wrt;
      {
        MiniJSONWriter jw(wrt);
        jw["NL_OBJECTIVE_index"] = i;
        jw["name"] = GetModel().obj_name(i);
        auto obj = GetModel().obj(i);
        jw["sense"] = (int)obj.type();
        fmt::MemoryWriter w2;
        w2 << (obj::MAX==obj.type() ? "maximize " : "minimize ");
        w2 << GetModel().obj_name(i) << ": ";
        WriteExpr<typename ProblemType::ExprTypes>(
              w2, obj.linear_expr(), obj.nonlinear_expr(),
              GetModel().GetVarNamer());
        jw["printed"] = w2.c_str();
      }
      wrt.write("\n");                     // EOL
      GetFlatCvt().GetFileAppender().Append(wrt);
    }
  }

  /// Export algebraic constraint \a i.
  void ExportAlgCon(int i) {
    if (GetFlatCvt().GetFileAppender().IsOpen()) {
      fmt::MemoryWriter wrt;
      {
        MiniJSONWriter jw(wrt);
        auto con = GetModel().algebraic_con(i);
        jw["NL_CON_TYPE"] = (con.nonlinear_expr() ? "nonlin" : "lin");
        jw["index"] = i;
        jw["name"] = GetModel().con_name(i);
        fmt::MemoryWriter w2;
        w2 << GetModel().con_name(i) << ": ";
        WriteAlgCon<typename ProblemType::ExprTypes>(
              w2, con, GetModel().GetVarNamer());
        jw["printed"] = w2.c_str();
      }
      wrt.write("\n");                     // EOL
      GetFlatCvt().GetFileAppender().Append(wrt);
    }
  }

  /// Export logical constraint \a i.
  void ExportLogCon(int i) {
    if (GetFlatCvt().GetFileAppender().IsOpen()) {
      fmt::MemoryWriter wrt;
      {
        MiniJSONWriter jw(wrt);
        auto con = GetModel().logical_con(i);
        jw["NL_CON_TYPE"] = "logical";
        int i_actual = GetModel().num_algebraic_cons() + i;
        jw["index"] = i_actual;
        jw["name"] = GetModel().con_name(i_actual);
        fmt::MemoryWriter w2;
        w2 << GetModel().con_name(i_actual) << ": ";
        WriteExpr<typename ProblemType::ExprTypes>(
              w2, con.expr(), GetModel().GetVarNamer());
        jw["printed"] = w2.c_str();
      }
      wrt.write("\n");                     // EOL
      GetFlatCvt().GetFileAppender().Append(wrt);
    }
  }

  /// Convert variables
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
    GetCopyLink().AddEntry({
          GetValuePresolver().GetSourceNodes().GetVarValues().MakeSingleKey().
                             Add(lbs.size()),
          vnr });
  }

  /// Convert a common expr
  void Convert(typename ProblemType::MutCommonExpr ) {
    /// Converting on demand, see VisitCommonExpr
  }

  /// Convert an objective
  void Convert(typename ProblemType::MutObjective obj) {
    auto obj_src =              // source value node for this obj
        GetValuePresolver().GetSourceNodes().GetObjValues()().Add();
    GetCopyLink().AddEntry(
          {
            obj_src,
            GetValuePresolver().GetTargetNodes().GetObjValues()().Add() });
    /// After the CopyLink, add One2ManyLink for converted expressions.
    /// When postsolving, CopyLink is executed last and copies obj values.
    /// This should resolve the issue of the "max-out"
    /// value conflict resolution.
    /// Tested by suf_common/funcpieces_01_01_obj.mod.
    pre::AutoLinkScope<FlatConverterType> auto_link_scope{
      GetFlatCvt(), obj_src
    };
    auto le = ToLinTerms(obj.linear_expr());
    NumericExpr e = obj.nonlinear_expr();
    EExpr eexpr;
    if (e) {
      eexpr=MP_DISPATCH( Visit(e) );
      le.add(eexpr.GetLinTerms());
      if (std::fabs(eexpr.constant_term())!=0.0) {
        /// Not using objective constant, should we?
        le.add_term(1.0, MakeFixedVar(eexpr.constant_term()));
      }
    }
    /// Sort/merge terms, otherwise we lose repeated terms
    /// in Gurobi where we just set 'obj attributes'
    /// to variables
    le.sort_terms();
    eexpr.GetQPTerms().sort_terms();
    /// Propagate context
    auto ctx = obj::MAX==obj.type() ? Context::CTX_POS : Context::CTX_NEG;
    GetFlatCvt().PropagateResult2LinTerms(le,
                                          GetFlatCvt().MinusInfty(),
                                          GetFlatCvt().Infty(), ctx);
    GetFlatCvt().PropagateResult2QuadTerms(eexpr.GetQPTerms(),
                                           GetFlatCvt().MinusInfty(),
                                           GetFlatCvt().Infty(), ctx);
    /// Add linear / quadratic obj
    LinearObjective lo { obj.type(),
          std::move(le.coefs()), std::move(le.vars()) };
    GetFlatCvt().AddObjective(
          QuadraticObjective{std::move(lo),
                             std::move(eexpr.GetQPTerms())});
  }

  /// Convert an algebraic constraint
  void ConvertAlgCon(int i) {
    pre::AutoLinkScope<FlatConverterType> auto_link_scope{
      GetFlatCvt(),
      GetValuePresolver().GetSourceNodes().GetConValues()().
          Add()           // Just add next node -
    };                    // assume the constraint order in NL
    AddAlgebraicConstraint( PrepareAlgConstraint(i) );
  }

  /// Algebraic constraint flattening preparation info
  struct AlgConPrepare {
    LinTerms lt;
    QuadTerms qt;
    double lb, ub;
    /// Complementarity
    double const_term;
    int compl_var;
  };

  /// Prepare alg constraint for moving into FlatCvt
  AlgConPrepare PrepareAlgConstraint(int i) {
    AlgConPrepare pre_result;
    auto con = GetModel().algebraic_con(i);
    pre_result.lt = ToLinTerms(con.linear_expr());
    EExpr ee;
    if (NumericExpr e = con.nonlinear_expr()) {
      ee=MP_DISPATCH( Visit(e) );
      pre_result.lt.add(ee.GetLinTerms());
    }
    pre_result.compl_var = GetModel().GetComplementarityVariable(i)-1;  // -1
    pre_result.const_term = ee.constant_term();
    if (pre_result.compl_var<0) {                // no complementarity
      pre_result.lb = con.lb() - pre_result.const_term;
      pre_result.ub = con.ub() - pre_result.const_term;
    } else {
      if (std::isfinite(con.lb())) {
        assert(!std::isfinite(con.ub()));
        pre_result.const_term -= con.lb();
      } else if (std::isfinite(con.ub())) {
        assert(!std::isfinite(con.lb()));
        pre_result.const_term -= con.ub();
      }
    }
    pre_result.lt.sort_terms();
    pre_result.qt = std::move(ee.GetQPTerms());        // quadratic terms, if any
    pre_result.qt.sort_terms();
    return pre_result;
  }

  /// Add algebraic constraint to FlatCvt
  void AddAlgebraicConstraint(AlgConPrepare&& pr) {
    if (pr.qt.empty()) {
      if (pr.compl_var<0)
        AddConstraint_AS_ROOT( LinConRange{ std::move(pr.lt),
                            { pr.lb, pr.ub }} );
      else
        AddConstraint_AS_ROOT(
              ComplementarityLinear{
                AffineExpr(std::move(pr.lt), pr.const_term), pr.compl_var } );
    } else {
      if (pr.compl_var<0)
        AddConstraint_AS_ROOT( QuadConRange{
                              { std::move(pr.lt), std::move(pr.qt) },
                              { pr.lb, pr.ub }} );
      else
        AddConstraint_AS_ROOT(
              ComplementarityQuadratic{
                QuadraticExpr
                { { std::move(pr.lt), std::move(pr.qt) }, pr.const_term },
                pr.compl_var } );
    }
  }

  /// Convert a logical constraint
  void ConvertLogicalCon(int i) {
    pre::AutoLinkScope<FlatConverterType> auto_link_scope{
      GetFlatCvt(),
      GetValuePresolver().GetSourceNodes().GetConValues()().
          Add()           // Just add next node -
    };                    // assume the constraint order in NL
    auto e = GetModel().logical_con(i);
    const auto resvar = MP_DISPATCH( Convert2Var(e.expr()) );
    if (GetFlatCvt().is_fixed(resvar)) {
      if (0==GetFlatCvt().fixed_value(resvar)) {
        MP_INFEAS("Logical constraint _slogcon["
                  + std::to_string(i+1)
                  + "] is false");
      }
    } else {
      GetFlatCvt().FixAsTrue(resvar);
      assert(GetFlatCvt().HasInitExpression(resvar));
    }
  }

  void CopyItemNames() {
    GetFlatCvt().AddVarNames(GetModel().var_names());
    GetFlatCvt().AddConNames(GetModel().con_names());
    GetFlatCvt().AddObjNames(GetModel().obj_names());
  }


protected:
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
  AffineExpr Convert2VarAsAffineExpr(EExpr&& ee) {
    return AffineExpr::Variable{Convert2Var(std::move(ee))};
  }
  AffineExpr Convert2AffineExpr(EExpr&& ee) {
    if (ee.is_affine())
      return MoveOutAffineExpr(std::move(ee));
    return Convert2VarAsAffineExpr(std::move(ee));
  }

  /// Generic functional expression array visitor.
  /// Assumes the arguments should be converted to variables.
  /// Can produce a new result variable/expression and
  /// specified constraints (normally, \a FuncConstraint) on it.
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

  /// Generic relational expression visitor.
  /// Can produce a new variable/expression and specified constraints on it.
  /// Produces a conditional linear/quadratic constraint
  /// b==1 <=/=> c'x [+ x'Qx] <=/= d.
  /// @tparam comp_kind: < (-2), <= (-1), == (0)
  /// @param ea: array of 2 expressions (comparison arguments lhs, rhs)
  template <int comp_kind, class ExprArray=std::initializer_list<Expr> >
  EExpr VisitRelationalExpression(ExprArray ea) {
    std::array<EExpr, 2> ee;
    Exprs2EExprs(ea, ee);
    ee[0].subtract(std::move(ee[1]));
    auto& lhs = ee[0];
    lhs.sort_terms();                             // to catch duplicates
    if (lhs.is_affine())                          // no QP terms
      return AssignResult2Args(                   // add conditional linear constraint
            ConditionalConstraint< LinConRhs<comp_kind> >
            { { std::move(lhs.GetLinTerms()),
                -lhs.constant_term() } } );
    return AssignResult2Args(                     // add conditional quadratic constraint
            ConditionalConstraint< QuadConRhs<comp_kind> >
            { { std::move(lhs.GetAlgConBody()),
                -lhs.constant_term() } } );
  }

  /// Disequality visitor.
  /// @param ea: array of 2 expressions (comparison arguments lhs, rhs)
  template <class ExprArray=std::initializer_list<Expr> >
  EExpr VisitDisequality(ExprArray ea) {
    std::array<EExpr, 2> ee;
    Exprs2EExprs(ea, ee);
    ee[0].subtract(std::move(ee[1]));
    auto& lhs = ee[0];
    lhs.sort_terms();                            // to catch duplicates
    if (lhs.is_affine()) {                       // no QP terms
      if (1==lhs.GetLinTerms().size()) {
        auto v = lhs.GetLinTerms().var(0);
        auto lb = GetFlatCvt().lb(v);
        auto ub = GetFlatCvt().ub(v);
        auto type = GetFlatCvt().var_type(v);
        if (var::Type::INTEGER==type
            && 1.0 == std::round(ub-lb)) {       // e.g., binary variable
          auto b = -lhs.constant_term()
              / lhs.GetLinTerms().coef(0);
          bool f1 = false;       // Could put this into Prepro for NEConstr
          if (lb == b) { b=ub; f1=true; }
          else if (ub == b) { b=lb; f1=true; }
          if (f1) {
            return AssignResult2Args(  // conditional linear constraint
                  ConditionalConstraint< LinConRhs<0> >
                    { { { {1.0}, {v} },
                        b } } );
          }  // else: should skip this constraint?
        }
      }
    }
    // General case: represent as: Not(lhs == 0)
    auto eq = (lhs.is_affine()) ?                // no QP terms
        AssignResult2Args(  // add conditional linear constraint
              ConditionalConstraint< LinConRhs<0> >
              { { std::move(lhs.GetLinTerms()),
                  -lhs.constant_term() } } ) :
        AssignResult2Args(  // add conditional quadratic constraint
              ConditionalConstraint< QuadConRhs<0> >
              { { std::move(lhs.GetAlgConBody()),
                  -lhs.constant_term() } } );
    assert(eq.is_variable());
    return AssignResult2Args(
          NotConstraint({eq.get_representing_variable()}));
  }

  /// Convert array of Expr's to array of EExpr's
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
public:          // need to be public due to CRTP
  EExpr VisitNumericConstant(NumericConstant n) {
    return EExpr::Constant{ n.value() };
  }

  EExpr VisitVariable(Reference r) {
    return EExpr::Variable{ r.index() };
  }

  EExpr VisitCommonExpr(Reference r) {
    const auto index = r.index();
    if (index >= (int)common_exprs_.size()) {
      assert(index < GetModel().num_common_exprs());
      common_exprs_.resize(
                  GetModel().num_common_exprs(),
                  -1);          // init by -1, "no variable"
    }
    if (common_exprs_[index]<0) {                 // not yet converted
      auto ce = MP_DISPATCH( GetModel() ).common_expr(index);
      EExpr eexpr( ToLinTerms(ce.linear_expr()) );
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

  EExpr VisitDiv(BinaryExpr e) {
    return VisitFunctionalExpression<DivConstraint>(
          { e.lhs(), e.rhs() });
  }

  EExpr VisitSum(typename BaseExprVisitor::SumExpr expr) {
    EExpr sum;              // Add up the elements' AffExpressions
    for (auto i =
         expr.begin(), end = expr.end(); i != end; ++i)
      sum.add( MP_DISPATCH( Convert2EExpr(*i) ) );
    return sum;
  }

  EExpr VisitMax(typename BaseExprVisitor::VarArgExpr e) {
    // Why need BaseExprVisitor:: here in g++ 9.2.1?
    return VisitFunctionalExpression<MaxConstraint>(e);
  }

  EExpr VisitMin(typename BaseExprVisitor::VarArgExpr e) {
    return VisitFunctionalExpression<MinConstraint>(e);
  }

  EExpr VisitAbs(UnaryExpr e) {
    return VisitFunctionalExpression<AbsConstraint>({ e.arg() });
  }

  EExpr VisitLogicalConstant(LogicalConstant c) {
    return EExpr::Constant{ double(c.value()) };
  }

  EExpr VisitEQ(RelationalExpr e) {
    return VisitRelationalExpression<0>({ e.lhs(), e.rhs() });
  }

  EExpr VisitNE(RelationalExpr e) {
    return VisitDisequality({ e.lhs(), e.rhs() });
  }

  EExpr VisitLE(RelationalExpr e) {
    return VisitRelationalExpression<-1>({ e.lhs(), e.rhs() });
  }

  EExpr VisitLT(RelationalExpr e) {
    return VisitRelationalExpression<-2>({ e.lhs(), e.rhs() });
  }

  EExpr VisitGE(RelationalExpr e) {
    return VisitRelationalExpression<1>({ e.lhs(), e.rhs() });
  }

  EExpr VisitGT(RelationalExpr e) {
    return VisitRelationalExpression<2>({ e.lhs(), e.rhs() });
  }

  EExpr VisitNot(NotExpr e) {
    return VisitFunctionalExpression<NotConstraint>({ e.arg() });
  }

  EExpr VisitAnd(BinaryLogicalExpr e) {
    return VisitFunctionalExpression<AndConstraint>({ e.lhs(), e.rhs() });
  }
  EExpr VisitForAll(IteratedLogicalExpr e) {
    return VisitFunctionalExpression<AndConstraint>(e);
  }

  EExpr VisitOr(BinaryLogicalExpr e) {
    return VisitFunctionalExpression<OrConstraint>({ e.lhs(), e.rhs() });
  }
  EExpr VisitExists(IteratedLogicalExpr e) {
    return VisitFunctionalExpression<OrConstraint>(e);
  }

  EExpr VisitIf(IfExpr e) {
    return VisitFunctionalExpression<IfThenConstraint>({
                e.condition(), e.then_expr(), e.else_expr() });
  }

  EExpr VisitImplication(ImplicationExpr e) {
    return VisitFunctionalExpression<ImplicationConstraint>({
                e.condition(), e.then_expr(), e.else_expr() });
  }

  EExpr VisitIff(BinaryLogicalExpr e) {
    return VisitRelationalExpression<0>({ e.lhs(), e.rhs() });
  }

  EExpr VisitAllDiff(PairwiseExpr e) {
    if (expr::ALLDIFF != e.kind())
      MP_RAISE("NOT_ALLDIFF NOT IMPLEMENTED");
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
          PLConstraint::Parameters{
                                  PLSlopes{
                                    breakpoints, slopes, 0.0, 0.0 } } ) );
  }


  ////////////////////////////////////////////////////
  /////////////// NONLINEAR FUNCTIONS ////////////////
  ////////////////////////////////////////////////////
  EExpr VisitPowConstExp(BinaryExpr e) {
    auto c = Cast<NumericConstant>(e.rhs()).value();
    if (2.0==c && IfQuadratizePow2()) {
      auto el = Convert2EExpr(e.lhs());
      return QuadratizeOrLinearize(el, el);
    }
    return AssignResult2Args( PowConstraint(
      PowConstraint::Arguments{ Convert2Var(e.lhs()) },
      PowConstraint::Parameters{ c } ) );
  }

  EExpr VisitPow2(UnaryExpr e) {
    if (IfQuadratizePow2()) {
      auto el = Convert2EExpr(e.arg());
      return QuadratizeOrLinearize(el, el);
    }
    return AssignResult2Args( PowConstraint(
      PowConstraint::Arguments{ Convert2Var(e.arg()) },
      PowConstraint::Parameters{ 2.0 } ) );
  }

  EExpr VisitPow(BinaryExpr e) {
    auto el = Convert2EExpr(e.lhs());
    auto er = Convert2EExpr(e.rhs());
    if (er.is_constant()) {
      if (2.0==er.constant_term() &&
          IfQuadratizePow2()) {
        return QuadratizeOrLinearize(el, el);
      }
      return AssignResult2Args(
            PowConstraint(
              PowConstraint::Arguments{ Convert2Var(std::move(el)) },
              PowConstraint::Parameters{ er.constant_term() } ) );
    }
    else if (el.is_constant())
      return VisitPowConstBase(e);
    else
      MP_RAISE("Unsupported: operator ^ with variable base and exponent");
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

  EExpr VisitAsin(UnaryExpr e) {
    return VisitFunctionalExpression<AsinConstraint>({ e.arg() });
  }

  EExpr VisitAcos(UnaryExpr e) {
    return VisitFunctionalExpression<AcosConstraint>({ e.arg() });
  }

  EExpr VisitAtan(UnaryExpr e) {
    return VisitFunctionalExpression<AtanConstraint>({ e.arg() });
  }

  EExpr VisitSinh(UnaryExpr e) {
    return VisitFunctionalExpression<SinhConstraint>({ e.arg() });
  }

  EExpr VisitCosh(UnaryExpr e) {
    return VisitFunctionalExpression<CoshConstraint>({ e.arg() });
  }

  EExpr VisitTanh(UnaryExpr e) {
    return VisitFunctionalExpression<TanhConstraint>({ e.arg() });
  }

  EExpr VisitAsinh(UnaryExpr e) {
    return VisitFunctionalExpression<AsinhConstraint>({ e.arg() });
  }

  EExpr VisitAcosh(UnaryExpr e) {
    return VisitFunctionalExpression<AcoshConstraint>({ e.arg() });
  }

  EExpr VisitAtanh(UnaryExpr e) {
    return VisitFunctionalExpression<AtanhConstraint>({ e.arg() });
  }


protected:         // More utilities
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
  /// Moreover if true, we assume they are for linearized PL
  void ConvertSOSCollection(ArrayRef<int> sosno, ArrayRef<double> ref,
                            bool fAllSOS2) {
    assert(sosno.size() == ref.size());
    std::map< int, std::multimap< double, int > > sos_map;
    for (auto i=ref.size(); i--; )
      if (sosno[i] &&
          ref[i]) {        // AMPL 2022 can produce SOS groups with all weights 0
        auto& sos_group = sos_map[sosno[i]];
        if (sos_group.end() != sos_group.find(ref[i])) {
          GetFlatCvt().AddWarning( "SOS_repeated_weight",
                                   "An SOS/SOS2 constraint has repeated weights, "
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
        AddConstraint(
              SOS2Constraint(vars, weights,
                             fAllSOS2 ?    // for linearized PL
                               SOSExtraInfo{{1.0, 1.0}} :
                               SOSExtraInfo{},
                             (fAllSOS2 ? "SOS2_PL_" : "SOS2_")
                             + std::to_string(group.first)
                             + '_'));
      else
        AddConstraint(
              SOS1Constraint(vars, weights,
                             "SOS1_"
                             + std::to_string(group.first)
                             + '_'));
    }
  }

  /// Depending on the target backend,
  /// can either convert to factor^2 (if el==er)
  /// or leave as quadratics.
  /// Currently only quadratize higher-order products.
  /// Can change arguments (move out).
  /// They could point to the same expr.
  /// PERFORMANCE WARNING: allows &el==&er, needed from Pow2
  EExpr QuadratizeOrLinearize(EExpr& el, EExpr& er) {
    if (!el.is_affine() && !er.is_constant())
      el = Convert2AffineExpr(std::move(el));      // will convert to a new var now
    if (!er.is_affine() && !el.is_constant())
      er = Convert2AffineExpr(std::move(er));
    if (!IfQuadratizePow2() &&
        !er.is_constant() && !el.is_constant() &&
        er.GetLinTerms().size() == el.GetLinTerms().size()) {
      const auto& ellt = el.GetLinTerms();
      const auto& erlt = er.GetLinTerms();
      if (1 == erlt.size() &&    // same variable in el and er
          0.0 == er.constant_term() && 0.0 == el.constant_term() &&
          ellt.var(0) == erlt.var(0)) {
        return Convert2Pow2(ellt, erlt);
      }
      el.sort_terms();
      er.sort_terms();
      if  (el == er) {            // Convert expr*expr to expr^2
        return Convert2Pow2(std::move(el));
      }
    }  // Otherwise, we proceed to store proper multiplication,
    // unless the result is affine
    if (!IfMultOutQPTerms() &&
        !er.is_constant() && !el.is_constant() ) {
      // Create a separate QC with this product.
      // This is handy if we are walking the objective,
      // as MIPFlatCvt only linearizes QC.
      return DontMultOut(std::move(el), std::move(er));
    }
    return MultiplyOut(el, er);   // Quadratize
  }

  /// (9*x) * x -> 9 x^2
  EExpr Convert2Pow2(const LinTerms& ellt, const LinTerms& erlt) {
    assert(1 == erlt.size() &&    // same variable in el and er
           ellt.var(0) == erlt.var(0));
    auto coef = ellt.coef(0) * erlt.coef(0);
    auto pow2var = GetFlatCvt().AssignResultVar2Args(
          PowConstraint{ {{ellt.var(0)}}, {2.0} });
    return { coef, pow2var };
  }

  /// aff_expr * aff_expr -> aff_expr^2
  EExpr Convert2Pow2(EExpr&& el) {
    auto affexpr2var = Convert2Var( std::move(el) );
    auto pow2var = GetFlatCvt().AssignResultVar2Args(
          PowConstraint{ {{affexpr2var}}, {2.0} });
    return EExpr::Variable{ pow2var };
  }

  /// Create product without multiplying out.
  /// Create a separate QC.
  EExpr DontMultOut(EExpr&& el, EExpr&& er) {
    const auto& ellt = el.GetLinTerms();
    const auto& erlt = er.GetLinTerms();
    if (1 == ellt.size() && 1 == erlt.size() &&   // a variable in el and er
        0.0 == er.constant_term() && 0.0 == el.constant_term()) {
      auto coef = ellt.coef(0) * erlt.coef(0);
      auto qc_res = GetFlatCvt().AssignResultVar2Args(
            QuadraticFunctionalConstraint{ { {      // = x*y+0
              LinTerms{},
              QuadTerms{ {1.0}, {ellt.var(0)}, {erlt.var(0)} }
            }, 0.0 } });
      return { coef, qc_res };
    }
    el.sort_terms();
    er.sort_terms();
    auto qc_res = GetFlatCvt().AssignResultVar2Args(
          QuadraticFunctionalConstraint{ { {      // = el*er+0
            LinTerms{},
            QuadTerms{ {1.0},
                       {Convert2Var( std::move(el) )},
                       {Convert2Var( std::move(er) )} }
          }, 0.0 } });
    return { 1.0, qc_res };
  }

  /// Multiply out two EEXprs
  EExpr MultiplyOut(const EExpr& el, const EExpr& er) {
    assert((el.is_affine() && er.is_affine()) ||
           (el.is_constant() || er.is_constant()));
    EExpr result;
    if (0.0!=std::fabs(er.constant_term())) {
      result.GetLinTerms().add(el.GetLinTerms());  // no const here
      result.GetLinTerms() *= er.constant_term();
      result.GetQPTerms().add(el.GetQPTerms());
      result.GetQPTerms() *= er.constant_term();
    }
    if (0.0!=std::fabs(el.constant_term())) {
      {
        auto ae2 = er.GetLinTerms();
        ae2 *= el.constant_term();
        result.GetLinTerms().add(ae2);
        result.constant_term(
              er.constant_term() * el.constant_term());
      }
      result.GetQPTerms().add(er.GetQPTerms());
      result.GetQPTerms() *= el.constant_term();
    }
    const auto& ae1 = el.GetLinTerms();
    const auto& ae2 = er.GetLinTerms();
    for (auto i1 = ae1.size(); i1--; ) {
      for (auto i2 = ae2.size(); i2--; ) {
        result.add_term(ae1.coef(i1) * ae2.coef(i2),
                           ae1.var(i1), ae2.var(i2) );
      }
    }
    result.sort_terms();      // eliminate 0's and duplicates
    return result;
  }


protected:
  //////////////////////// ADD CUSTOM CONSTRAINT ///////////////////////
  //////////////////////// Takes ownership /////////////////////////////

  /// Add constraint and propagate result into arguments.
  /// Use this when any variables can be the result of an expression
  /// (flat constraint)
  template <class Constraint>
  pre::NodeRange AddConstraint_AS_ROOT(Constraint&& con) {
    return GetFlatCvt().AddConstraint_AS_ROOT(std::move(con));
  }

  /// Add constraint, do not propagate result into arguments
  template <class Constraint>
  pre::NodeRange AddConstraint(Constraint&& con) {
    return GetFlatCvt().AddConstraint(std::move(con));
  }


  //////////////////////////// UTILITIES /////////////////////////////////
  ///
protected:
  pre::ValuePresolver& GetValuePresolver()
  { return GetFlatCvt().GetValuePresolver(); }


private:
  std::vector<int> common_exprs_;               // should be in FlatModel

  int ifFltCon_ = -1;   // -1: undefined, 0: walking an expr tree in an objective,
                        // 1: in a constraint

protected:
  /// Whether flattening a constraint vs an objective
  bool IfFlatteningAConstraint() const
  { assert(ifFltCon_>=0); return ifFltCon_; }

  //////////////////////////// CREATE OR FIND A FIXED VARIABLE //////////////////////////////
  int MakeFixedVar(double value)
  { return int(GetFlatCvt().MakeFixedVar(value)); }

  /// Presolve link just copying values between model items
  pre::CopyLink& GetCopyLink() { return GetFlatCvt().GetCopyLink(); }


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

  /// Distinguish between constraints and objectives.
  /// What about common expressions?
  int IfMultOutQPTerms() const {
    return IfFlatteningAConstraint() ?
          GetFlatCvt().IfPassQuadCon() : GetFlatCvt().IfPassQuadObj();
  }

  /// Quadratize Pow2 exactly when we pass QP terms
  bool IfQuadratizePow2() const { return IfMultOutQPTerms(); }


public:
  void InitOptions() override {
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
  const ProblemType& GetModel() const override { return model_; }
  /// The model as input from NL
  ProblemType& GetModel() override { return model_; }
  /// The model as input from NL
  const ProblemType& GetInputModel() const { return GetModel(); }
  /// The model as input from NL
  ProblemType& GetInputModel() { return GetModel(); }

  const FlatConverter& GetFlatCvt() const { return flat_cvt_; }
  FlatConverter& GetFlatCvt() { return flat_cvt_; }


private:
  ProblemType model_;
  FlatConverter flat_cvt_;
};


/// A 'final' ProblemFlattener in a hierarchy
template <template <typename, typename, typename> class ProblemFlt,
          class Problem, class FlatCvt>
class ProblemFltImpl :
    public ProblemFlt<
      ProblemFltImpl<ProblemFlt, Problem, FlatCvt>, Problem, FlatCvt> {
  /// Typedef Base
  using Base = ProblemFlt<
    ProblemFltImpl<ProblemFlt, Problem, FlatCvt>, Problem, FlatCvt>;
public:
  /// Construct
  ProblemFltImpl(Env& e) : Base(e) { }
};


} // namespace mp

#endif // PROBLEM_FLATTENER_H
