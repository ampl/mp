#ifndef XPRESSMPMODELAPI_H
#define XPRESSMPMODELAPI_H

#include <memory>

#include "mp/env.h"
#include "xpresscommon.h"
#include "mp/flat/nl_expr/model_api_base.h"
#include <variant>

namespace mp {


  /// Simple class to store expressions in the format Xpress expects them
  /// Note that, for the expression based API, as we visit the expression tree, 
  /// we accumulate nodes reversed in respect to what Xpress needs.
  /// For that, the function "reverse" is provided
  class NLParams {
    std::vector<int> types_;
    std::vector<double> values_;
    int resultVar_;
  public:
    /// No result variable (for generic expressions)
    NLParams() : resultVar_(-1) {}
    /// resvar = expression
    NLParams(int resultVar) : resultVar_(resultVar) {}
    
    const int* resultVar() const { return &resultVar_; }
    void addMember(int tokentype, double value) {
      types_.push_back(tokentype);
      values_.push_back(value);
    }
    void addMember(std::pair<int, double> exp) {
      types_.push_back(exp.first);
      values_.push_back(exp.second);
    }
    void addMembers(const NLParams& p) {
      types_.insert(types_.end(), p.types_.begin(), p.types_.end());
      values_.insert(values_.end(), p.values_.begin(), p.values_.end());
    }
    static std::pair<int, double> constant(double v) {
      return { XPRS_TOK_CON, v };
    }
    static std::pair<int, double> var(int index) {
      return { XPRS_TOK_COL, index };
    }
    static std::pair<int, double> op(int op) {
      return { XPRS_TOK_OP, op };
    }
    static std::pair<int, double> func(int func) {
      return { XPRS_TOK_IFUN, func };
    }
    const int* types() const { return types_.data(); }
    const double* values() const { return values_.data(); }
    int size() const { return static_cast<int>(types_.size()); }
    void reverse() {
      std::reverse(types_.begin(), types_.end());
      std::reverse(values_.begin(), values_.end());
    }
  };

class XpressmpModelAPI :
    public XpressmpCommon, public EnvKeeper,
    public BasicExprModelAPI<XpressmpModelAPI, NLParams >
{
  using BaseModelAPI = BasicExprModelAPI<XpressmpModelAPI, NLParams >;
public:
  /// Construct
  XpressmpModelAPI(Env& e) : EnvKeeper(e) { }

  /// Class name
  static const char* GetTypeName() { return "XpressmpModelAPI"; }

  /// Called before problem input.
  /// Chanve to allocate storage
  void InitProblemModificationPhase(const FlatModelInfo*);
  /// After
  void FinishProblemModificationPhase();

  /// TODO Implement the following functions using the solver's API
  void AddVariables(const VarArrayDef& );
  void SetLinearObjective( int iobj, const LinearObjective& lo );
  /// Whether accepting quadratic objectives:
  /// 0 - no, 1 - convex, 2 - nonconvex
  static int AcceptsQuadObj() { return 1; }
  void SetQuadraticObjective(int iobj, const QuadraticObjective& qo);
  
  //////////////////////////// GENERAL CONSTRAINTS ////////////////////////////
  USE_BASE_CONSTRAINT_HANDLERS(BaseModelAPI)

  /// TODO For each suppoted constraint type, add the ACCEPT_CONSTRAINT macro
  /// and the relative AddConstraint function.
  /// Below some typical constraint handlers of a MIP solver.
  /// Further constraint types which could be handled natively by some solvers:
  /// - IndicatorConstraint(Lin/Quad)(LE/EQ/GE)
  /// - Multidirectional indicators Cond(Lin/Quad)Con(LT/LE/EQ/GE/GT), where
  ///   the implication direction (</==/>) depends in the context
  /// - Complementarity
  /// - Logical, counting, piecewise-linear constraints.
  /// See \a constraints_std.h and other drivers.

  /// Ask if the solver accepts non-convex quadratic constraints
  static constexpr bool AcceptsNonconvexQC() { return true; }

  /// If cvt:prod=7 (and not 5) default.
  /// Recommendation to return the opposite value as
  /// AcceptsNonconvexQC().
  static constexpr bool WantLogicalizedProd2Bin()
  { return !AcceptsNonconvexQC(); }

  /// Ask if the solver can recognize SOCP corner cases
  /// (non-std representations such as xy>=1, see tests)
  /// from quadratic representations
  static constexpr bool CanSOCPCornerCasesFromQC() { return true; }


  /// LinCon(LE/EQ/GE) should have 'Recommended' for all backends
  /// and have an implementation,
  /// or a conversion rule is needed in a derived FlatConverter
  ACCEPT_CONSTRAINT(LinConLE, Recommended, CG_Linear)
  void AddConstraint(const LinConLE& lc);
  ACCEPT_CONSTRAINT(LinConEQ, Recommended, CG_Linear)
  void AddConstraint(const LinConEQ& lc);
  ACCEPT_CONSTRAINT(LinConGE, Recommended, CG_Linear)
  void AddConstraint(const LinConGE& lc);

  ACCEPT_CONSTRAINT(QuadConLE, Recommended, CG_Quadratic)
  void AddConstraint(const QuadConLE& qc);
  ACCEPT_CONSTRAINT(QuadConEQ, Recommended, CG_Quadratic)
  void AddConstraint(const QuadConEQ& qc);
  ACCEPT_CONSTRAINT(QuadConGE, Recommended, CG_Quadratic)
  void AddConstraint(const QuadConGE& qc);
  void AddLinTerms(XPRSprob lp, const LinTerms& lt, double rhsc, const char typec); // for quadratics

  /// Linear indicator constraints can be used as
  /// auxiliary constraints for logical conditions.
  /// If not handled, the compared expressions need
  /// deducible finite bounds for a big-M redefinition.
  ACCEPT_CONSTRAINT(IndicatorConstraintLinLE, Recommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinLE& mc);
  ACCEPT_CONSTRAINT(IndicatorConstraintLinEQ, Recommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinEQ& mc);
  ACCEPT_CONSTRAINT(IndicatorConstraintLinGE, Recommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinGE& mc);

  /// SOS constraints can be used by AMPL for redefinition of
  /// piecewise-linear expressions.
  /// Set ``option pl_linearize 0;`` in AMPL if the solver
  /// supports PL natively.
  ACCEPT_CONSTRAINT(SOS1Constraint, Recommended, CG_SOS)
  void AddConstraint(const SOS1Constraint& cc);
  ACCEPT_CONSTRAINT(SOS2Constraint, Recommended, CG_SOS)
  void AddConstraint(const SOS2Constraint& cc);

  // GENERAL CONSTRAINTS
  // Helper function for general constraints
  template <class Args, class Params, class NumOrLogic, class Id>
  void addGenCon(
    const CustomFunctionalConstraint<Args, Params, NumOrLogic, Id>& c,
      int xpressConType, bool fMarkVarsBinary = false);
  ACCEPT_CONSTRAINT(AbsConstraint, Recommended, CG_General)
  void AddConstraint(const AbsConstraint& ac);
  ACCEPT_CONSTRAINT(MaxConstraint, Recommended, CG_General)
  void AddConstraint(const MaxConstraint& ac);
  ACCEPT_CONSTRAINT(MinConstraint, Recommended, CG_General)
  void AddConstraint(const MinConstraint& ac);
  ACCEPT_CONSTRAINT(OrConstraint, Recommended, CG_General)
  void AddConstraint(const OrConstraint& ac);
  /// @todo bug in 9.7.0
  ACCEPT_CONSTRAINT(AndConstraint, AcceptedButNotRecommended, CG_General)
  void AddConstraint(const AndConstraint& ac);
  
  #define GLOBAL_LEVEL Recommended // Since v9.5.0
  ACCEPT_CONSTRAINT(DivConstraint, GLOBAL_LEVEL, CG_General)
  void AddConstraint(const DivConstraint& cc);
  ACCEPT_CONSTRAINT(SinConstraint, GLOBAL_LEVEL, CG_General)
  void AddConstraint(const SinConstraint& cc);
  ACCEPT_CONSTRAINT(CosConstraint, GLOBAL_LEVEL, CG_General)
  void AddConstraint(const CosConstraint& cc); 
  ACCEPT_CONSTRAINT(TanConstraint, GLOBAL_LEVEL, CG_General)
  void AddConstraint(const TanConstraint& cc);
  ACCEPT_CONSTRAINT(AsinConstraint, GLOBAL_LEVEL, CG_General)
  void AddConstraint(const AsinConstraint& cc);
  ACCEPT_CONSTRAINT(AcosConstraint, GLOBAL_LEVEL, CG_General)
  void AddConstraint(const AcosConstraint& cc);
  ACCEPT_CONSTRAINT(AtanConstraint, GLOBAL_LEVEL, CG_General)
  void AddConstraint(const AtanConstraint& cc);

  ACCEPT_CONSTRAINT(PLConstraint, GLOBAL_LEVEL, CG_General)
  void AddConstraint(const PLConstraint& plc);

  ACCEPT_CONSTRAINT(LogConstraint, GLOBAL_LEVEL, CG_General)
  void AddConstraint(const LogConstraint& cc);
  void AddConstraint(const LogAConstraint& cc);

  ACCEPT_CONSTRAINT(PowConstExpConstraint, GLOBAL_LEVEL, CG_General)
    void AddConstraint(const PowConstExpConstraint& cc);
  ACCEPT_CONSTRAINT(ExpConstraint, GLOBAL_LEVEL, CG_General)
  void AddConstraint(const ExpConstraint& cc);
  ACCEPT_CONSTRAINT(ExpAConstraint, GLOBAL_LEVEL, CG_General)
  void AddConstraint(const ExpAConstraint& cc);

  ACCEPT_CONSTRAINT(SinhConstraint, GLOBAL_LEVEL, CG_General)
  void AddConstraint(const SinhConstraint& cc);
  ACCEPT_CONSTRAINT(CoshConstraint, GLOBAL_LEVEL, CG_General)
  void AddConstraint(const CoshConstraint& cc);
  ACCEPT_CONSTRAINT(TanhConstraint, GLOBAL_LEVEL, CG_General)
  void AddConstraint(const TanhConstraint& cc);

  ACCEPT_CONSTRAINT(AsinhConstraint, GLOBAL_LEVEL, CG_General)
  void AddConstraint(const AsinhConstraint& cc);
  ACCEPT_CONSTRAINT(AcoshConstraint, GLOBAL_LEVEL, CG_General)
  void AddConstraint(const AcoshConstraint& cc);
  ACCEPT_CONSTRAINT(AtanhConstraint, GLOBAL_LEVEL, CG_General)
  void AddConstraint(const AtanhConstraint& cc);


  /// Add an explicifier constraint of the form
  /// -resulvar + expr  (>=|<=|=) 0
  void AddGlobalConstraint(const NLParams& params, char type);
  void AddGlobalConstraint(int resultVar, int argumentVar, int functionId);


  //////////////////////////// EXPRESSION TREES ////////////////////////////
  /// Handle expression trees: inherit basic API
  USE_BASE_EXPRESSION_HANDLERS(BaseModelAPI)
  ACCEPT_EXPRESSION_INTERFACE(Recommended)

  /// Whether accepts NLObjective
  static int AcceptsNLObj() { return 0; }
  /// Make a constant expression.
  static Expr MakeConstantExpr(double v) {
    NLParams p;
    p.addMember(NLParams::constant(v));
    return p;
  }

  /// Make an empty expression.
  static Expr MakeEmptyExpr() { return MakeConstantExpr(0.0); }

  /// Make an expression representing variable \a v.
  static Expr MakeVarExpr(int v) {
    NLParams p;
    p.addMember(NLParams::var(v));
    return p;
  }
    /// GetVarExpression(\a i): expression representing variable 0<=i<n_var.
    /// Only called for 'nonlinear' variables.
    Expr GetVarExpression(int i) { return MakeVarExpr(i); }

  /// GetZeroExpr(): constant 0.0 expression.
  /// Can be used to represent empty expression in an NLConstraint.
    Expr GetZeroExpression() { return MakeConstantExpr(0.0); }

  ACCEPT_CONSTRAINT(NLConstraint, Recommended, CG_General)
    void AddConstraint(const NLConstraint& nl);
    /// NLAssignEQ: algebraic expression expicifier.
    /// Meaning: var == expr.
    /// @note Accessors: GetName(), GetExpression(nle), GetVariable(nle).
    ACCEPT_CONSTRAINT(NLAssignEQ, Recommended, CG_General)
    void AddConstraint(const NLAssignEQ& nle);
  /// NLAssignLE: algebraic expression expicifier in positive context.
  /// Meaning: var <= expr.
  /// @note Accessors: GetName(), GetExpression(nle), GetVariable(nle).
  ACCEPT_CONSTRAINT(NLAssignLE, Recommended, CG_General)
    void AddConstraint(const NLAssignLE& nle);
  /// NLAssignGE: algebraic expression expicifier in negative context.
  /// Meaning: var >= expr.
  /// @note Accessors: GetName(), GetExpression(nle), GetVariable(nle).
  ACCEPT_CONSTRAINT(NLAssignGE, Recommended, CG_General)
    void AddConstraint(const NLAssignGE& nle);

  template <class MPExpr> void 
    AppendLinAndConstTerms(Expr& ff, const MPExpr& nla);

  /// Add a row to the linear coefficients matrix, return the last added row number, 
  /// to be used in XPRSnlpaddformulas to add a formula to the latest added row
  template <class MPExpr> int addLinearRow(MPExpr nl, char type, double rhs, double *prange) {
    int start = 0;

    int size = GetLinSize(nl);
    std::vector<int> cols(size);
    std::vector<double> coefs(size);
    for (int i = 0; i < size; ++i)
    {
      cols[i] = GetLinVar(nl, i);
      coefs[i] = GetLinCoef(nl, i);
    }
    XPRESSMP_CCALL(XPRSaddrows(lp(), 1, size, &type, &rhs, prange, &start, cols.data(), coefs.data()));
    return NumLinCons() - 1;
  }

  /// Create an expression with one argument (e.g. sin(exp(x)))
  template <class MPExpr> NLParams CreateExpressionOneArg(MPExpr expr, int xpressfunc) {
    auto ex = GetArgExpression(expr, 0);
    NLParams exp;
    exp.addMember(NLParams::func(xpressfunc));
    exp.addMembers(ex);
    exp.addMember(XPRS_TOK_RB, 0);
    return exp;
  }
  /// Create an expression with multiple arguments (e.g. min(x,sin(y),exp(z)))
  template <class MPExpr> NLParams CreateExpressionNArgs(MPExpr expr, int xpressfunc) {
    int n = GetNumArguments(expr);
    
    NLParams exp;
    exp.addMember(NLParams::func(xpressfunc));
    for (int i = 0; i < n; i++)
      exp.addMembers(GetArgExpression(expr, 1));
    exp.addMember(XPRS_TOK_RB, 0);
    return exp;
  }

  /// @brief Accept NLAffineExpr.
  /// @note Use accessors, not methods;
  /// - GetLinSize(le), GetLinCoef(le, i), GetLinTerm(le, i);
  ///   GetConstTerm(le).
  ACCEPT_EXPRESSION(NLAffineExpression, Recommended);
  Expr AddExpression(const NLAffineExpression& le);

  /// Accept NLQuadExpr.
  /// @note Use accessors, not methods;
  /// - GetLinSize(le), GetLinCoef(le, i), GetLinTerm(le, i);
  ///   GetQuadSize(le), GetQuadCoef(le, i),
  ///   GetQuadTerm1(le, i), GetQuadTerm2(le, i);
  ///   GetConstTerm(le).
  ACCEPT_EXPRESSION(NLQuadExpression, Recommended);
  Expr AddExpression(const NLQuadExpression& le);

  /// Each expression can be accepted as a proper expression,
  /// or as a flat functional constraint var <=/==/>= expr
  /// (in this case, with variables as arguments).
  /// The equality/inequality type of the flat constraint is
  /// determied by GetContext().
  ///
  /// @note Use accessor: GetArgExpression(ee, 0)
  /// - don't use ...Expression's methods.
  ///
  /// Similar for other expression types.

  ACCEPT_EXPRESSION(DivExpression, Recommended)
    Expr AddExpression(const DivExpression&);

  ACCEPT_EXPRESSION(MinExpression, Recommended)
    Expr AddExpression(const MinExpression&);
  ACCEPT_EXPRESSION(MaxExpression, Recommended)
    Expr AddExpression(const MaxExpression&);
  ACCEPT_EXPRESSION(AbsExpression, Recommended)
    Expr AddExpression(const AbsExpression&);

  ACCEPT_EXPRESSION(PowConstExpExpression, Recommended)
    Expr AddExpression(const PowConstExpExpression&);

  ACCEPT_EXPRESSION(ExpExpression, Recommended)
    Expr AddExpression(const ExpExpression&);
  ACCEPT_EXPRESSION(LogExpression, Recommended)
    Expr AddExpression(const LogExpression&);
  ACCEPT_EXPRESSION(LogAExpression, Recommended)
    Expr AddExpression(const LogAExpression&);

  ACCEPT_EXPRESSION(SinExpression, Recommended)
    Expr AddExpression(const SinExpression&);
  ACCEPT_EXPRESSION(CosExpression, Recommended)
    Expr AddExpression(const CosExpression&);
  ACCEPT_EXPRESSION(TanExpression, Recommended)
    Expr AddExpression(const TanExpression&);
  ACCEPT_EXPRESSION(AsinExpression, Recommended)
    Expr AddExpression(const AsinExpression&);
  ACCEPT_EXPRESSION(AcosExpression, Recommended)
    Expr AddExpression(const AcosExpression&);
  ACCEPT_EXPRESSION(AtanExpression, Recommended)
    Expr AddExpression(const AtanExpression&);
private:
  std::vector<int> obj_ind_save_, qobj_ind1_save_, qobj_ind2_save_;
};

} // namespace mp

#endif // XPRESSMPMODELAPI_H
