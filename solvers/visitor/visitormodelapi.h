#ifndef VISITORMODELAPI_H
#define VISITORMODELAPI_H

#include "mp/env.h"
#include "visitorcommon.h"
// Include the following if you don't need expressions support
// #include "mp/flat/model_api_base.h"
// else include this:
#include "mp/flat/nl_expr/model_api_base.h"

namespace mp {


/// VisitorModelAPI.
/// @note For expression tree solvers,
///       inherit from public `BasicExprModelAPI<ModelAPIclass, expressionclass>` 
///       (see below),
///       if you don't need expressions inherit from `BasicFlatModelAPI`
class VisitorModelAPI :
    public VisitorCommon, public EnvKeeper,
    public BasicExprModelAPI<VisitorModelAPI, Solver::VExpr>
{

  /// If you don't need expression support, use this typedef
  // using BaseModelAPI = BasicFlatModelAPI;
  /// else, typedef main base class for expressions support
  using BaseModelAPI = BasicExprModelAPI<VisitorModelAPI, Solver::VExpr>;

  
public:
  /// Construct
  VisitorModelAPI(Env& e) : EnvKeeper(e) { }

  /// Class name
  static const char* GetTypeName() { return "VisitorModelAPI"; }

  /// If any driver options added from here
  void InitCustomOptions() { }

  /// Called before problem modification.
  /// @param fmi: current problem information.
  ///   Always available via GetFlatModelInfo().
  /// @note this is called before each phase of model modification
  ///   which can happen during iterative solving.
  void InitProblemModificationPhase(const FlatModelInfo* fmi);
  /// After
  void FinishProblemModificationPhase();

  /// TODO Implement adding variables
  void AddVariables(const VarArrayDef& );
  /// TODO Implement setting (also changing) a linear (part of the) objective
  void SetLinearObjective( int iobj, const LinearObjective& lo );
  /// Whether accepting quadratic objectives:
  /// 0 - no, 1 - convex, 2 - nonconvex
  static int AcceptsQuadObj() { return 1; }
  /// TODO Implement setting (also changing) a quadratic objective
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
  /// See \a constr_std.h and other drivers.


  /// The linear range constraint, if fully supported with basis info etc.
  ACCEPT_CONSTRAINT(LinConRange, Recommended, CG_Linear)
  void AddConstraint(const LinConRange& lc);

  /// LinCon(LE/EQ/GE) should have 'Recommended' for all backends
  /// and have an implementation,
  /// or a conversion rule is needed in a derived FlatConverter
  ACCEPT_CONSTRAINT(LinConLE, Recommended, CG_Linear)
  void AddConstraint(const LinConLE& lc);
  ACCEPT_CONSTRAINT(LinConEQ, Recommended, CG_Linear)
  void AddConstraint(const LinConEQ& lc);
  ACCEPT_CONSTRAINT(LinConGE, Recommended, CG_Linear)
  void AddConstraint(const LinConGE& lc);

  /// Ask if the solver accepts non-convex quadratic constraints
  static constexpr bool AcceptsNonconvexQC() { return false; }

  /// If cvt:prod=7 (and not 5) default.
  /// Recommendation to return the opposite value as
  /// AcceptsNonconvexQC().
  static constexpr bool WantLogicalizedProd2Bin()
  { return !AcceptsNonconvexQC(); }

  /// QuadConRange is optional.
  ACCEPT_CONSTRAINT(QuadConRange, Recommended, CG_Quadratic)
  void AddConstraint(const QuadConRange& qc);

  /// If using quadratics,
  /// QuadCon(LE/EQ/GE) should have 'Recommended'
  /// and have an implementation.
  ACCEPT_CONSTRAINT(QuadConLE, Recommended, CG_Quadratic)
  void AddConstraint(const QuadConLE& qc);
  ACCEPT_CONSTRAINT(QuadConEQ, Recommended, CG_Quadratic)
  void AddConstraint(const QuadConEQ& qc);
  ACCEPT_CONSTRAINT(QuadConGE, Recommended, CG_Quadratic)
  void AddConstraint(const QuadConGE& qc);

  /// SOCP flags
  ///////////////////////////////////////////////////////////////////////
  /// Ask if the solver can mix conic quadratic
  /// (entered via dedicated API) and direct quadratic constraints
  static constexpr bool CanMixConicQCAndQC() { return false; }

  /// Ask if the solver can recognize SOCP corner cases
  /// (non-std representations such as xy>=1, see tests)
  /// from quadratic representations
  static constexpr bool CanSOCPCornerCasesFromQC() { return false; }

  /// Cones: SOCP
  ACCEPT_CONSTRAINT(QuadraticConeConstraint, Recommended, CG_Conic)
  void AddConstraint(const QuadraticConeConstraint& qc);
  ACCEPT_CONSTRAINT(RotatedQuadraticConeConstraint, Recommended, CG_Conic)
  void AddConstraint(const RotatedQuadraticConeConstraint& qc);
  // Other cones
  ACCEPT_CONSTRAINT(ExponentialConeConstraint, Recommended, CG_Conic)
  void AddConstraint(const ExponentialConeConstraint& qc);


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


  /// Some non linear constraints.
  /// See constr_std.h for more.
  ACCEPT_CONSTRAINT(PLConstraint, Recommended, CG_General)
  void AddConstraint(const PLConstraint& cc);
  ACCEPT_CONSTRAINT(MaxConstraint, Recommended, CG_General)
    void AddConstraint(const MaxConstraint& mc);
  ACCEPT_CONSTRAINT(MinConstraint, Recommended, CG_General)
    void AddConstraint(const MinConstraint& mc);
  ACCEPT_CONSTRAINT(AbsConstraint, Recommended, CG_General)
    void AddConstraint(const AbsConstraint& absc);
  ACCEPT_CONSTRAINT(AndConstraint, Recommended, CG_General)
    void AddConstraint(const AndConstraint& cc);
  ACCEPT_CONSTRAINT(OrConstraint, Recommended, CG_General)
    void AddConstraint(const OrConstraint& mc);
  ACCEPT_CONSTRAINT(ExpConstraint, Recommended, CG_General)
    void AddConstraint(const ExpConstraint& cc);
  ACCEPT_CONSTRAINT(ExpAConstraint, Recommended, CG_General)
    void AddConstraint(const ExpAConstraint& cc);
  ACCEPT_CONSTRAINT(LogConstraint, Recommended, CG_General)
    void AddConstraint(const LogConstraint& cc);
  ACCEPT_CONSTRAINT(LogAConstraint, Recommended, CG_General)
    void AddConstraint(const LogAConstraint& cc);
  ACCEPT_CONSTRAINT(PowConstExpConstraint, Recommended, CG_General)
    void AddConstraint(const PowConstExpConstraint& cc);
  ACCEPT_CONSTRAINT(PowConstraint, Recommended, CG_General)
  void AddConstraint(const PowConstraint& cc);
  ACCEPT_CONSTRAINT(SinConstraint, Recommended, CG_General)
    void AddConstraint(const SinConstraint& cc);
  ACCEPT_CONSTRAINT(CosConstraint, Recommended, CG_General)
    void AddConstraint(const CosConstraint& cc);  
  ACCEPT_CONSTRAINT(TanConstraint, Recommended, CG_General)
    void AddConstraint(const TanConstraint& cc);

  /// For the expression API, see below
  /// Overall switch
  /// Handle expression trees: inherit basic API

  USE_BASE_EXPRESSION_HANDLERS(BaseModelAPI)
  ACCEPT_EXPRESSION_INTERFACE(Recommended);

  /// GetVarExpression(\a i): expression representing variable 0<=i<n_var.
/// Only called for 'nonlinear' variables.
  Expr GetVarExpression(int i) {
    // Expr is a typedef to the solver-specific expression class, in this case
    // Solver::VExpr
    return  Expr::makeVariable(i); }

  /// Can be used to represent empty expression in an NLConstraint.
  Expr GetZeroExpression() { return  Solver::VExpr::makeConstant(0); }

  /// Whether accepts NLObjective
  static int AcceptsNLObj() { return 1; }
  void SetNLObjective(int, const NLObjective&);

  /// Above which reference count,
  /// an algebraic formula node should be assigned to a variable.
  /// Should normally be INT_MAX for solvers
  /// using pointers to store expressions (SCIP)
  /// or introducing defined variables (MP2NL).
  /// Should be small for solvers
  /// using strings to represent formulas.
  /// Note: 0 means all nodes outlined.
  static int NLAssignLevelDefault() { return 1; }

  /// Same for logical expressions.
  static int NLReifLevelDefault() { return 1; }


  /// Top level general NL constraint
  /// lhs <= f(x) <= rhs
  ACCEPT_CONSTRAINT(NLConstraint, Recommended, CG_Algebraic)
  void AddConstraint(const NLConstraint& nlc);
  
  

  // TODO: Mandatory:
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



  /// NL logical constraint: expression = true.
  ///
  /// @note Should be 'Recommended'
  ///   whenever logical expressions are accepted.
  /// @note Use GetExpression(nll) to access the expression.
  /// @note Constraint group: CG_Nonlinear for SCIP,
  ///   because using SCIPcreateConsBasicNonlinear().
  ACCEPT_CONSTRAINT(NLLogical, Recommended, CG_Nonlinear)
    void AddConstraint(const NLLogical& nll);

  /// NL equivalence: expression <==> var.
  /// This is an 'expression explicifier'.
  ///
  /// @note Should be 'Recommended'
  ///   whenever logical expressions are accepted.
  /// @note Accessors: GetExpression(nle), GetVariable(nle).
  ACCEPT_CONSTRAINT(NLReifEquiv, Recommended, CG_Nonlinear)
  void AddConstraint(const NLReifEquiv& nle);
  /// NL implication: var==1 ==> expression.
  /// This is an expression explicifier in positive context.
  ///
  /// @note Should be 'Recommended'
  ///   whenever logical expressions are accepted.
  /// @note Accessors: GetExpression(nle), GetVariable(nle).
  ACCEPT_CONSTRAINT(NLReifImpl, Recommended, CG_Nonlinear)
  void AddConstraint(const NLReifImpl& nle);
  /// NL reverse implication: expression ==> var==1.
  /// This is an expression explicifier in negative context.
  ///
  /// @note Should be 'Recommended'
  ///   whenever logical expressions are accepted.
  /// @note Accessors: GetExpression(nle), GetVariable(nle).
  ACCEPT_CONSTRAINT(NLReifRimpl, Recommended, CG_Nonlinear)
  void AddConstraint(const NLReifRimpl& nle);

  

  /// @brief Accept NLAffineExpr.
  /// @note Use accessors, not methods;
  /// - GetLinSize(le), GetLinCoef(le, i), GetLinTerm(le, i);
  ///   GetConstTerm(le).
  ACCEPT_EXPRESSION(NLAffineExpression, Recommended);
  Expr AddExpression(const NLAffineExpression& le);

  ACCEPT_EXPRESSION(NLQuadExpression, Recommended);
  Expr AddExpression(const NLQuadExpression& le);

  ACCEPT_EXPRESSION(DivExpression, Recommended)
  Expr AddExpression(const DivExpression&);


  ACCEPT_EXPRESSION(LogExpression, Recommended)
  Expr AddExpression(const LogExpression&);
  ACCEPT_EXPRESSION(LogAExpression, Recommended)
  Expr AddExpression(const LogAExpression&);
  ACCEPT_EXPRESSION(ExpExpression, Recommended)
  Expr AddExpression(const ExpExpression&);
  ACCEPT_EXPRESSION(ExpAExpression, Recommended)
  Expr AddExpression(const ExpAExpression&);

  ACCEPT_EXPRESSION(PowExpression, Recommended)
  Expr AddExpression(const PowExpression&);
  ACCEPT_EXPRESSION(PowConstExpExpression, Recommended)
  Expr AddExpression(const PowConstExpExpression&);
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

  ACCEPT_EXPRESSION(SinhExpression, Recommended)
  Expr AddExpression(const SinhExpression&);
  ACCEPT_EXPRESSION(CoshExpression, Recommended)
  Expr AddExpression(const CoshExpression&);
  ACCEPT_EXPRESSION(TanhExpression, Recommended)
  Expr AddExpression(const TanhExpression&);

  ACCEPT_EXPRESSION(AsinhExpression, Recommended)
  Expr AddExpression(const AsinhExpression&);
  ACCEPT_EXPRESSION(AcoshExpression, Recommended)
  Expr AddExpression(const AcoshExpression&);
  ACCEPT_EXPRESSION(AtanhExpression, Recommended)
  Expr AddExpression(const AtanhExpression&);

  private:
    // Helper function to fetch variable names from the solver API stub,
    // used to print the model variable names
    std::function<std::string_view(int)> GetVarName;

    // Helper function to append an affine expression (linear + constant terms)
    // coming from MP (nla) to an expression in terms of this driver (ff)
    template <class MPExpr>
    void AppendLinAndConstTerms(Expr& ff, const MPExpr& nla);

    template <class MPExpr>
    void AppendQuadTerms(Expr&, const MPExpr&);
};

} // namespace mp

#endif // VISITORMODELAPI_H
