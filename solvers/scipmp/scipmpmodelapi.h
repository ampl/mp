#ifndef SCIPMODELAPI_H
#define SCIPMODELAPI_H

#include <climits>

#include "mp/env.h"
#include "scipmpcommon.h"
#include "mp/flat/nl_expr/model_api_base.h"

namespace mp {

class ScipModelAPI
    : public ScipCommon, public EnvKeeper,
      public BasicExprModelAPI<ScipModelAPI, SCIP_EXPR*>
{
  using BaseModelAPI = BasicExprModelAPI<ScipModelAPI, SCIP_EXPR*>;

protected:
  void linearHelper(const int* pvars, const double* pcoefs, const size_t size, const char* name, const double lb, const double ub);
  void helpIndicatorLin(const int* pvars, const double* pcoefs, const size_t size, const char* name, const double rhs, const int binary_var, const int binary_value, SCIP_Bool lessthanineq);
  void helpQuad(const char* name, const int* lt_pvars, const double* lt_pcoefs, const size_t lt_size, const int* qt_pvars1, const int* qt_pvars2, const double* qt_pcoefs, const size_t qt_size, const double lb, const double ub);

public:
  /// Construct
  ScipModelAPI(Env& e) : EnvKeeper(e) { }

  /// Class name
  static const char* GetTypeName() { return "ScipModelAPI"; }

  /// If any driver options added from the ModelAPI
  /// (usually only in the Backend)
  void InitCustomOptions() { }

  /// Called before problem input.
  /// Model info can be used to preallocate memory.
  void InitProblemModificationPhase(const FlatModelInfo*);
  /// After
  void FinishProblemModificationPhase();

  /// Implement the following functions using the solver's API
  void AddVariables(const VarArrayDef& );
  void SetLinearObjective( int iobj, const LinearObjective& lo );
  /// Whether accepting quadratic objectives:
  /// 0 - no, 1 - convex, 2 - nonconvex
  static int AcceptsQuadObj() { return 0; }
  void SetQuadraticObjective(int iobj, const QuadraticObjective& qo);

  /// Above which reference count,
  /// a formula node should be assigned to a variable.
  /// Should normally be INT_MAX for solvers
  /// using pointers to store expressions (SCIP),
  /// and positive for solvers using strings to represent formulas.
  static int NLAssignLevelDefault() { return INT_MAX; }

  //////////////////////////// GENERAL CONSTRAINTS ////////////////////////////
  /// Handle flat constraints: inherit basic API
  USE_BASE_CONSTRAINT_HANDLERS(BaseModelAPI)

  //////////////////////////// EXPRESSION TREES ////////////////////////////
  /// Hanlde expression trees: inherit basic API
  USE_BASE_EXPRESSION_HANDLERS(BaseModelAPI)

  /// ACCEPT_EXPRESSION_INTERFACE():
  /// 'NotAccepted' or
  /// 'AcceptedButNotRecommended' would resort to flat constraints uniformly
  /// which outline each expression with an auxiliary variable:
  /// e.g., ExpConstraint(a, b), meaning a = exp(b), with a, b variables.
  /// But the latter option would enable the solver option acc:_expr,
  /// which, when set to 1, switches on the expression trees.
  ///
  /// See also per-expression and per-constraint type switches
  /// (ACCEPT_EXPRESSION and ACCEPT_CONSTRAINT.)
  /// When both are enabled for certain function
  /// (e.g., CosExpression and CosConstraint),
  /// solver option acc:cos allows switching between them.
  ///
  /// Expression trees are submitted to the solver via
  /// the high-level constraints
  /// NLConstraint, NLAssignLE, NLAssignEQ, NLAssignGE,
  /// NLComplementarity,
  /// NLLogical, NLEquivalence, NLImpl, NLRimpl, and NLObjective.
  ACCEPT_EXPRESSION_INTERFACE(Recommended);

  /// Whether accepts NLObjective.
  /// No, as of SCIP 9.1.
  static int AcceptsNLObj() { return 0; }

  /// Once expressions are supported, need the following
  /// helper methods.
  ///
  /// GetVarExpression(\a i): expression representing variable 0<=i<n_var.
  /// Only called for 'nonlinear' variables.
  SCIP_EXPR* GetVarExpression(int i);

  /// GetZeroExpr(): constant 0.0 expression.
  /// Can be used to represent empty expression in an NLConstraint.
  SCIP_EXPR* GetZeroExpression();

  /// For each supported top-level constraint type,
  /// add the ACCEPT_CONSTRAINT macro
  /// and the corresponding AddConstraint function.
  /// Below some typical constraint handlers of a MIP solver.
  /// Further constraint types which could be handled natively by some solvers:
  /// - IndicatorConstraint(Lin/Quad)(LE/EQ/GE)
  /// - Multidirectional indicators Cond(Lin/Quad)Con(LT/LE/EQ/GE/GT), where
  ///   the implication direction (</==/>) depends on the context
  /// - Complementarity
  /// - Logical, counting, piecewise-linear constraints.
  /// See \a constr_std.h.

  /// Below are high-level flat algebraic constraints.

  /// The linear range constraint, if fully supported with basis info etc.
  ACCEPT_CONSTRAINT(LinConRange, Recommended, CG_Linear)
  void AddConstraint(const LinConRange& lc);

  /// LinCon(LE/EQ/GE) should have 'Recommended' for all non-expression backends
  /// and have an implementation,
  /// or a conversion rule is needed in a derived FlatConverter.
  /// Even for expression backends, they can be implemented for efficiency.
  ACCEPT_CONSTRAINT(LinConLE, Recommended, CG_Linear)
  void AddConstraint(const LinConLE& lc);
  ACCEPT_CONSTRAINT(LinConEQ, Recommended, CG_Linear)
  void AddConstraint(const LinConEQ& lc);
  ACCEPT_CONSTRAINT(LinConGE, Recommended, CG_Linear)
  void AddConstraint(const LinConGE& lc);

  /// Ask if the solver accepts non-convex quadratic constraints
  static constexpr bool AcceptsNonconvexQC() { return true; }

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
  /// Even for expression backends, they can be implemented for efficiency,
  /// if the solver supports pure variable-quadratics.
  ACCEPT_CONSTRAINT(QuadConLE, Recommended, CG_Quadratic)
  void AddConstraint(const QuadConLE& qc);
  ACCEPT_CONSTRAINT(QuadConEQ, Recommended, CG_Quadratic)
  void AddConstraint(const QuadConEQ& qc);
  ACCEPT_CONSTRAINT(QuadConGE, Recommended, CG_Quadratic)
  void AddConstraint(const QuadConGE& qc);

  /// If NLConstraint is accepted, with expression trees,
  /// then top-level nonlinear algebraic constraints
  /// are submitted to the solver via NLConstraint.
  ///
  /// @note If NLConstraint is not accepted, then LinCon(LE/GE/EQ/[Range])
  ///   should be. In such case, top-level algebraic expressions
  ///   are explicified via NLAssign(LE/GE/EQ),
  ///   for example, for Gurobi 12.
  ///
  /// @note It is still recommended that pure-linear
  ///   and pure-quadratic constraints are accepted, then they are
  ///   used to submit the corresponding constraint types
  ///   without expressions.
  ///
  /// @note To access information from an NLConstraint,
  ///   use the following accessors (don't use methods of NLConstraint itself):
  ///   - GetLinSize(nlc), GetLinCoef(nlc, i), GetLinVar(nlc, i),
  ///     GetExpression(nlc), GetLower(nlc), GetUpper(nlc).
  ACCEPT_CONSTRAINT(NLConstraint, Recommended, CG_Nonlinear)
  void AddConstraint(const NLConstraint& nlc);

  /// NLAssignEQ: algebraic expression expicifier.
  /// Meaning: var == expr.
  ///
  /// @note Should be 'Recommended' for expression solvers.
  /// @note Example: GRBaddgenconstrNL in Gurobi 12.
  ///   In other API types can be implemented using
  ///   NL algebraic constraint.
  ///
  /// @note Accessors: GetExpression(nle), GetVariable(nle).
  ACCEPT_CONSTRAINT(NLAssignEQ, Recommended, CG_Nonlinear)
  void AddConstraint(const NLAssignEQ& nle);
  /// NLAssignLE: algebraic expression expicifier in positive context.
  /// Meaning: var <= expr.
  ///
  /// @note Should be 'Recommended' in expression solvers.
  /// @note Can be implemented as NLAssignEQ,
  ///   but this may lose convexity.
  /// @note Accessors: GetExpression(nle), GetVariable(nle).
  ACCEPT_CONSTRAINT(NLAssignLE, Recommended, CG_Nonlinear)
  void AddConstraint(const NLAssignLE& nle);
  /// NLAssignGE: algebraic expression expicifier in negative context.
  /// Meaning: var >= expr.
  ///
  /// @note Should be 'Recommended' in expression solvers.
  /// @note Can be implemented as NLAssignEQ,
  ///   but this may lose convexity.
  /// @note Accessors: GetExpression(nle), GetVariable(nle).
  ACCEPT_CONSTRAINT(NLAssignGE, Recommended, CG_Nonlinear)
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
  
  /// Moreover, once algebraic expressions are accepted
  /// via NLConstraint, subexpressions might be submitted via
  /// NLAffineExpr and NLQuadExpr.

  /// @brief NLAffineExpr.
  /// @note Use accessors, not methods;
  /// - GetLinSize(le), GetLinCoef(le, i), GetLinTerm(le, i);
  ///   GetConstTerm(le).
  ACCEPT_EXPRESSION(NLAffineExpression, Recommended);
  SCIP_EXPR* AddExpression(const NLAffineExpression& le);

  /// NLQuadExpr.
  /// @note Use accessors, not methods;
  /// - GetLinSize(le), GetLinCoef(le, i), GetLinTerm(le, i);
  ///   GetQuadSize(le), GetQuadCoef(le, i),
  ///   GetQuadTerm1(le, i), GetQuadTerm2(le, i);
  ///   GetConstTerm(le).
  ACCEPT_EXPRESSION(NLQuadExpression, Recommended);
  SCIP_EXPR* AddExpression(const NLQuadExpression& le);

  /// Each expression can be accepted as a proper expression,
  /// or as a flat functional constraint var <=/==/>= expr
  /// (in this case, with variables as arguments).
  /// The uequality/inqeuality type of the flat constraint is
  /// determied by GetContext().
  ///
  /// @note For each expression,
  /// say ACCEPT_EXPRESSION(Recommended)
  /// and/or ACCEPT_EXPRESSION(AcceptedButNotRecommended).
  /// This can be user-configured via solver options 'acc:exp' etc.
  ///
  /// @note Use accessor: GetArgExpression(ee, 0)
  /// - don't ExpExpression's methods.
  ///
  /// Similar for other expression types.
  ACCEPT_EXPRESSION(AbsExpression, Recommended)
  SCIP_EXPR* AddExpression(const AbsExpression& absc);

  /// For each flat constraint type,
  /// say ACCEPT_CONSTRAINT(Recommended)
  /// and/or ACCEPT_CONSTRAINT(AcceptedButNotRecommended).
  /// This can be user-configured via solver options 'acc:exp' etc.
  ACCEPT_CONSTRAINT(AbsConstraint, Recommended, CG_General)
  void AddConstraint(const AbsConstraint& absc);

  // SCIP 9 has AND/OR as constraints only:
  // ACCEPT_EXPRESSION(AndExpression, AcceptedButNotRecommended)
  // void AddExpression(const AndExpression& cc);
  // ACCEPT_EXPRESSION(OrExpression, Recommended)
  // void AddExpression(const OrExpression& dc);
  ACCEPT_CONSTRAINT(AndConstraint, Recommended, CG_General)
  void AddConstraint(const AndConstraint& cc);
	ACCEPT_CONSTRAINT(OrConstraint, Recommended, CG_General)
  void AddConstraint(const OrConstraint& dc);

  /// Linear indicator constraints can be used as
  /// auxiliary constraints for logical conditions.
  /// If not handled, the compared expressions need
  /// deducible finite bounds for a big-M redefinition.
  ///
  /// @note Used to be 'AcceptedButNotRecommended':
  ///   not sure about performance (vs big-M) before SCIP 9.2.2.
  ACCEPT_CONSTRAINT(IndicatorConstraintLinLE, Recommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinLE& mc);
  ACCEPT_CONSTRAINT(IndicatorConstraintLinEQ, Recommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinEQ& mc);
  ACCEPT_CONSTRAINT(IndicatorConstraintLinGE, Recommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinGE& mc);

  /// Cones
	ACCEPT_CONSTRAINT(QuadraticConeConstraint, AcceptedButNotRecommended, CG_Conic)
	void AddConstraint(const QuadraticConeConstraint& qc);

  /// SOS constraints can be used by AMPL for redefinition of
  /// piecewise-linear expressions.
  /// @note Set ``option pl_linearize 0;`` in AMPL if the solver
  ///   supports PL natively.
  ACCEPT_CONSTRAINT(SOS1Constraint, Recommended, CG_SOS)
  void AddConstraint(const SOS1Constraint& cc);
	ACCEPT_CONSTRAINT(SOS2Constraint, Recommended, CG_SOS)
  void AddConstraint(const SOS2Constraint& cc);

  /// SCIP nonlinear general constraints and expressions.

  ACCEPT_EXPRESSION(ExpExpression, Recommended)
  SCIP_EXPR* AddExpression(const ExpExpression& );
  ACCEPT_CONSTRAINT(ExpConstraint, Recommended, CG_General)
  void AddConstraint(const ExpConstraint& cc);

  ACCEPT_EXPRESSION(LogExpression, Recommended)
  SCIP_EXPR* AddExpression(const LogExpression& );
  ACCEPT_EXPRESSION(LogAExpression, Recommended)
  SCIP_EXPR* AddExpression(const LogAExpression&);

  ACCEPT_CONSTRAINT(LogAConstraint, Recommended, CG_General)
  void AddConstraint(const LogAConstraint& cc);
  ACCEPT_CONSTRAINT(LogConstraint, Recommended, CG_General)
  void AddConstraint(const LogConstraint& cc);

  /// @note Use accessor: GetParameter(pe, 0)
  ///   - don't use PowConstExpExpression's methods.
  ACCEPT_EXPRESSION(PowConstExpExpression, Recommended)
  SCIP_EXPR* AddExpression(const PowConstExpExpression& );
  ACCEPT_CONSTRAINT(PowConstExpConstraint, Recommended, CG_General)
  void AddConstraint(const PowConstExpConstraint& cc);

  ACCEPT_EXPRESSION(SinExpression, Recommended)
  SCIP_EXPR* AddExpression(const SinExpression& );
  ACCEPT_CONSTRAINT(SinConstraint, Recommended, CG_General)
  void AddConstraint(const SinConstraint& cc);

  ACCEPT_EXPRESSION(CosExpression, Recommended)  //pretty slow in SCIP 8/9
  SCIP_EXPR* AddExpression(const CosExpression& );
  ACCEPT_CONSTRAINT(CosConstraint, Recommended, CG_General)
  void AddConstraint(const CosConstraint& cc);

  // TODO Div; PowVarVar;
  // CondLin... - not really, reader_nl.cpp only handles bool args
};

} // namespace mp

#endif // SCIPMODELAPI_H
