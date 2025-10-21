#ifndef COPTMODELAPI_H
#define COPTMODELAPI_H

#include <memory>

#include "mp/env.h"
#include "coptcommon.h"
#include "mp/flat/nl_expr/model_api_base.h"


namespace mp {
    /// Simple class to store expressions in the format Copt expects them
      /// Note that, for the expression based API, as we visit the expression tree, 
      /// we accumulate nodes reversed in respect to what Copt needs.
      /// For that, the function "reverse" is provided
    class NLParams {
        std::vector<int> typesOrVars_;
        std::vector<double> values_;
        std::vector<int> linearIndices_;
        std::vector<double> linearCoeffs_;
       
        int resultVar_;
    public:
        /// No result variable (for generic expressions)
        NLParams() : resultVar_(-1) {}
        /// resvar = expression
        NLParams(int resultVar) : resultVar_(resultVar) {}

        const int* resultVar() const { return &resultVar_; }
        
        void reserveLinear(int size) {
            linearCoeffs_.reserve(size);
            linearIndices_.reserve(size);
        }
        void addLinear(int index, double coeff) {
            linearCoeffs_.push_back(coeff);
            linearIndices_.push_back(index);
        }
        void addMembers(const NLParams& p) {
            typesOrVars_.insert(typesOrVars_.end(), p.typesOrVars_.begin(), p.typesOrVars_.end());
            values_.insert(values_.end(), p.values_.begin(), p.values_.end());
        }
        void addConstant(double v) {
            typesOrVars_.push_back(COPT_NL_GET);
            values_.push_back(v);
        }
        void addVar(int index) {
            assert(index >= 0);
            typesOrVars_.push_back(index);
        }
        void addVar(const NLParams& exp) {
            typesOrVars_.push_back(exp.typesOrVars_[0]);
        }
        void addOp(int op) {
            assert(op < 0);
            typesOrVars_.push_back(op);
        }
    
        int index() const { return linearIndices_[0]; }
        const int* linearIndices() const { return linearIndices_.data(); }
        const double* linearCoeffs() const { return linearCoeffs_.data(); }
        const int nLinear() const { return linearIndices_.size(); }
        const int* tokens() const { return typesOrVars_.data(); }
        const double* tokenElements() const { return values_.data(); }
        int nTokens() const { return static_cast<int>(typesOrVars_.size()); }
        int nTokenElements() const { return static_cast<int>(values_.size());}
        void reverse() {
            std::reverse(typesOrVars_.begin(), typesOrVars_.end());
            std::reverse(values_.begin(), values_.end());
        }
    };
class CoptModelAPI :
    public CoptCommon, public EnvKeeper,
    public BasicExprModelAPI<CoptModelAPI, NLParams >
{
    using BaseModelAPI = BasicExprModelAPI<CoptModelAPI, NLParams >;

public:
  /// Construct
  CoptModelAPI(Env& e) : EnvKeeper(e) { }

  /// Class name
  static const char* GetTypeName() { return "CoptModelAPI"; }

  /// Called before problem input
  void InitProblemModificationPhase(const FlatModelInfo*);
  /// After
  void FinishProblemModificationPhase();

  void AddVariables(const VarArrayDef& );
  void SetLinearObjective( int iobj, const LinearObjective& lo );
  /// Whether accepting quadratic objectives:
  /// 0 - no, 1 - convex, 2 - nonconvex
  static int AcceptsQuadObj() { return 1; }
  void SetQuadraticObjective(int iobj, const QuadraticObjective& qo);

  //////////////////////////// GENERAL CONSTRAINTS ////////////////////////////
  USE_BASE_CONSTRAINT_HANDLERS(BaseModelAPI)

  /// COPT prefers linear ranges
  ACCEPT_CONSTRAINT(LinConRange, Recommended, CG_Linear)
  void AddConstraint(const LinConRange& lc);

  /// LinCon(LE/EQ/GE) should have 'Recommended' for all backends
  /// and have an implementation
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

  ACCEPT_CONSTRAINT(IndicatorConstraintLinLE, Recommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinLE& mc);
  ACCEPT_CONSTRAINT(IndicatorConstraintLinEQ, Recommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinEQ& mc);
  ACCEPT_CONSTRAINT(IndicatorConstraintLinGE, Recommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinGE& mc);

  ACCEPT_CONSTRAINT(SOS1Constraint, Recommended, CG_SOS)
  void AddConstraint(const SOS1Constraint& cc);
  ACCEPT_CONSTRAINT(SOS2Constraint, Recommended, CG_SOS)
  void AddConstraint(const SOS2Constraint& cc);

  //////////////////////////// EXPRESSION TREES ////////////////////////////
  /// Handle expression trees: inherit basic API
  USE_BASE_EXPRESSION_HANDLERS(BaseModelAPI)
  ACCEPT_EXPRESSION_INTERFACE(AcceptedButNotRecommended)

  /// Whether accepts NLObjective
  static int AcceptsNLObj() { return 0; }
  /// Make a constant expression.
  static Expr MakeConstantExpr(double v) {
      NLParams p;
      p.addConstant(v);
      return p;
  }

  /// Make an empty expression.
  static Expr MakeEmptyExpr() { return MakeConstantExpr(0.0); }

  /// Make an expression representing variable \a v.
  static Expr MakeVarExpr(int v) {
      NLParams p;
      p.addVar(v);
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

  
  void AddGlobalConstraint(const NLParams& params, char type);
  /// Create an expression with one argument (e.g. sin(exp(x)))
  template <class MPExpr> NLParams CreateExpressionOneArg(MPExpr expr, 
        int coptop) {
      auto ex = GetArgExpression(expr, 0);
      NLParams exp;
      exp.addOp(coptop);
      exp.addMembers(ex);
      return exp;
  }
  /// Create an expression with multiple arguments (e.g. min(x,sin(y),exp(z)))
  template <class MPExpr> NLParams CreateExpressionNArgs(MPExpr expr, int coptop) {
      int n = GetNumArguments(expr);
      NLParams exp;
      exp.addOp(coptop);
      exp.addVar(n); // todo this is clearly not a var
      for (int i = 0; i < n; i++)
          exp.addMembers(GetArgExpression(expr, 1));
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


  ACCEPT_EXPRESSION(LogExpression, Recommended)
      Expr AddExpression(const LogExpression&);
  ACCEPT_EXPRESSION(LogAExpression, Recommended)
      Expr AddExpression(const LogAExpression&);
  ACCEPT_EXPRESSION(ExpExpression, Recommended)
      Expr AddExpression(const ExpExpression&);
  ACCEPT_EXPRESSION(ExpAExpression, Recommended)
      Expr AddExpression(const ExpAExpression&);


  ACCEPT_EXPRESSION(PowConstExpExpression, Recommended)
      Expr AddExpression(const PowConstExpExpression&);


  ACCEPT_EXPRESSION(SinExpression, Recommended)
      Expr AddExpression(const SinExpression&);
  ACCEPT_EXPRESSION(CosExpression, Recommended)
      Expr AddExpression(const CosExpression&);
  ACCEPT_EXPRESSION(TanExpression, Recommended)
      Expr AddExpression(const TanExpression&);

};

} // namespace mp

#endif // COPTMODELAPI_H
