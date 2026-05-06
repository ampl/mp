#ifndef KNITROMPMODELAPI_H
#define KNITROMPMODELAPI_H

#include <optional>
#include <map>
#include <unordered_map>

#include "mp/env.h"
#include "knitrompcommon.h"
#include "knitromputils.h"

#include "mp/flat/nl_expr/model_api_base.h"

#define NDEBUG // To disable is_nan check in cppad, that made initializing independent variables a nightmare
#include "cppad/cppad.hpp"


namespace mp {

    // Helper to build CppAD expression from AST
    CppAD::AD<double> buildAD(const ExpressionData& d, const std::map<int, CppAD::AD<double>>& varMap);
    // Create a CppAD tape from this expression tree
    std::pair<CppAD::ADFun<double>, std::vector<int>> createTape(const ExpressionData& d);


    // Data for nonlinear constraints
    class NonlinearConstraintData {
    public:
        double constant=0;
        
        std::shared_ptr<ExpressionData> originalExpr; // Keep the expression if printing is needed
        std::shared_ptr<CppAD::ADFun<double>> tape;  // CppAD tape
        int resultVar=-1;               // result variable index
        std::vector<int> argVars;    // argument variable indices
        int knitroConIndex=-1;          // Knitro constraint index

        // For Hessian computation
        std::vector<std::pair<int, int>> hessianSparsity;  // (row, col) pairs
        // Maps each local sparsity entry to its position in the global Hessian array
        std::vector<size_t> hessianGlobalIndices;

        std::vector<int> linearIndices_;
        std::vector<double> linearCoeffs_;

        // Default constructor
        NonlinearConstraintData() = default;
        
        // Copy constructor (explicitly defaulted - needed by framework)
        NonlinearConstraintData(const NonlinearConstraintData&) = default;
        
        // Copy assignment (explicitly defaulted - needed by framework)
        NonlinearConstraintData& operator=(const NonlinearConstraintData&) = default;
        
        // Move constructor
        NonlinearConstraintData(NonlinearConstraintData&& other) noexcept = default;
        
        // Move assignment
        NonlinearConstraintData& operator=(NonlinearConstraintData&& other) noexcept = default;

        void addLinear(int index, double coeff) {
            linearCoeffs_.push_back(coeff);
            linearIndices_.push_back(index);
        }
    };


class KnitrompModelAPI :
    public KnitrompCommon, public EnvKeeper,
    public BasicExprModelAPI<KnitrompModelAPI, ExpressionData >
{
    using BaseModelAPI = BasicExprModelAPI<KnitrompModelAPI, ExpressionData >;


public:
  /// Construct
  KnitrompModelAPI(Env& e) : EnvKeeper(e) { 
     
  }

  /// Class name
  static const char* GetTypeName() { return "KnitrompModelAPI"; }

  /// If any driver options added from here
  void InitCustomOptions() { }

  /// Called before problem modification.
  /// @param fmi: current problem information.
  /// @note this is called before each phase of model modification
  ///   which can happen during iterative solving.
  void InitProblemModificationPhase(const FlatModelInfo* fmi);
  /// After
  void FinishProblemModificationPhase();

  void AddVariables(const VarArrayDef& );

  void SetLinearObjective( int iobj, const LinearObjective& lo );
  /// Whether accepting quadratic objectives:
  /// 0 - no, 1 - convex, 2 - nonconvex
  static int AcceptsQuadObj() { return 2; }
  /// Implement setting (also changing) a quadratic objective
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
  /// -W Complementarity
  /// - Logical, counting, piecewise-linear constraints.
  /// See \a constr_std.h and other drivers.

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
  static constexpr bool AcceptsNonconvexQC() { return true; }

  /// If cvt:prod=7 (and not 5) default.
  /// Recommendation to return the opposite value as
  /// AcceptsNonconvexQC().
  static constexpr bool WantLogicalizedProd2Bin()
  { return !AcceptsNonconvexQC(); }

  /// If using quadratics,
  /// QuadCon(LE/EQ/GE) should have 'Recommended'
  /// and have an implementation.
  ACCEPT_CONSTRAINT(QuadConLE, Recommended, CG_Quadratic)
  void AddConstraint(const QuadConLE& qc);
  ACCEPT_CONSTRAINT(QuadConEQ, Recommended, CG_Quadratic)
  void AddConstraint(const QuadConEQ& qc);
  ACCEPT_CONSTRAINT(QuadConGE, Recommended, CG_Quadratic)
  void AddConstraint(const QuadConGE& qc);

  ACCEPT_CONSTRAINT(SinConstraint, Recommended, CG_Nonlinear)
  void AddConstraint(const SinConstraint& c);
  ACCEPT_CONSTRAINT(CosConstraint, Recommended, CG_Nonlinear)
  void AddConstraint(const CosConstraint& c);
  ACCEPT_CONSTRAINT(TanConstraint, Recommended, CG_Nonlinear)
  void AddConstraint(const TanConstraint& c);
  ACCEPT_CONSTRAINT(AsinConstraint, Recommended, CG_Nonlinear)
  void AddConstraint(const AsinConstraint& c);
  ACCEPT_CONSTRAINT(AcosConstraint, Recommended, CG_Nonlinear)
  void AddConstraint(const AcosConstraint& c);
  ACCEPT_CONSTRAINT(AtanConstraint, Recommended, CG_Nonlinear)
  void AddConstraint(const AtanConstraint& c);
  ACCEPT_CONSTRAINT(SinhConstraint, Recommended, CG_Nonlinear)
  void AddConstraint(const SinhConstraint& c);
  ACCEPT_CONSTRAINT(CoshConstraint, Recommended, CG_Nonlinear)
  void AddConstraint(const CoshConstraint& c);
  ACCEPT_CONSTRAINT(TanhConstraint, Recommended, CG_Nonlinear)
  void AddConstraint(const TanhConstraint& c);
  ACCEPT_CONSTRAINT(AsinhConstraint, Recommended, CG_Nonlinear)
  void AddConstraint(const AsinhConstraint& c);
  ACCEPT_CONSTRAINT(AcoshConstraint, Recommended, CG_Nonlinear)
  void AddConstraint(const AcoshConstraint& c);
  ACCEPT_CONSTRAINT(AtanhConstraint, Recommended, CG_Nonlinear)
  void AddConstraint(const AtanhConstraint& c);

  ACCEPT_CONSTRAINT(LogConstraint, Recommended, CG_Nonlinear)
  void AddConstraint(const LogConstraint& c);
  ACCEPT_CONSTRAINT(LogAConstraint, Recommended, CG_Nonlinear)
  void AddConstraint(const LogAConstraint& c);
  ACCEPT_CONSTRAINT(PowConstraint, Recommended, CG_Nonlinear)
  void AddConstraint(const PowConstraint& c);
  ACCEPT_CONSTRAINT(ExpConstraint, Recommended, CG_Nonlinear)
  void AddConstraint(const ExpConstraint& c);
  ACCEPT_CONSTRAINT(DivConstraint, Recommended, CG_Nonlinear)
  void AddConstraint(const DivConstraint& c);

  //////////////////////////// EXPRESSION TREES ////////////////////////////
    /// Handle expression trees: inherit basic API
  USE_BASE_EXPRESSION_HANDLERS(BaseModelAPI)
  ACCEPT_EXPRESSION_INTERFACE(Recommended)


  /// Whether accepts NLObjective
  static int AcceptsNLObj() { return 1; }
  void SetNLObjective(int, const NLObjective&);

 
  /// GetVarExpression(\a i): expression representing variable 0<=i<n_var.
  /// Only called for 'nonlinear' variables.
  Expr GetVarExpression(int i) { return ExpressionData::MakeVarExpr(i); }

  /// GetZeroExpr(): constant 0.0 expression.
  /// Can be used to represent empty expression in an NLConstraint.
  Expr GetZeroExpression() { return ExpressionData::MakeConstantExpr(0.0); }

  template <class MPExpr> ExpressionData AddUnaryExpression(ExpressionData::OpType, const MPExpr& expr);


  ACCEPT_CONSTRAINT(NLConstraint, Recommended, CG_General)
  void AddConstraint(const NLConstraint& nl);

  ACCEPT_CONSTRAINT(NLAssignGE, Recommended, CG_General)
  void AddConstraint(const NLAssignGE& nl);
  ACCEPT_CONSTRAINT(NLAssignLE, Recommended, CG_General)
  void AddConstraint(const NLAssignLE& nl);
  ACCEPT_CONSTRAINT(NLAssignEQ, Recommended, CG_General)
  void AddConstraint(const NLAssignEQ& nl);

  ACCEPT_EXPRESSION(NLAffineExpression, Recommended);
  Expr AddExpression(const NLAffineExpression& le);
  
  ACCEPT_EXPRESSION(NLQuadExpression, Recommended);
  Expr AddExpression(const NLQuadExpression& le);
  
  ACCEPT_EXPRESSION(SinExpression, Recommended)
  Expr AddExpression(const SinExpression&);

  ACCEPT_EXPRESSION(CosExpression, Recommended)
      Expr AddExpression(const CosExpression&);

  ACCEPT_EXPRESSION(TanExpression, Recommended)
      Expr AddExpression(const TanExpression&);

  ACCEPT_EXPRESSION(SinhExpression, Recommended)
      Expr AddExpression(const SinhExpression&);

  ACCEPT_EXPRESSION(CoshExpression, Recommended)
      Expr AddExpression(const CoshExpression&);

  ACCEPT_EXPRESSION(TanhExpression, Recommended)
      Expr AddExpression(const TanhExpression&);

  ACCEPT_EXPRESSION(AsinExpression, Recommended)
      Expr AddExpression(const AsinExpression&);

  ACCEPT_EXPRESSION(AcosExpression, Recommended)
      Expr AddExpression(const AcosExpression&);

  ACCEPT_EXPRESSION(AtanExpression, Recommended)
      Expr AddExpression(const AtanExpression&);

  ACCEPT_EXPRESSION(AsinhExpression, Recommended)
      Expr AddExpression(const AsinhExpression&);

  ACCEPT_EXPRESSION(AcoshExpression, Recommended)
      Expr AddExpression(const AcoshExpression&);

  ACCEPT_EXPRESSION(AtanhExpression, Recommended)
      Expr AddExpression(const AtanhExpression&);


  ACCEPT_EXPRESSION(LogExpression, Recommended)
      Expr AddExpression(const LogExpression&);

  ACCEPT_EXPRESSION(LogAExpression, Recommended)
      Expr AddExpression(const LogAExpression&);

  ACCEPT_EXPRESSION(DivExpression, Recommended)
      Expr AddExpression(const DivExpression&);

  ACCEPT_EXPRESSION(PowExpression, Recommended)
      Expr AddExpression(const PowExpression&);


  ACCEPT_EXPRESSION(PowConstExpExpression, Recommended)
      Expr AddExpression(const PowConstExpExpression&);

  ACCEPT_EXPRESSION(AbsExpression, Recommended)
	  Expr AddExpression(const AbsExpression&);


  private:

      template<typename LinConType>
      void AddLinearConstraintHelper(const LinConType& lc,
          double lb, double ub);

      template<typename QuadConType>
      void AddQuadraticConstraintHelper(const QuadConType& qc,
          double lb, double ub);

    // Static callback for evaluating nonlinear constraints
    static int evalNonlinearConstraint(
        KN_context_ptr kc,
        CB_context_ptr cb,
        KN_eval_request_ptr const evalRequest,
        KN_eval_result_ptr const evalResult,
        void* const userParams);

    static int evalHessian(
        KN_context_ptr kc,
        CB_context_ptr cb,
        KN_eval_request_ptr const evalRequest,
        KN_eval_result_ptr const evalResult,
        void* const userParams);


    void markLinearVariables();
    void addUnaryNonLinearConstraint(int resultVar, int argVar, ExpressionData::OpType op, const char* name = nullptr);


    enum ConstraintType {
        EQ, LE, GE
    };

        // Helper to add a constraint of the type
     // resultVar =,<=,>= f(argVars)
    void addAssignNonlinearConstraint(
        ConstraintType type, double rhs, int resultVar, const ExpressionData& exp, const char* name);
   

    // Store the nonlinear constraint data for evaluation.
    // To be called after adding the constraint to Knitro.
    void storeNonLinearData(int conIndex, const ExpressionData& exp,
        int resultVar, const char* name = nullptr);

    void formatModel(fmt::MemoryWriter &w) const;

    std::optional<NonlinearConstraintData> nlObjective_;
    std::vector<NonlinearConstraintData> nlConstraints_;

    // Set of variables that appear in nonlinear constraints, to not to maek them as 'linear' in Knitro
    // Note: marking the vars is not quite stable now
	std::unordered_set<int> nonlinearVars_; 
	// Keep track of which constraint is nonlinear, useful when printing the model 
    // to know if to print from printConstraints_ or from nlConstraints_.
    std::vector<int> nonLinearConstraintsIndex_;
    // Map from Knitro constraint index to nlConstraints_ index useful for printing
	std::unordered_map<int, int> knitroToNLConstraintIndex_; 

    std::unordered_map<int, PrintConstraint> printConstraints_;

    std::map<std::pair<int, int>, size_t> hessSparsityMap_;

};

} // namespace mp

#endif // KNITROMPMODELAPI_H
