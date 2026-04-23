#ifndef KNITROMPMODELAPI_H
#define KNITROMPMODELAPI_H

#include "mp/env.h"
#include "knitrompcommon.h"
#include "mp/flat/model_api_base.h"

#include "cppad/cppad.hpp"


namespace mp {


class KnitrompModelAPI :
    public KnitrompCommon, public EnvKeeper,
    public BasicFlatModelAPI
{

  /// If you don't need expression support, use this typedef
   using BaseModelAPI = BasicFlatModelAPI;

  
public:
  /// Construct
  KnitrompModelAPI(Env& e) : EnvKeeper(e) { }

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

  private:
    // Helper function to fetch variable names from the solver API stub,
    // used to print the model variable names
    std::function<std::string_view(int)> GetVarName;
    int nc = 0; // number of constraints used during model buildup

    // CppAD tape data for nonlinear constraints
    struct NonlinearConstraintData {
        CppAD::ADFun<double> tape;  // CppAD tape
        int resultVar;               // result variable index
        std::vector<int> argVars;    // argument variable indices
        int knitroConIndex;          // Knitro constraint index
    };
    std::vector<NonlinearConstraintData> nlConstraints_;
    std::map<int, int> knitroConIndexToNLConIndex_;  // Map Knitro con index to nlConstraints_ index

    // Static callback for evaluating nonlinear constraints
    static int evalNonlinearConstraint(
        KN_context_ptr kc,
        CB_context_ptr cb,
        KN_eval_request_ptr const evalRequest,
        KN_eval_result_ptr const evalResult,
        void* const userParams);

    // Helper to add nonlinear constraint with CppAD tape
    template<typename UnaryOp>
    void AddUnaryNonLinearConstraint(int resultVar, int argVar, UnaryOp op);

    void AddNonlinearConstraint(
        int resultVar,
        const std::vector<int>& argVars,
        CppAD::ADFun<double>&& tape);

};

} // namespace mp

#endif // KNITROMPMODELAPI_H
