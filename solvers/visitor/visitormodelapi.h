#ifndef VISITORMODELAPI_H
#define VISITORMODELAPI_H

#include <memory>

#include "mp/env.h"
#include "visitorcommon.h"
#include "mp/flat/model_api_base.h"
#include "mp/flat/constr_std.h"

namespace mp {

class VisitorModelAPI :
    public VisitorCommon, public EnvKeeper,
    public BasicFlatModelAPI
{
  using BaseModelAPI = BasicFlatModelAPI;

public:
  /// Construct
  VisitorModelAPI(Env& e) : EnvKeeper(e) { }

  /// Class name
  static const char* GetTypeName() { return "VisitorModelAPI"; }

  /// If any driver options added from here
  void InitCustomOptions() { }

  /// Called before problem input.
  /// Model info can be used to preallocate memory.
  void InitProblemModificationPhase(const FlatModelInfo*);
  /// After
  void FinishProblemModificationPhase();

  /// TODO Implement the following functions using the solver's API
  void AddVariables(const VarArrayDef& );
  void SetLinearObjective( int iobj, const LinearObjective& lo );
  /// Whether accepting quadratic objectives:
  /// 0 - no, 1 - convex, 2 - nonconvex
  static int AcceptsQuadObj() { return 0; }
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
  ACCEPT_CONSTRAINT(IndicatorConstraintLinLE, AcceptedButNotRecommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinLE& mc);
  ACCEPT_CONSTRAINT(IndicatorConstraintLinEQ, AcceptedButNotRecommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinEQ& mc);
  ACCEPT_CONSTRAINT(IndicatorConstraintLinGE, AcceptedButNotRecommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinGE& mc);

  /// SOS constraints can be used by AMPL for redefinition of
  /// piecewise-linear expressions.
  /// Set ``option pl_linearize 0;`` in AMPL if the solver
  /// supports PL natively.
  ACCEPT_CONSTRAINT(SOS1Constraint, Recommended, CG_SOS)
  void AddConstraint(const SOS1Constraint& cc);
  ACCEPT_CONSTRAINT(SOS2Constraint, Recommended, CG_SOS)
  void AddConstraint(const SOS2Constraint& cc);

  ACCEPT_CONSTRAINT(PLConstraint, Recommended, CG_General)
  void AddConstraint(const PLConstraint& cc);



  // Non linear constraints
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
  ACCEPT_CONSTRAINT(PowConstraint, Recommended, CG_General)
    void AddConstraint(const PowConstraint& cc);
  ACCEPT_CONSTRAINT(SinConstraint, Recommended, CG_General)
    void AddConstraint(const SinConstraint& cc);
  ACCEPT_CONSTRAINT(CosConstraint, Recommended, CG_General)
    void AddConstraint(const CosConstraint& cc);  
  ACCEPT_CONSTRAINT(TanConstraint, Recommended, CG_General)
    void AddConstraint(const TanConstraint& cc);

};

} // namespace mp

#endif // VISITORMODELAPI_H
