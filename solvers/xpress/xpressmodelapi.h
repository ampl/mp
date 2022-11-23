#ifndef XPRESSMPMODELAPI_H
#define XPRESSMPMODELAPI_H

#include <memory>

#include "mp/env.h"
#include "xpresscommon.h"
#include "mp/flat/model_api_base.h"
#include "mp/flat/constr_std.h"

namespace mp {

class XpressmpModelAPI :
    public XpressmpCommon, public EnvKeeper,
    public BasicFlatModelAPI
{
  using BaseModelAPI = BasicFlatModelAPI;

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
  template <class Args, class Params, class NumOrLogic, class Id> void addGenCon(
    const CustomFunctionalConstraint<Args, Params, NumOrLogic, Id>& c, int xpressConType);
  ACCEPT_CONSTRAINT(AbsConstraint, Recommended, CG_General)
  void AddConstraint(const AbsConstraint& ac);
  ACCEPT_CONSTRAINT(MaxConstraint, Recommended, CG_General)
    void AddConstraint(const MaxConstraint& ac);
  ACCEPT_CONSTRAINT(MinConstraint, Recommended, CG_General)
    void AddConstraint(const MinConstraint& ac);
  ACCEPT_CONSTRAINT(OrConstraint, Recommended, CG_General)
    void AddConstraint(const OrConstraint& ac);
  ACCEPT_CONSTRAINT(AndConstraint, Recommended, CG_General)
    void AddConstraint(const AndConstraint& ac);
};

} // namespace mp

#endif // XPRESSMPMODELAPI_H
