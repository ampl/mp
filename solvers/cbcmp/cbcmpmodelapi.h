#ifndef CBCMPMODELAPI_H
#define CBCMPMODELAPI_H

#include <memory>

#include "mp/env.h"
#include "cbcmpcommon.h"
#include "mp/flat/model_api_base.h"
#include "mp/flat/constr_std.h"

namespace mp {

class CbcmpModelAPI :
    public CbcmpCommon, public EnvKeeper,
    public BasicFlatModelAPI
{
  using BaseModelAPI = BasicFlatModelAPI;

public:
  /// Constructor
  CbcmpModelAPI(Env& e) : EnvKeeper(e) { }

  /// Class name
  static const char* GetTypeName() { return "CbcmpModelAPI"; }

  /// Called before problem input.
  void InitProblemModificationPhase(const FlatModelInfo*);
  /// After problem input.
  void FinishProblemModificationPhase();

  void AddVariables(const VarArrayDef& );
  void SetLinearObjective( int iobj, const LinearObjective& lo );
  /// Whether accepting quadratic objectives:
  /// 0 - no, 1 - convex, 2 - nonconvex
  static int AcceptsQuadObj() { return 0; }


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
  ACCEPT_CONSTRAINT(LinConRange, NotAccepted, CG_Linear)

  /// LinCon(LE/EQ/GE) should have 'Recommended' for all backends
  /// and have an implementation,
  /// or a conversion rule is needed in a derived FlatConverter
  ACCEPT_CONSTRAINT(LinConLE, Recommended, CG_Linear)
  void AddConstraint(const LinConLE& lc);
  ACCEPT_CONSTRAINT(LinConEQ, Recommended, CG_Linear)
  void AddConstraint(const LinConEQ& lc);
  ACCEPT_CONSTRAINT(LinConGE, Recommended, CG_Linear)
  void AddConstraint(const LinConGE& lc);
  /// SOS constraints can be used by AMPL for redefinition of
  /// piecewise-linear expressions.
  /// Set ``option pl_linearize 0;`` in AMPL if the solver
  /// supports PL natively.
  ACCEPT_CONSTRAINT(SOS1Constraint, Recommended, CG_SOS)
  void AddConstraint(const SOS1Constraint& cc);
  ACCEPT_CONSTRAINT(SOS2Constraint, Recommended, CG_SOS)
  void AddConstraint(const SOS2Constraint& cc);

};

} // namespace mp

#endif // CBCMPMODELAPI_H
