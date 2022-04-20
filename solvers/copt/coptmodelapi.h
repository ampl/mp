#ifndef COPTMODELAPI_H
#define COPTMODELAPI_H

#include <memory>

#include "mp/env.h"
#include "coptcommon.h"
#include "mp/flat/backend_model_api_base.h"
#include "mp/flat/constraints_std.h"

namespace mp {

class CoptModelAPI :
    public CoptCommon, public EnvKeeper,
    public BasicBackendFlatModelAPI
{
  using BaseModelAPI = BasicBackendFlatModelAPI;

public:
  /// Construct
  CoptModelAPI(Env& e) : EnvKeeper(e) { }

  /// Class name
  static const char* GetName() { return "CoptModelAPI"; }

  /// Called before problem input
  void InitProblemModificationPhase();
  /// After
  void FinishProblemModificationPhase();

  void AddVariables(const VarArrayDef& );
  void SetLinearObjective( int iobj, const LinearObjective& lo );
  void SetQuadraticObjective(int iobj, const QuadraticObjective& qo);

  //////////////////////////// GENERAL CONSTRAINTS ////////////////////////////
  USE_BASE_CONSTRAINT_HANDLERS(BaseModelAPI)

  ACCEPT_CONSTRAINT(RangeLinCon, Recommended, CG_Linear)
  void AddConstraint(const RangeLinCon& lc);
  ACCEPT_CONSTRAINT(LinConLE, Recommended, CG_Linear)
  void AddConstraint(const LinConLE& lc);
  ACCEPT_CONSTRAINT(LinConEQ, Recommended, CG_Linear)
  void AddConstraint(const LinConEQ& lc);
  ACCEPT_CONSTRAINT(LinConGE, Recommended, CG_Linear)
  void AddConstraint(const LinConGE& lc);
  ACCEPT_CONSTRAINT(IndicatorConstraintLinLE, AcceptedButNotRecommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinLE& mc);
  ACCEPT_CONSTRAINT(IndicatorConstraintLinEQ, AcceptedButNotRecommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinEQ& mc);

  ACCEPT_CONSTRAINT(QuadraticConstraint, Recommended, CG_Quadratic)
  void AddConstraint(const QuadraticConstraint& qc);

  ACCEPT_CONSTRAINT(SOS1Constraint, Recommended, CG_SOS)
  void AddConstraint(const SOS1Constraint& cc);
  ACCEPT_CONSTRAINT(SOS2Constraint, Recommended, CG_SOS)
  void AddConstraint(const SOS2Constraint& cc);

};

} // namespace mp

#endif // COPTMODELAPI_H
