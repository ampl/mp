#ifndef CPLEXMODELAPI_H
#define CPLEXMODELAPI_H

#include <memory>

#include "mp/env.h"
#include "cplexcommon.h"
#include "mp/flat/backend_model_api_base.h"
#include "mp/flat/constraints_std.h"

namespace mp {

class CplexModelAPI :
    public CplexCommon, public EnvKeeper,
    public BasicBackendFlatModelAPI
{
  using BaseModelAPI = BasicBackendFlatModelAPI;

public:
  /// Construct
  CplexModelAPI(Env& e) : EnvKeeper(e) { }

  /// Class name
  static const char* GetName() { return "CplexModelAPI"; }

  /// Called before problem input
  void InitProblemModificationPhase();
  /// After
  void FinishProblemModificationPhase();

  void AddVariables(const VarArrayDef& );
  void SetLinearObjective( int iobj, const LinearObjective& lo );

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
  /// Enabling built-in indicator for infinite bounds,
  /// but not recommended otherwise --- may be slow as of 12.10
  ACCEPT_CONSTRAINT(IndicatorConstraintLinLE, AcceptedButNotRecommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinLE& mc);
  ACCEPT_CONSTRAINT(IndicatorConstraintLinEQ, AcceptedButNotRecommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinEQ& mc);


};

} // namespace mp

#endif // CPLEXMODELAPI_H
