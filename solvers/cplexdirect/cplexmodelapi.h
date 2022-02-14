#ifndef CPLEXMODELAPI_H
#define CPLEXMODELAPI_H

#include <memory>

#include "cplexcommon.h"
#include "mp/model_mgr_std_prob.h"

namespace mp {

/// Create Cplex Model Manager from a Backend
std::unique_ptr<BasicModelManager>
CreateCplexModelMgr(CplexCommon&);

class CplexModelAPI : public CplexCommon
{
  using BaseModelAPI = BasicFlatModelAPI;

public:
  CplexModelAPI();

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
  /// but not recommended otherwise --- may be slow
  ACCEPT_CONSTRAINT(IndicatorConstraintLinLE, AcceptedButNotRecommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinLE& mc);
  ACCEPT_CONSTRAINT(IndicatorConstraintLinEQ, AcceptedButNotRecommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinEQ& mc);


};

} // namespace mp

#endif // CPLEXMODELAPI_H
