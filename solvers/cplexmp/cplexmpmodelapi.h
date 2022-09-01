#ifndef CPLEXMODELAPI_H
#define CPLEXMODELAPI_H

#include <memory>

#include "mp/env.h"
#include "cplexmpcommon.h"
#include "mp/flat/model_api_base.h"
#include "mp/flat/constr_std.h"

namespace mp {

class CplexModelAPI :
    public CplexCommon, public EnvKeeper,
    public BasicFlatModelAPI
{
  using BaseModelAPI = BasicFlatModelAPI;

public:
  /// Construct
  CplexModelAPI(Env& e) : EnvKeeper(e) { }

  /// Class name
  static const char* GetTypeName() { return "CplexModelAPI"; }

  /// Called before problem input
  void InitProblemModificationPhase(const FlatModelInfo*);
  /// After
  void FinishProblemModificationPhase();

  void AddVariables(const VarArrayDef& );
  void SetLinearObjective( int iobj, const LinearObjective& lo );

  //////////////////////////// GENERAL CONSTRAINTS ////////////////////////////
  USE_BASE_CONSTRAINT_HANDLERS(BaseModelAPI)

  /// LinConRange is optional
  ACCEPT_CONSTRAINT(LinConRange, Recommended, CG_Linear)
  void AddConstraint(const LinConRange& lc);

  /// Every Backend should accept LinCon(LE/EQ/GE),
  /// or add a conversion rule in a derived FlatConverter
  ACCEPT_CONSTRAINT(LinConLE, Recommended, CG_Linear)
  void AddConstraint(const LinConLE& lc);
  ACCEPT_CONSTRAINT(LinConEQ, Recommended, CG_Linear)
  void AddConstraint(const LinConEQ& lc);
  ACCEPT_CONSTRAINT(LinConGE, Recommended, CG_Linear)
  void AddConstraint(const LinConGE& lc);

  /// Discrete generals constraints

  /// Enabling built-in indicator for infinite bounds,
  /// but not recommended otherwise --- may be slow as of CPLEX 12.10
  ACCEPT_CONSTRAINT(IndicatorConstraintLinLE, AcceptedButNotRecommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinLE& mc);
  ACCEPT_CONSTRAINT(IndicatorConstraintLinEQ, AcceptedButNotRecommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinEQ& mc);
  ACCEPT_CONSTRAINT(IndicatorConstraintLinGE, AcceptedButNotRecommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinGE& mc);

  ACCEPT_CONSTRAINT(PLConstraint, Recommended, CG_General)
  void AddConstraint(const PLConstraint& cc);


};

} // namespace mp

#endif // CPLEXMODELAPI_H
