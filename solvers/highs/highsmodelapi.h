#ifndef HIGHSMODELAPI_H
#define HIGHSMODELAPI_H

#include <memory>

#include "mp/env.h"
#include "highscommon.h"
#include "mp/flat/model_api_base.h"
#include "mp/flat/constraints_std.h"

namespace mp {

class HighsModelAPI :
    public HighsCommon, public EnvKeeper,
    public BasicFlatModelAPI
{
  using BaseModelAPI = BasicFlatModelAPI;

public:
  /// Construct
  HighsModelAPI(Env& e) : EnvKeeper(e) { }

  /// Class name
  static const char* GetTypeName() { return "HighsModelAPI"; }

  /// Called before problem input
  void InitProblemModificationPhase();
  /// After
  void FinishProblemModificationPhase();

  /// TODO Implement the following functions using the solver's API
  void AddVariables(const VarArrayDef& );
  void SetLinearObjective( int iobj, const LinearObjective& lo );
  void SetQuadraticObjective(int iobj, const QuadraticObjective& qo);

  //////////////////////////// GENERAL CONSTRAINTS ////////////////////////////
  USE_BASE_CONSTRAINT_HANDLERS(BaseModelAPI)

  /// TODO For each suppoted constraint type, add the ACCEPT_CONSTRAINT macro
  /// and the relative AddConstraint function
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



};

} // namespace mp

#endif // HIGHSMODELAPI_H
