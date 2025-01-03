#ifndef CUOPTMODELAPI_H
#define CUOPTMODELAPI_H

#include "mp/env.h"
#include "cuoptcommon.h"
#include "mp/flat/model_api_base.h"

namespace mp {


/// CuoptModelAPI.
class CuoptModelAPI :
    public CuoptCommon, public EnvKeeper,
    public BasicFlatModelAPI
{
  using BaseModelAPI = BasicFlatModelAPI;

public:
  /// Construct
  CuoptModelAPI(Env& e) : EnvKeeper(e) { }

  /// Class name
  static const char* GetTypeName() { return "CuoptModelAPI"; }

  /// If any driver options added from here
  void InitCustomOptions() { }

  /// Called before problem modification.
  /// @param fmi: current problem information.
  /// @note this is called before each phase of model modification
  ///   which can happen during iterative solving.
  void InitProblemModificationPhase(const FlatModelInfo* fmi);
  /// After
  void FinishProblemModificationPhase();

  /// Implement adding variables
  void AddVariables(const VarArrayDef& );
  /// Implement setting (also changing) a linear (part of the) objective
  void SetLinearObjective( int iobj, const LinearObjective& lo );
  /// Whether accepting quadratic objectives:
  /// 0 - no, 1 - convex, 2 - nonconvex
  static int AcceptsQuadObj() { return 0; }
  /// Implement setting (also changing) a quadratic objective
  void SetQuadraticObjective(int iobj, const QuadraticObjective& qo);

  //////////////////////////// GENERAL CONSTRAINTS ////////////////////////////
  USE_BASE_CONSTRAINT_HANDLERS(BaseModelAPI)

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

private:
  int nvars_;
};

} // namespace mp

#endif // CUOPTMODELAPI_H
