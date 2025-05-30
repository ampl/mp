#ifndef CUOPTLPMODELAPI_H
#define CUOPTLPMODELAPI_H

#include "mp/env.h"
#include "cuoptmpcommon.h"
#include "mp/flat/nl_expr/model_api_base.h"

namespace mp {


/// CuoptmpModelAPI.
class CuoptmpModelAPI :
    public CuoptmpCommon, public EnvKeeper,
    public BasicFlatModelAPI
{

  using BaseModelAPI = BasicFlatModelAPI;


public:
  /// Construct
  CuoptmpModelAPI(Env& e) : EnvKeeper(e) { }

  /// Class name
  static const char* GetTypeName() { return "CuoptmpModelAPI"; }

  /// If any driver options added from here
  void InitCustomOptions() { }

  /// Called before problem modification.
  /// @param fmi: current problem information.
  /// @note this is called before each phase of model modification
  ///   which can happen during iterative solving.
  void InitProblemModificationPhase(const FlatModelInfo* fmi);
  /// After
  void FinishProblemModificationPhase();

  void AddVariables(const VarArrayDef& );
  void SetLinearObjective( int iobj, const LinearObjective& lo );
  /// Whether accepting quadratic objectives:
  /// 0 - no, 1 - convex, 2 - nonconvex
  static int AcceptsQuadObj() { return 0; }
  void SetQuadraticObjective(int iobj, const QuadraticObjective& qo) {}

  //////////////////////////// GENERAL CONSTRAINTS ////////////////////////////
  USE_BASE_CONSTRAINT_HANDLERS(BaseModelAPI)

  ACCEPT_CONSTRAINT(LinConLE, Recommended, CG_Linear)
  void AddConstraint(const LinConLE& lc);
  ACCEPT_CONSTRAINT(LinConEQ, Recommended, CG_Linear)
  void AddConstraint(const LinConEQ& lc);
  ACCEPT_CONSTRAINT(LinConGE, Recommended, CG_Linear)
  void AddConstraint(const LinConGE& lc);

  void cuOptAddConstraint(size_t num_coefficients, const double* coefficients, const int* variables, char sense, double rhs);
};

} // namespace mp

#endif // CUOPTLPMODELAPI_H
