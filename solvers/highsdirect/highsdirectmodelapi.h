#ifndef HIGHSMODELAPI_H
#define HIGHSMODELAPI_H

#include <memory>

#include "mp/env.h"
#include "highsdirectcommon.h"
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

  void AddVariables(const VarArrayDef& );
  void SetLinearObjective( int iobj, const LinearObjective& lo );
  void SetQuadraticObjective(int iobj, const QuadraticObjective& qo);

  //////////////////////////// GENERAL CONSTRAINTS ////////////////////////////
  USE_BASE_CONSTRAINT_HANDLERS(BaseModelAPI)


  struct AccConstraints{
      /* This is to accumulate the constraints in a format suitable for
      Highs_addRows(...). Adding them one by one was killing performance
      unacceptably. */
      std::vector<double> lb, ub, coeffs;
      std::vector<HighsInt> starts, indices;

      void add(const LinConRange& lc)
      {
        lb.push_back(lc.lb());
        ub.push_back(lc.ub());
        coeffs.insert(coeffs.end(), lc.coefs().begin(), lc.coefs().end());
        indices.insert(indices.end(), lc.vars().begin(), lc.vars().end());
        starts.push_back(indices.size());
      }
      void add(const LinConLE& lc)
      {
        lb.push_back(lc.lb());
        ub.push_back(lc.ub());
        coeffs.insert(coeffs.end(), lc.coefs().begin(), lc.coefs().end());
        indices.insert(indices.end(), lc.vars().begin(), lc.vars().end());
        starts.push_back(indices.size());
      }
      void add(const LinConEQ& lc)
      {
        lb.push_back(lc.lb());
        ub.push_back(lc.ub());
        coeffs.insert(coeffs.end(), lc.coefs().begin(), lc.coefs().end());
        indices.insert(indices.end(), lc.vars().begin(), lc.vars().end());
        starts.push_back(indices.size());
      }
      void add(const LinConGE& lc)
      {
        lb.push_back(lc.lb());
        ub.push_back(lc.ub());
        coeffs.insert(coeffs.end(), lc.coefs().begin(), lc.coefs().end());
        indices.insert(indices.end(), lc.vars().begin(), lc.vars().end());
        starts.push_back(indices.size());
      }
      AccConstraints() {
        starts.push_back(0);
      }
  } AccConstraints;

  ACCEPT_CONSTRAINT(LinConRange, Recommended, CG_Linear)
  void AddConstraint(const LinConRange& lc);
  ACCEPT_CONSTRAINT(LinConLE, Recommended, CG_Linear)
  void AddConstraint(const LinConLE& lc);
  ACCEPT_CONSTRAINT(LinConEQ, Recommended, CG_Linear)
  void AddConstraint(const LinConEQ& lc);
  ACCEPT_CONSTRAINT(LinConGE, Recommended, CG_Linear)
  void AddConstraint(const LinConGE& lc);



};

} // namespace mp

#endif // HIGHSMODELAPI_H
