#ifndef GUROBIMODELAPI_H
#define GUROBIMODELAPI_H

#include <memory>

#include "mp/env.h"
#include "mp/flat/model_api_base.h"
#include "mp/flat/constr_std.h"

#include "gurobicommon.h"

namespace mp {

class GurobiModelAPI :
    public GurobiCommon,
    public EnvKeeper,
    public BasicFlatModelAPI {
  using BaseModelAPI = BasicFlatModelAPI;

public:
  /// Model API name
  static const char* GetTypeName();
  /// Reserved
  static const char* GetLongName() { return nullptr; }

  /// Construct
  GurobiModelAPI(Env& e) : EnvKeeper(e) { }

  /// This is called before model is pushed to the Backend
  void InitProblemModificationPhase(const FlatModelInfo*);
  /// Chance to call GRBupdatemodel()
  void FinishProblemModificationPhase();

  /////////////////////////////////////////////////////////////////////////////
  //////////////////////////// MODELING ACCESSORS /////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  //////////////////////////// VARIABLES //////////////////////////////////////
  void AddVariables(const VarArrayDef& );

  //////////////////////////// OBJECTIVES /////////////////////////////////////
  void SetLinearObjective( int iobj, const LinearObjective& lo );
  void SetQuadraticObjective( int iobj, const QuadraticObjective& qo );

  //////////////////////////// GENERAL CONSTRAINTS ////////////////////////////
  USE_BASE_CONSTRAINT_HANDLERS(BaseModelAPI)

  /// Gurobi does not properly support range linear constraints.
  /// Gurobi 9.5: GRBaddrangeconstr() creates a slack variable but
  /// does not account for it in the basis information, etc.

  /// LinCon(LE/EQ/GE) should have 'Recommended' for all backends
  /// and have an implementation.
  ACCEPT_CONSTRAINT(LinConLE, Recommended, CG_Linear)
  void AddConstraint(const LinConLE& lc);
  ACCEPT_CONSTRAINT(LinConEQ, Recommended, CG_Linear)
  void AddConstraint(const LinConEQ& lc);
  ACCEPT_CONSTRAINT(LinConGE, Recommended, CG_Linear)
  void AddConstraint(const LinConGE& lc);

  /// If using quadratics,
  /// QuadCon(LE/EQ/GE) should have 'Recommended'
  /// and have an implementation.
  /// QuadConRange is optional.
  ACCEPT_CONSTRAINT(QuadConLE, Recommended, CG_Quadratic)
  void AddConstraint(const QuadConLE& qc);
  ACCEPT_CONSTRAINT(QuadConEQ, Recommended, CG_Quadratic)
  void AddConstraint(const QuadConEQ& qc);
  ACCEPT_CONSTRAINT(QuadConGE, Recommended, CG_Quadratic)
  void AddConstraint(const QuadConGE& qc);

  /// Discrete general constraints
  ACCEPT_CONSTRAINT(MaxConstraint, Recommended, CG_General)
  void AddConstraint(const MaxConstraint& mc);
  ACCEPT_CONSTRAINT(MinConstraint, Recommended, CG_General)
  void AddConstraint(const MinConstraint& mc);
  ACCEPT_CONSTRAINT(AbsConstraint, Recommended, CG_General)
  void AddConstraint(const AbsConstraint& absc);
  ACCEPT_CONSTRAINT(AndConstraint, Recommended, CG_General)
  void AddConstraint(const AndConstraint& cc);
  ACCEPT_CONSTRAINT(OrConstraint, Recommended, CG_General)
  void AddConstraint(const OrConstraint& mc);
  ACCEPT_CONSTRAINT(IndicatorConstraintLinLE, Recommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinLE& mc);
  ACCEPT_CONSTRAINT(IndicatorConstraintLinEQ, Recommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinEQ& mc);
  ACCEPT_CONSTRAINT(IndicatorConstraintLinGE, Recommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinGE& mc);
  ACCEPT_CONSTRAINT(PLConstraint, Recommended, CG_General)
  void AddConstraint(const PLConstraint& cc);

  /// Gurobi SOS1/2
  ACCEPT_CONSTRAINT(SOS1Constraint, Recommended, CG_SOS)
  void AddConstraint(const SOS1Constraint& cc);
  ACCEPT_CONSTRAINT(SOS2Constraint, Recommended, CG_SOS)
  void AddConstraint(const SOS2Constraint& cc);

  /// Gurobi nonlinear generals
  ACCEPT_CONSTRAINT(ExpConstraint, Recommended, CG_General)
  void AddConstraint(const ExpConstraint& cc);
  ACCEPT_CONSTRAINT(ExpAConstraint, Recommended, CG_General)
  void AddConstraint(const ExpAConstraint& cc);
  ACCEPT_CONSTRAINT(LogConstraint, Recommended, CG_General)
  void AddConstraint(const LogConstraint& cc);
  ACCEPT_CONSTRAINT(LogAConstraint, Recommended, CG_General)
  void AddConstraint(const LogAConstraint& cc);
  ACCEPT_CONSTRAINT(PowConstraint, Recommended, CG_General)
  void AddConstraint(const PowConstraint& cc);
  ACCEPT_CONSTRAINT(SinConstraint, Recommended, CG_General)
  void AddConstraint(const SinConstraint& cc);
  ACCEPT_CONSTRAINT(CosConstraint, Recommended, CG_General) // y = cos(x)
  void AddConstraint(const CosConstraint& cc);  // GRBaddgenconstrCos(x, y);
  ACCEPT_CONSTRAINT(TanConstraint, Recommended, CG_General)
  void AddConstraint(const TanConstraint& cc);

  /// Init GurobiModelAPI driver options
  void InitOptions();


protected:
  /// First objective's sense
  void NoteGurobiMainObjSense(obj::Type s);
  obj::Type GetGurobiMainObjSense() const;


private:
  /// The sense of the main objective
  obj::Type main_obj_sense_;

  /// These options are stored in the class as variables
  /// for direct access
  struct Options {
  } storedOptions_;
};

} // namespace mp

#endif // GUROBIMODELAPI_H
