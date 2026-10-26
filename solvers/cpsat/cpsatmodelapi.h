#ifndef CPSATMODELAPI_H
#define CPSATMODELAPI_H

#include <memory>
#include <vector>

#include "mp/env.h"
#include "cpsatcommon.h"
#include "mp/flat/model_api_base.h"
#include "mp/flat/constr_std.h"

namespace mp {

/// Maps the flat model to a CP-SAT model (CpModelBuilder).
///
/// CP-SAT works on integer variables with integer coefficients only.
/// The constraints handled natively are the linear ones plus
/// abs, min, max, and, or, alldiff, count, numberof, if-then-else,
/// and indicators. Everything else is reformulated by MP, if possible.
class CpsatModelAPI :
    public CpsatCommon, public EnvKeeper,
    public BasicFlatModelAPI
{
  using BaseModelAPI = BasicFlatModelAPI;

public:
  /// Construct
  CpsatModelAPI(Env& e) : EnvKeeper(e) { }

  /// Class name
  static const char* GetTypeName() { return "CpsatModelAPI"; }

  /// CP-SAT has no continuous variables:
  /// keep implied-integer auxiliaries integer (cvt:pre:auxint)
  static int AuxIntDefault() { return 1; }

  /// Products are native: keep a product of sums as one product
  /// of two auxiliary variables (cvt:multoutcard)
  static double MultOutCardDefault() { return 1.0; }

  /// Called before problem input.
  void InitProblemModificationPhase(const FlatModelInfo*);
  /// After
  void FinishProblemModificationPhase();

  void AddVariables(const VarArrayDef& );
  void SetLinearObjective( int iobj, const LinearObjective& lo );
  void SetQuadraticObjective(int iobj, const QuadraticObjective& qo);

  //////////////////////////// GENERAL CONSTRAINTS ////////////////////////////
  USE_BASE_CONSTRAINT_HANDLERS(BaseModelAPI)

  /// Linear constraints
  ACCEPT_CONSTRAINT(LinConRange, Recommended, CG_Linear)
  void AddConstraint(const LinConRange& lc);
  ACCEPT_CONSTRAINT(LinConLE, Recommended, CG_Linear)
  void AddConstraint(const LinConLE& lc);
  ACCEPT_CONSTRAINT(LinConEQ, Recommended, CG_Linear)
  void AddConstraint(const LinConEQ& lc);
  ACCEPT_CONSTRAINT(LinConGE, Recommended, CG_Linear)
  void AddConstraint(const LinConGE& lc);

  /// Arithmetic: r = max/min/abs(...)
  ACCEPT_CONSTRAINT(MaxConstraint, Recommended, CG_General)
  void AddConstraint(const MaxConstraint& mc);
  ACCEPT_CONSTRAINT(MinConstraint, Recommended, CG_General)
  void AddConstraint(const MinConstraint& mc);
  ACCEPT_CONSTRAINT(AbsConstraint, Recommended, CG_General)
  void AddConstraint(const AbsConstraint& ac);
  /// Products: r = x*y (from quadratic terms, see acc:mul),
  /// r = x^k for integer k >= 0
  ACCEPT_CONSTRAINT(MulConstraint, Recommended, CG_General)
  void AddConstraint(const MulConstraint& mc);
  ACCEPT_CONSTRAINT(PowConstExpConstraint, Recommended, CG_General)
  void AddConstraint(const PowConstExpConstraint& pc);
  /// r = x div y = trunc(x/y), r = x mod y (sign of x), as in AMPL
  ACCEPT_CONSTRAINT(IntDivConstraint, Recommended, CG_General)
  void AddConstraint(const IntDivConstraint& dc);
  ACCEPT_CONSTRAINT(ModConstraint, Recommended, CG_General)
  void AddConstraint(const ModConstraint& mc);

  /// Logic: r = and/or(...), indicators
  ACCEPT_CONSTRAINT(AndConstraint, Recommended, CG_General)
  void AddConstraint(const AndConstraint& cc);
  ACCEPT_CONSTRAINT(OrConstraint, Recommended, CG_General)
  void AddConstraint(const OrConstraint& cc);
  ACCEPT_CONSTRAINT(IndicatorConstraintLinLE, Recommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinLE& ic);
  ACCEPT_CONSTRAINT(IndicatorConstraintLinEQ, Recommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinEQ& ic);
  ACCEPT_CONSTRAINT(IndicatorConstraintLinGE, Recommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinGE& ic);

  /// Constraint programming:
  /// alldiff, count, numberof, if-then-else
  ACCEPT_CONSTRAINT(AllDiffConstraint, Recommended, CG_General)
  void AddConstraint(const AllDiffConstraint& ac);
  ACCEPT_CONSTRAINT(CountConstraint, Recommended, CG_General)
  void AddConstraint(const CountConstraint& cc);
  ACCEPT_CONSTRAINT(NumberofConstConstraint, Recommended, CG_General)
  void AddConstraint(const NumberofConstConstraint& nc);
  ACCEPT_CONSTRAINT(NumberofVarConstraint, Recommended, CG_General)
  void AddConstraint(const NumberofVarConstraint& nc);
  ACCEPT_CONSTRAINT(IfThenConstraint, Recommended, CG_General)
  void AddConstraint(const IfThenConstraint& ic);
  ACCEPT_CONSTRAINT(ImplicationConstraint, Recommended, CG_General)
  void AddConstraint(const ImplicationConstraint& ic);

private:
  using IntVar = operations_research::sat::IntVar;
  using BoolVar = operations_research::sat::BoolVar;
  using LinearExpr = operations_research::sat::LinearExpr;
  using Domain = operations_research::Domain;

  IntVar Var(int i) const { return data()->vars.at(i); }
  /// Variable used as a literal; must be in [0,1]
  BoolVar Lit(int i) const;
  std::vector<IntVar> Vars(const std::vector<int>& idx) const;
  /// Domain from bounds for an integer-valued linear expression
  static Domain MakeDomain(double lb, double ub);
  bool IsFixedTo(int var, int64_t value) const;
  /// The divisor y, or a copy with 0 removed from its domain
  /// (CP-SAT requires that for division)
  IntVar NonzeroDivisor(int y);
  /// Add r <==> and(lits)
  void AddAnd(int r, const std::vector<BoolVar>& lits);
  template <class Con>
  void AddIndicator(const Con& ic);
};

} // namespace mp

#endif // CPSATMODELAPI_H
