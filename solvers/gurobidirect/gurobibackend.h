#ifndef MP_GUROBI_BACKEND_H_
#define MP_GUROBI_BACKEND_H_

#if __clang__
# pragma clang diagnostic push
# pragma clang diagnostic ignored "-Wunused-parameter"
# pragma clang diagnostic ignored "-Wunused-private-field"
#elif _MSC_VER
# pragma warning(push)
# pragma warning(disable: 4244)
#endif

extern "C" {
  #include "gurobi_c.h"
}

#if __clang__
# pragma clang diagnostic pop
#elif _MSC_VER
# pragma warning(pop)
#endif

#include <string>

#include "mp/backend.h"
#include "mp/convert/std_constr.h"

namespace mp {

class GurobiBackend : public BasicBackend<GurobiBackend>
{
  using BaseBackend = BasicBackend<GurobiBackend>;

  //////////////////// [[ The public interface ]] //////////////////////
public:
  GurobiBackend();
  ~GurobiBackend();
  static const char* GetBackendName();

  /// [[ Prototype an incremental interface ]]
  void InitProblemModificationPhase();
  void FinishProblemModificationPhase();

  void AddVariable(Variable var);
  void AddLinearObjective( const LinearObjective& lo );

  //////////////////////////// GENERAL CONSTRAINTS ////////////////////////////
  USE_BASE_CONSTRAINT_HANDLERS(BaseBackend)

  ACCEPT_CONSTRAINT(LinearConstraint, Recommended)
  /// TODO Attributes (lazy/user cut, etc)
  void AddConstraint(const LinearConstraint& lc);
  ACCEPT_CONSTRAINT(MaximumConstraint, AcceptedButNotRecommended)
  void AddConstraint(const MaximumConstraint& mc);
  ACCEPT_CONSTRAINT(MinimumConstraint, AcceptedButNotRecommended)
  void AddConstraint(const MinimumConstraint& mc);
  ACCEPT_CONSTRAINT(AbsConstraint, AcceptedButNotRecommended)
  void AddConstraint(const AbsConstraint& absc);
  ACCEPT_CONSTRAINT(ConjunctionConstraint, AcceptedButNotRecommended)
  void AddConstraint(const ConjunctionConstraint& cc);
  ACCEPT_CONSTRAINT(DisjunctionConstraint, AcceptedButNotRecommended)
  void AddConstraint(const DisjunctionConstraint& mc);
  /// Enabling built-in indicator for infinite bounds,
  /// but not recommended otherwise --- may be slow
  ACCEPT_CONSTRAINT(IndicatorConstraintLinLE, AcceptedButNotRecommended)
  void AddConstraint(const IndicatorConstraintLinLE& mc);

  /// Nonlinear
  ACCEPT_CONSTRAINT(ExpConstraint, Recommended)
  void AddConstraint(const ExpConstraint& cc);
  ACCEPT_CONSTRAINT(ExpAConstraint, Recommended)
  void AddConstraint(const ExpAConstraint& cc);
  ACCEPT_CONSTRAINT(LogConstraint, Recommended)
  void AddConstraint(const LogConstraint& cc);
  ACCEPT_CONSTRAINT(LogAConstraint, Recommended)
  void AddConstraint(const LogAConstraint& cc);
  ACCEPT_CONSTRAINT(PowConstraint, Recommended)
  void AddConstraint(const PowConstraint& cc);
  ACCEPT_CONSTRAINT(SinConstraint, Recommended)
  void AddConstraint(const SinConstraint& cc);
  ACCEPT_CONSTRAINT(CosConstraint, Recommended)
  void AddConstraint(const CosConstraint& cc);
  ACCEPT_CONSTRAINT(TanConstraint, Recommended)
  void AddConstraint(const TanConstraint& cc);


  ////////////////////////////////////////// Model attributes
  bool IsMIP() const;
  bool IsQCP() const;

  int NumberOfConstraints() const;
  int NumberOfVariables() const;
  int NumberOfObjectives() const;

  void ExportModel(const std::string& file);


  //////////////////////////// SOLVING ///////////////////////////////
  void SetInterrupter(mp::Interrupter* inter);
  void DoOptimize();
  std::string ConvertSolutionStatus(
      const mp::Interrupter &interrupter, int &solve_code);

  /// Solution values
  void PrimalSolution(std::vector<double>& x);
  void DualSolution(std::vector<double>& pi);
  double ObjectiveValue() const;

  /// Solution attributes
  double NodeCount() const;
  double Niterations() const;


  //////////////////// [[ Implementation details ]] //////////////////////
  ///////////////////////////////////////////////////////////////////////////////
private:
  GRBenv   *env   = NULL;
  GRBmodel *model = NULL;

public:
  void InitBackend();
  void CloseBackend();

  static double Infinity() { return GRB_INFINITY; }
  static double MinusInfinity() { return -GRB_INFINITY; }

  int GetGrbIntAttribute(const char* attr_id) const;
  double GetGrbDblAttribute(const char* attr_id) const;


protected:
  void InitOptions(); //////////////////////////// OPTIONS ////////////////

private:
  /// These options are stored in the class
  struct Options {
    std::string exportFile_;
  };
  Options storedOptions_;

  /// These options are passed to the solver
  SolverOptionAccessor<GurobiBackend, int, const char*> slvOptInt_ = *this;
  SolverOptionAccessor<GurobiBackend, double, const char*> slvOptDouble_ = *this;
  SolverOptionAccessor<GurobiBackend, std::string, const char*> slvOptString_ = *this;

public:
  void GetSolverOption(const char* key, int& value) const;
  void SetSolverOption(const char* key, int value);
  void GetSolverOption(const char* key, double& value) const;
  void SetSolverOption(const char* key, double value);
  void GetSolverOption(const char* key, std::string& value) const;
  void SetSolverOption(const char* key, const std::string& value);

};

} // namespace mp

#endif  // MP_GUROBI_BACKEND_H_
