#ifndef MP_CPLEX_BACKEND_H_
#define MP_CPLEX_BACKEND_H_

#ifdef __APPLE__
#include <limits.h>
#include <string.h>
#endif

#if __clang__
# pragma clang diagnostic push
# pragma clang diagnostic ignored "-Wunused-parameter"
# pragma clang diagnostic ignored "-Wunused-private-field"
#elif _MSC_VER
# pragma warning(push)
# pragma warning(disable: 4244)
#endif

extern "C" {
  #include <ilcplex/cplex.h>
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

class CplexBackend : public BasicBackend<CplexBackend>
{
  using BaseBackend = BasicBackend<CplexBackend>;

  //////////////////// [[ The public interface ]] //////////////////////
public:
  CplexBackend();
  ~CplexBackend();

  /// [[ Prototype the incremental interface ]]
  void InitProblemModificationPhase();
  void FinishProblemModificationPhase();

  void AddVariables(int n, double* lbs, double* ubs, var::Type* types);
  void AddLinearObjective( obj::Type sense, int nnz,
                           const double* c, const int* v);
  void AddLinearConstraint(int nnz, const double* c, const int* v,
                           double lb, double ub);

  //////////////////////////// GENERAL CONSTRAINTS ////////////////////////////
  USE_BASE_CONSTRAINT_HANDLERS(BaseBackend)

  /// Enabling built-in indicator for infinite bounds,
  /// but not recommended otherwise --- may be slow
  ACCEPT_CONSTRAINT(IndicatorConstraintLinLE, AcceptedButNotRecommended)
  void AddConstraint(const IndicatorConstraintLinLE& mc);


  /////////////////////////// Model attributes
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
  CPXENVptr     env = NULL;
  CPXLPptr      lp = NULL;

public:  // public for static polymorphism
  void InitBackend();
  void CloseBackend();

  static double Infinity() { return CPX_INFBOUND; }
  static double MinusInfinity() { return -CPX_INFBOUND; }

protected:
  void InitOptions();  /////////////////////////////// OPTIONS ///////////////////////

private:
  /// These options are stored in the class
  struct Options {
    std::string exportFile_;
  };
  Options storedOptions_;

  /// These options are passed to the solver
  SolverOptionAccessor<CplexBackend, int, int> slvOptInt_ = *this;
  SolverOptionAccessor<CplexBackend, double, int> slvOptDouble_ = *this;
  SolverOptionAccessor<CplexBackend, std::string, int> slvOptString_ = *this;

public:
  void GetSolverOption(int key, int& value) const;
  void SetSolverOption(int key, int value);
  void GetSolverOption(int key, double& value) const;
  void SetSolverOption(int key, double value);
  void GetSolverOption(int key, std::string& value) const;
  void SetSolverOption(int key, const std::string& value);

};

}  // namespace mp

#endif  // MP_CPLEX_BACKEND_H_
