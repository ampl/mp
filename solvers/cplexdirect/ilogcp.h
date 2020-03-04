/*
 IBM/ILOG CP solver for AMPL.

 Copyright (C) 2012 AMPL Optimization Inc

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization Inc disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Victor Zverovich
 */

#ifndef MP_SOLVERS_ILOGCP_H_
#define MP_SOLVERS_ILOGCP_H_

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

#include <ilcp/cp.h>
#include <ilcplex/ilocplex.h>

#if __clang__
# pragma clang diagnostic pop
#elif _MSC_VER
# pragma warning(pop)
#endif

#include <string.h> /* This and -fpermissive seem to be needed for MacOSX, */
                    /* at least with g++ 4.6.  Otherwise there are errors */
                    /* with iloconcert/iloenv.h . */
#include <limits.h> /* Needed for g++ -m32 on MacOSX. */
#include <string>

#include "mp/clock.h"
#include "mp/problem.h"
#include "mp/solver.h"

#include "mp/backend.h"

namespace mp {

class MPToConcertConverter;

template <typename T>
struct ParamTraits;

template <>
struct ParamTraits<int> {
  typedef IloCP::IntParam Type;
};

template <>
struct ParamTraits<double> {
  typedef IloCP::NumParam Type;
};

// IlogCP solver.
class IlogCPSolver : public SolverImpl<Problem>,
    public BasicBackend<IlogCPSolver>
{
 private:
  IloEnv env_;
  IloCP cp_;
  IloCplex cplex_;

  FMT_DISALLOW_COPY_AND_ASSIGN(IlogCPSolver);

 public:
  // Integer options.
  enum Option {
    DEBUGEXPR,
    USENUMBEROF,
    SOLUTION_LIMIT,
    NUM_OPTIONS
  };

  enum Optimizer { AUTO, CP, CPLEX };

 private:
  Optimizer optimizer_;
  Optimizer optimizer;
  int options_[NUM_OPTIONS];

  unsigned converter_flags_ = 0;
  std::unique_ptr<MPToConcertConverter> converter_;


  enum FileKind {
    DUMP_FILE,
    EXPORT_FILE,
    NUM_FILES
  };

  std::string filenames_[NUM_FILES];

  std::string GetOptimizer(const SolverOption &) const;
  void SetOptimizer(const SolverOption &opt, fmt::StringRef value);

  int DoGetIntOption(const SolverOption &, Option id) const {
    return options_[id];
  }
  void SetBoolOption(const SolverOption &opt, int value, Option id);
  void DoSetIntOption(const SolverOption &opt, int value, Option id);

  // Returns an option of the constraint programming optimizer.
  template <typename T>
  T GetCPOption(const SolverOption &opt,
                typename ParamTraits<T>::Type param) const;

  // Sets an option of the constraint programming optimizer.
  template <typename T>
  void SetCPOption(const SolverOption &opt, T value,
                   typename ParamTraits<T>::Type param);

  std::string GetFile(const SolverOption &, FileKind kind) const {
    assert(kind < NUM_FILES);
    return filenames_[kind];
  }
  void SetFile(const SolverOption &, fmt::StringRef filename, FileKind kind) {
    assert(kind < NUM_FILES);
    filenames_[kind] = filename.to_string();
  }

  // Returns an integer option of the CPLEX optimizer.
  int GetCPLEXIntOption(const SolverOption &opt, int param) const;

  // Sets an integer option of the CPLEX optimizer.
  void SetCPLEXIntOption(const SolverOption &opt, int value, int param);

  struct Stats {
    steady_clock::time_point time;
    double setup_time;
    double solution_time;
  };
  Stats stats;

  void SolveWithCP(Problem &p, const MPToConcertConverter &converter_,
                   Stats &stats, SolutionHandler &sh);
  void SolveWithCPLEX(Problem &p, const MPToConcertConverter &converter_,
                      Stats &stats, SolutionHandler &sh);

 public:
  IlogCPSolver();
  virtual ~IlogCPSolver();

  IloEnv env() const { return env_; }
  IloCP cp() const { return cp_; }
  IloCplex cplex() const { return cplex_; }

  int GetOption(Option id) const {
    assert(id >= 0 && id < NUM_OPTIONS);
    return options_[id];
  }

  void use_numberof(bool use = true) { options_[USENUMBEROF] = use; }

  void Solve(Problem &p, SolutionHandler &sh);

 public:                    // [[ The interface ]]
  void Convert(Problem& p);
  void Resolve(Problem& p, SolutionHandler &sh);

  /// [[ Surface the incremental interface ]]
  void InitProblemModificationPhase(const Problem& p);
  void AddVariables(int n, double* lbs, double* ubs, var::Type* types);
  /// Supporting linear stuff for now
  void AddLinearObjective( obj::Type sense, int nnz,
                           const double* c, const int* v);
  void AddLinearConstraint(int nnz, const double* c, const int* v,
                           double lb, double ub);
  void FinishProblemModificationPhase();
};
}

#endif  // MP_SOLVERS_ILOGCP_H_
