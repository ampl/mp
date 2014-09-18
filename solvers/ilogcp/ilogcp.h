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
#include "asl/aslsolver.h"

namespace mp {

class NLToConcertConverter;

// IlogCP solver.
class IlogCPSolver : public ASLSolver {
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
    MULTIOBJ,
    OBJNO,
    NUM_OPTIONS
  };

  enum Optimizer { AUTO, CP, CPLEX };

 private:
  Optimizer optimizer_;
  int options_[NUM_OPTIONS];

  std::string GetOptimizer(const SolverOption &) const;
  void SetOptimizer(const SolverOption &opt, const char *value);

  int DoGetIntOption(const SolverOption &, Option id) const {
    return options_[id];
  }
  void SetBoolOption(const SolverOption &opt, int value, Option id);
  void DoSetIntOption(const SolverOption &opt, int value, Option id);

  // Returns a double option of the constraint programming optimizer.
  double GetCPDblOption(const SolverOption &opt, IloCP::NumParam param) const;

  // Sets a double option of the constraint programming optimizer.
  void SetCPDblOption(
      const SolverOption &opt, double value, IloCP::NumParam param);

  // Returns an integer option of the CPLEX optimizer.
  int GetCPLEXIntOption(const SolverOption &opt, int param) const;

  // Sets an integer option of the CPLEX optimizer.
  void SetCPLEXIntOption(const SolverOption &opt, int value, int param);

  struct Stats {
    steady_clock::time_point time;
    double setup_time;
    double solution_time;
  };

  void SolveWithCP(Problem &p, const NLToConcertConverter &converter,
                   Stats &stats, SolutionHandler &sh);
  void SolveWithCPLEX(Problem &p, const NLToConcertConverter &converter,
                      Stats &stats, SolutionHandler &sh);

 protected:

  void DoSolve(Problem &p, SolutionHandler &sh);

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
};
}

#endif  // MP_SOLVERS_ILOGCP_H_
