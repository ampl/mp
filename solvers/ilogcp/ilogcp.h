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

#ifndef SOLVERS_ILOGCP_ILOGCP_H_
#define SOLVERS_ILOGCP_ILOGCP_H_

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

#include "util/clock.h"
#include "util/solver.h"

namespace ampl {

class NLToConcertConverter;

// IlogCP solver.
class IlogCPSolver : private Noncopyable, public Solver {
 private:
  IloEnv env_;
  IloCP cp_;
  IloCplex cplex_;

 public:
  // Integer options.
  enum Option {
    DEBUGEXPR,
    USENUMBEROF,
    SOLUTION_LIMIT,
    MULTIOBJ,
    NUM_OPTIONS
  };

  enum Optimizer { AUTO, CP, CPLEX };

 private:
  Optimizer optimizer_;
  int options_[NUM_OPTIONS];

  std::string GetOptimizer(const char *) const;
  void SetOptimizer(const char *name, const char *value);

  int DoGetIntOption(const char *, Option opt) const { return options_[opt]; }
  void SetBoolOption(const char *name, int value, Option opt);
  void DoSetIntOption(const char *name, int value, Option opt);

  // Returns a double option of the constraint programming optimizer.
  double GetCPDblOption(const char *name, IloCP::NumParam param) const;

  // Sets a double option of the constraint programming optimizer.
  void SetCPDblOption(const char *name, double value, IloCP::NumParam param);

  // Returns an integer option of the CPLEX optimizer.
  int GetCPLEXIntOption(const char *name, int param) const;

  // Sets an integer option of the CPLEX optimizer.
  void SetCPLEXIntOption(const char *name, int value, int param);

  struct Stats {
    steady_clock::time_point time;
    double setup_time;
    double solution_time;
  };

  void SolveWithCP(Problem &p,
      const NLToConcertConverter &converter, Stats &stats);
  void SolveWithCPLEX(Problem &p,
      const NLToConcertConverter &converter, Stats &stats);

 protected:

  std::string GetOptionHeader();
  void DoSolve(Problem &p);

 public:
  IlogCPSolver();
  virtual ~IlogCPSolver();

  IloEnv env() const { return env_; }
  IloCP cp() const { return cp_; }
  IloCplex cplex() const { return cplex_; }

  int GetOption(Option opt) const {
    assert(opt >= 0 && opt < NUM_OPTIONS);
    return options_[opt];
  }

  void use_numberof(bool use = true) { options_[USENUMBEROF] = use; }
};
}

#endif  // SOLVERS_ILOGCP_ILOGCP_H_
