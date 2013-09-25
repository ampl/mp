/*
 IBM/ILOG CPLEX solver for AMPL.

 Copyright (C) 2013 AMPL Optimization Inc

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

#ifndef SOLVERS_ILOGCP_CPLEX_H_
#define SOLVERS_ILOGCP_CPLEX_H_

#if __clang__
# pragma clang diagnostic push
# pragma clang diagnostic ignored "-Wunused-parameter"
# pragma clang diagnostic ignored "-Wunused-private-field"
#elif _MSC_VER
# pragma warning(push)
# pragma warning(disable: 4244)
#endif

#include <ilcplex/ilocplex.h>

#if __clang__
# pragma clang diagnostic pop
#elif _MSC_VER
# pragma warning(pop)
#endif

#include "solvers/util/solver.h"

namespace ampl {

class CPLEXSolver :
    private Interruptible, private Noncopyable, public Solver<CPLEXSolver> {
 private:
  IloEnv env_;
  IloCplex cplex_;
  IloCplex::Aborter aborter_;

 public:
  // Boolean options.
  enum Option {
    DEBUGEXPR,
    USENUMBEROF,
    NUM_OPTIONS
  };

 private:
  int options_[NUM_OPTIONS];

  int GetBoolOption(const char *, Option opt) const { return options_[opt]; }
  void SetBoolOption(const char *name, int value, Option opt);

  // Returns an integer option of the CPLEX optimizer.
  int GetCPLEXIntOption(const char *name, int param) const;

  // Sets an integer option of the CPLEX optimizer.
  void SetCPLEXIntOption(const char *name, int value, int param);

  void Interrupt() { aborter_.abort(); }

 protected:

  std::string GetOptionHeader();

 public:
  CPLEXSolver();
  virtual ~CPLEXSolver();

  IloEnv env() const { return env_; }
  IloCplex optimizer() const { return cplex_; }

  int GetOption(Option opt) const {
    assert(opt >= 0 && opt < NUM_OPTIONS);
    return options_[opt];
  }

  void use_numberof(bool use = true) { options_[USENUMBEROF] = use; }

  void Solve(Problem &p);
};
}

#endif  // SOLVERS_ILOGCP_CPLEX_H_
