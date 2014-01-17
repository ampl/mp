/*
 A solver used for testing the C API.

 Copyright (C) 2014 AMPL Optimization Inc

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

#include "solvers/util/solver.h"

#undef getenv

#include <cstdlib>
#include <stdexcept>

namespace ampl {

class TestSolver : public Solver {
 protected:
  void DoSolve(Problem &) {
    const char *fail = std::getenv("ASL_FAIL");
    if (fail && std::strcmp(fail, "1") == 0)
      throw std::runtime_error("epic fail");
  }

  std::string GetOption(const SolverOption &) const { return ""; }
  void SetOption(const SolverOption &, const char * ) {}

 public:
  TestSolver() : Solver("testsolver") {
    set_option_header("Options rock!");
    AddStrOption("opt1", "desc1",
        &TestSolver::GetOption, &TestSolver::SetOption);
    static const OptionValueInfo VALUES[] = {
        {"val1", "valdesc1"},
        {"val2", "valdesc2"},
        {"val3", "valdesc3"}
    };
    AddStrOption("opt2", "desc2",
        &TestSolver::GetOption, &TestSolver::SetOption, VALUES);
    Problem p;
    DoSolve(p);
  }
};

SolverPtr CreateSolver() { return SolverPtr(new TestSolver()); }
}
