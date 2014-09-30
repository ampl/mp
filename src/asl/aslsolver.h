/*
 ASL solver.

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

#ifndef MP_ASL_SOLVER_H_
#define MP_ASL_SOLVER_H_

#include "mp/solver.h"
#include "asl/aslbuilder.h"

namespace mp {

class ASLSolver : public SolverImpl<internal::ASLBuilder> {
 private:
  void RegisterSuffixes(ASL *asl);

 public:
  ASLSolver(fmt::StringRef name, fmt::StringRef long_name = 0,
            long date = 0, int flags = 0);

  Problem::Proxy GetProblemBuilder(fmt::StringRef stub);

  // Solves a problem and report solutions via the solution handler.
  int Solve(Problem &problem, SolutionHandler &sh);

  int Solve(internal::ASLBuilder &builder, SolutionHandler &sh) {
    Problem problem(builder.GetProblem());
    return Solve(problem, sh);
  }
};
}

#endif  // MP_ASL_SOLVER_H_
