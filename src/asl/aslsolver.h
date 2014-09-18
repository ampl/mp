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

class SolutionWriter : public SolutionHandler {
 private:
  Solver &solver_;

  // The number of feasible solutions found so far.
  int num_solutions_;

 public:
  explicit SolutionWriter(Solver &s) : solver_(s), num_solutions_(0) {}

  void HandleFeasibleSolution(Problem &p, fmt::StringRef message,
        const double *values, const double *dual_values, double);
  void HandleSolution(Problem &p, fmt::StringRef message,
        const double *values, const double *dual_values, double);
};

class ASLSolver : public SolverImpl<internal::ASLBuilder> {
 private:
  void RegisterSuffixes(Problem &p);

 public:
  ASLSolver(fmt::StringRef name, fmt::StringRef long_name = 0,
            long date = 0, int flags = 0);

  // Processes command-line arguments, reads a problem from an .nl file
  // if the file name (stub) is specified and parses solver options.
  // Returns true if the arguments contain the file name and options were
  // parsed successfully; false otherwise.
  // If there was an  error parsing arguments or reading the problem
  // ProcessArgs will print an error message and call std::exit (this is
  // likely to change in the future version).
  bool ProcessArgs(char **&argv, Problem &p, unsigned flags = 0);

  // Solves a problem and report solutions via the solution handler.
  void Solve(Problem &p, SolutionHandler &sh);

  // Runs the solver.
  int Run(char **argv);
};
}

#endif  // MP_ASL_SOLVER_H_
