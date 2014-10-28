/*
 LocalSolver tests

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

#include "localsolver/localsolver.h"
#include "feature.h"

typedef mp::LocalSolver Solver;
enum { FEATURES = ~feature::PLTERM };

#define MP_THREAD_INTERRUPT 1

// Demo version of LocalSolver cannot handle a TSP problem with n > 9.
#define MP_TSP_SIZE 9

#include "nl-solver-test.h"

// Creates and solves a test problem.
void SolveTestProblem(mp::LocalSolver &s) {
  mp::LSProblemBuilder pb(s.GetProblemBuilder(""));
  MakeAllDiffProblem(pb);
  TestSolutionHandler sh;
  s.Solve(pb, sh);
}

// Creates a test LocalSolver model.
void CreateTestModel(localsolver::LocalSolver &s) {
  s.getParam().setVerbosity(0);
  localsolver::LSModel model = s.getModel();
  model.addObjective(model.createConstant(localsolver::lsint(0)),
                     localsolver::OD_Minimize);
  model.close();
}

TEST(LocalSolverTest, SeedOption) {
  struct TestLocalSolver : mp::LocalSolver {
    void DoSolve(localsolver::LocalSolver &s) {
      EXPECT_EQ(100, s.getParam().getSeed());
    }
  } solver;

  // Check that the default value and range agrees with LocalSolver.
  localsolver::LocalSolver ls;
  CreateTestModel(ls);
  localsolver::LSParam param = ls.getParam();
  EXPECT_EQ(param.getSeed(), TestLocalSolver().GetIntOption("seed"));
  param.setSeed(0);
  param.setSeed(INT_MAX);
  EXPECT_THROW(param.setSeed(-1), localsolver::LSException);

  solver.SetIntOption("seed", 100);
  EXPECT_EQ(100, solver.GetIntOption("seed"));
  SolveTestProblem(solver);
  solver.SetIntOption("seed", 0);
  EXPECT_THROW(solver.SetIntOption("seed", -1), mp::InvalidOptionValue);
  EXPECT_THROW(solver.SetStrOption("seed", "oops"), mp::OptionError);
}

TEST(LocalSolverTest, ThreadsOption) {
  struct TestLocalSolver : mp::LocalSolver {
    void DoSolve(localsolver::LocalSolver &s) {
      EXPECT_EQ(10, s.getParam().getNbThreads());
    }
  } solver;

  // Check that the default value agrees with LocalSolver.
  localsolver::LocalSolver ls;
  CreateTestModel(ls);
  localsolver::LSParam param = ls.getParam();
  EXPECT_EQ(param.getNbThreads(), TestLocalSolver().GetIntOption("threads"));
  param.setNbThreads(1);
  param.setNbThreads(1024);
  EXPECT_THROW(param.setNbThreads(0), localsolver::LSException);
  EXPECT_THROW(param.setNbThreads(1025), localsolver::LSException);

  solver.SetIntOption("threads", 10);
  EXPECT_EQ(10, solver.GetIntOption("threads"));
  SolveTestProblem(solver);
  solver.SetIntOption("threads", 1);
  EXPECT_THROW(solver.SetIntOption("threads", 0), mp::InvalidOptionValue);
  EXPECT_THROW(solver.SetStrOption("threads", "oops"), mp::OptionError);
}

// TODO: test solver options annealing_level, verbosity,
//       time_between_displays, logfile, timelimit, iterlimit
