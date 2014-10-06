/*
 ASLSolver tests.

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

#include "gtest/gtest.h"
#include "asl/aslsolver.h"

struct TestSolver : mp::ASLSolver {
  TestSolver() : ASLSolver("testsolver") {
    AddSuffix("answer", 0, mp::suf::VAR | mp::suf::OUTONLY, 0);
  }
  int DoSolve(mp::Problem &, mp::SolutionHandler &) { return 0; }
};

TEST(ASLSolverTest, RegisterSuffixes) {
  TestSolver s;
  mp::internal::ASLBuilder builder(s.GetProblemBuilder(""));
  mp::Problem p(builder.GetProblem());
  EXPECT_TRUE(p.suffixes(mp::suf::VAR).Find("answer"));
}
