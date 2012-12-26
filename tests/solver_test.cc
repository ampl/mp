/*
 Solver tests.

 Copyright (C) 2012 AMPL Optimization LLC

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization LLC disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Victor Zverovich
 */

#include <fstream>

#include "gtest/gtest.h"
#include "solvers/util/solver.h"
#include "tests/args.h"

using ampl::SolverBase;

struct TestSolver : SolverBase {
  TestSolver(const char *name) : SolverBase(name) {}
};

TEST(SolverTest, SolverBaseCtor) {
  TestSolver d("testsolver");
  EXPECT_EQ(0, d.problem().num_vars());
  EXPECT_STREQ("testsolver", d.name());
  EXPECT_STREQ("testsolver", d.long_name());
  EXPECT_STREQ("testsolver_options", d.options_var_name());
  EXPECT_STREQ("testsolver", d.version());
  EXPECT_EQ(0, d.date());
  EXPECT_EQ(0, d.flags());
  EXPECT_EQ(0, d.wantsol());
}

TEST(SolverTest, SolverNameInUsage) {
  FILE *saved_stderr = Stderr;
  Stderr = fopen("out", "w");

  TestSolver d("solver-name");
  Args args("program-name");
  char **argv = args;
  d.ReadProblem(argv);

  fclose(Stderr);
  Stderr = saved_stderr;

  std::ifstream ifs("out");
  enum { BUFFER_SIZE = 4096 };
  char buffer[BUFFER_SIZE];
  ifs.read(buffer, BUFFER_SIZE);
  std::string output(buffer, static_cast<std::string::size_type>(ifs.gcount()));
  std::string expected = "usage: solver-name ";
  EXPECT_EQ(expected, output.substr(0, expected.size()));
}

class DtorTestSolver : public SolverBase {
 private:
  bool &destroyed_;

 public:
  DtorTestSolver(bool &destroyed) : SolverBase("test"), destroyed_(destroyed) {}
  ~DtorTestSolver() { destroyed_ = true; }
};

TEST(SolverTest, SolverBaseVirtualDtor) {
  bool destroyed = false;
  (DtorTestSolver(destroyed));
  EXPECT_TRUE(destroyed);
}

// TODO: test version, date, setters, solution handler, etc.
