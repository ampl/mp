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

namespace {
std::string ReadFile(const char *name) {
  std::string data;
  std::ifstream ifs(name);
  enum { BUFFER_SIZE = 4096 };
  char buffer[BUFFER_SIZE];
  do {
    ifs.read(buffer, BUFFER_SIZE);
    data.append(buffer, static_cast<std::string::size_type>(ifs.gcount()));
  } while (ifs);
  return data;
}
}

struct TestSolver : SolverBase {
  TestSolver(const char *name, const char *long_name = 0, long date = 0)
  : SolverBase(name, long_name, date) {}

  void set_long_name(const char *name) {
    SolverBase::set_long_name(name);
  }

  void set_version(const char *version) {
    SolverBase::set_version(version);
  }

  bool ParseOptions(char **argv) {
    return SolverBase::ParseOptions(argv);
  }
};

TEST(SolverTest, SolverBaseCtor) {
  TestSolver s("testsolver");
  EXPECT_EQ(0, s.problem().num_vars());
  EXPECT_STREQ("testsolver", s.name());
  EXPECT_STREQ("testsolver", s.long_name());
  EXPECT_STREQ("testsolver_options", s.options_var_name());
  EXPECT_STREQ("testsolver", s.version());
  EXPECT_EQ(0, s.date());
  EXPECT_EQ(0, s.flags());
  EXPECT_EQ(0, s.wantsol());
}

class DtorTestSolver : public SolverBase {
 private:
  bool &destroyed_;

 public:
  DtorTestSolver(bool &destroyed)
  : SolverBase("test", 0, 0), destroyed_(destroyed) {}
  ~DtorTestSolver() { destroyed_ = true; }
};

TEST(SolverTest, SolverBaseVirtualDtor) {
  bool destroyed = false;
  (DtorTestSolver(destroyed));
  EXPECT_TRUE(destroyed);
}

TEST(SolverTest, NameInUsage) {
  FILE *saved_stderr = Stderr;
  Stderr = fopen("out", "w");
  TestSolver s("solver-name", "long-solver-name");
  s.set_version("solver-version");
  Args args("program-name");
  char **argv = args;
  s.ReadProblem(argv);
  fclose(Stderr);
  Stderr = saved_stderr;
  std::string usage = "usage: solver-name ";
  EXPECT_EQ(usage, ReadFile("out").substr(0, usage.size()));
}

TEST(SolverTest, LongName) {
  EXPECT_STREQ("solver-name", TestSolver("solver-name").long_name());
  EXPECT_STREQ("long-solver-name",
      TestSolver("solver-name", "long-solver-name").long_name());
  TestSolver s("solver-name");
  s.set_long_name("another-name");
  EXPECT_STREQ("another-name", s.long_name());
}

TEST(SolverTest, Version) {
  TestSolver s("testsolver", "Test Solver");
  Args args("program-name", "-v");
  char **argv = args;
  EXPECT_EXIT({
    freopen("out", "w", stdout);
    s.ReadProblem(argv);
  }, ::testing::ExitedWithCode(0), "");
  fmt::Formatter format;
  format("Test Solver ({0}), ASL({1})\n") << sysdetails_ASL << ASLdate_ASL;
  EXPECT_EQ(format.str(), ReadFile("out"));
}

TEST(SolverTest, VersionWithDate) {
  TestSolver s("testsolver", "Test Solver", 20121227);
  Args args("program-name", "-v");
  char **argv = args;
  EXPECT_EXIT({
    freopen("out", "w", stdout);
    s.ReadProblem(argv);
  }, ::testing::ExitedWithCode(0), "");
  fmt::Formatter format;
  format("Test Solver ({0}), driver(20121227), ASL({1})\n")
    << sysdetails_ASL << ASLdate_ASL;
  EXPECT_EQ(format.str(), ReadFile("out"));
}

TEST(SolverTest, SetVersion) {
  TestSolver s("testsolver", "Test Solver");
  const char *VERSION = "Solver Version 3.0";
  s.set_version(VERSION);
  EXPECT_STREQ(VERSION, s.version());
  Args args("program-name", "-v");
  char **argv = args;
  EXPECT_EXIT({
    freopen("out", "w", stdout);
    s.ReadProblem(argv);
  }, ::testing::ExitedWithCode(0), "");
  fmt::Formatter format;
  format("{0} ({1}), ASL({2})\n") << VERSION << sysdetails_ASL << ASLdate_ASL;
  EXPECT_EQ(format.str(), ReadFile("out"));
}

TEST(SolverTest, OptionsVar) {
  TestSolver s("testsolver");
  char options[] = "testsolver_options=wantsol=9";
  putenv(options);
  EXPECT_EQ(0, s.wantsol());
  EXPECT_TRUE(s.ParseOptions(Args("program-name")));
  EXPECT_EQ(9, s.wantsol());
}

// TODO: test setters, solution handler, etc.
