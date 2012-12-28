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

#ifdef WIN32
# define putenv _putenv
#endif

using ampl::SolverBase;
using ampl::Solver;

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

const void *const INFO = "";

char *TestKeywordFunc(Option_Info *oi, keyword *kw, char *value) {
  EXPECT_STREQ("testsolver", oi->sname);
  EXPECT_TRUE(kw->info == INFO);
  EXPECT_TRUE(kw->kf == TestKeywordFunc);
  EXPECT_STREQ("testopt", kw->name);
  EXPECT_STREQ("A Test Option", kw->desc);
  EXPECT_EQ("42", std::string(value, 2));
  return value + 2;
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
    return DoParseOptions(argv);
  }

  void EnableOptionEcho() {
    SolverBase::EnableOptionEcho();
  }

  void DisableOptionEcho() {
    SolverBase::DisableOptionEcho();
  }

  void AddKeyword() {
    SolverBase::AddKeyword("testopt", "A Test Option", TestKeywordFunc, INFO);
  }

  std::string FormatDescription(const char *description) {
    return SolverBase::FormatDescription(description);
  }
};

// Redirects Stderr to a file.
class StderrRedirect {
 private:
  FILE *saved_stderr;

 public:
  StderrRedirect(const char *filename) : saved_stderr(Stderr) {
    Stderr = fopen(filename, "w");
  }

  ~StderrRedirect() {
    fclose(Stderr);
    Stderr = saved_stderr;
  }
};
}

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

TEST(SolverTest, SolverBaseVirtualDtor) {
  bool destroyed = false;
  class DtorTestSolver : public SolverBase {
   private:
    bool &destroyed_;

   public:
    DtorTestSolver(bool &destroyed)
    : SolverBase("test", 0, 0), destroyed_(destroyed) {}
    ~DtorTestSolver() { destroyed_ = true; }
  };
  (DtorTestSolver(destroyed));
  EXPECT_TRUE(destroyed);
}

TEST(SolverTest, NameInUsage) {
  {
    StderrRedirect redirect("out");
    TestSolver s("solver-name", "long-solver-name");
    s.set_version("solver-version");
    Args args("program-name");
    char **argv = args;
    s.ReadProblem(argv);
  }
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
    FILE *f = freopen("out", "w", stdout);
    s.ReadProblem(argv);
    fclose(f);
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
    FILE *f = freopen("out", "w", stdout);
    s.ReadProblem(argv);
    fclose(f);
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
    FILE *f = freopen("out", "w", stdout);
    s.ReadProblem(argv);
    fclose(f);
  }, ::testing::ExitedWithCode(0), "");
  fmt::Formatter format;
  format("{0} ({1}), ASL({2})\n") << VERSION << sysdetails_ASL << ASLdate_ASL;
  EXPECT_EQ(format.str(), ReadFile("out"));
}

TEST(SolverTest, OptionsVar) {
  TestSolver s("testsolver");
  char options[] = "testsolver_options=wantsol=9";
  putenv(options);
  s.DisableOptionEcho();
  EXPECT_EQ(0, s.wantsol());
  EXPECT_TRUE(s.ParseOptions(Args(0)));
  EXPECT_EQ(9, s.wantsol());
}

TEST(SolverTest, SolutionHandler) {
  struct TestSolutionHandler : ampl::SolutionHandler {
    SolverBase *solver;
    std::string message;
    const double *primal;
    const double *dual;
    double obj_value;

    TestSolutionHandler() : solver(0), primal(0), dual(0), obj_value(0) {}

    void HandleSolution(SolverBase &s, fmt::StringRef message,
            const double *primal, const double *dual, double obj_value) {
      solver = &s;
      this->message = message;
      this->primal = primal;
      this->dual = dual;
      this->obj_value = obj_value;
    }
  };

  TestSolutionHandler sh;
  TestSolver s("test");
  s.set_solution_handler(&sh);
  EXPECT_TRUE(&sh == s.solution_handler());
  double primal = 0, dual = 0, obj = 42;
  s.HandleSolution("test message", &primal, &dual, obj);
  EXPECT_EQ(&s, sh.solver);
  EXPECT_EQ("test message", sh.message);
  EXPECT_EQ(&primal, sh.primal);
  EXPECT_EQ(&dual, sh.dual);
  EXPECT_EQ(42.0, sh.obj_value);
}

TEST(SolverTest, ReadProblem) {
  Args args("testprogram", "data/objconst.nl");
  char **argv = args;
  TestSolver s("test");
  EXPECT_EQ(0, s.problem().num_vars());
  EXPECT_TRUE(s.ReadProblem(argv));
  EXPECT_EQ(1, s.problem().num_vars());
}

TEST(SolverTest, ReadProblemNoStub) {
  StderrRedirect redirect("out");
  Args args("testprogram");
  char **argv = args;
  TestSolver s("test");
  EXPECT_EQ(0, s.problem().num_vars());
  EXPECT_FALSE(s.ReadProblem(argv));
  EXPECT_EQ(0, s.problem().num_vars());
}

TEST(SolverTest, ReadProblemError) {
  Args args("testprogram", "nonexistent");
  char **argv = args;
  TestSolver s("test");
  EXPECT_EXIT({
    Stderr = stderr;
    s.ReadProblem(argv);
  }, ::testing::ExitedWithCode(1), "testprogram: can't open nonexistent.nl");
}

TEST(SolverTest, ReportError) {
  TestSolver s("test");
  EXPECT_EXIT({
    s.ReportError("File not found: {0}") << "somefile";
    exit(0);
  }, ::testing::ExitedWithCode(0), "File not found: somefile");
}

TEST(SolverTest, OptionEcho) {
  EXPECT_EXIT({
    TestSolver s("test");
    FILE *f = freopen("out", "w", stdout);
    s.DisableOptionEcho();
    s.ParseOptions(Args("wantsol=3"));
    s.EnableOptionEcho();
    s.ParseOptions(Args("wantsol=5"));
    s.DisableOptionEcho();
    s.ParseOptions(Args("wantsol=9"));
    fclose(f);
    exit(0);
  }, ::testing::ExitedWithCode(0), "");
  EXPECT_EQ("wantsol=5\n", ReadFile("out"));
}

TEST(SolverTest, ParseOptions) {
  TestSolver s("testsolver");
  s.DisableOptionEcho();
  s.AddKeyword();
  s.ParseOptions(Args("wantsol=5", "testopt=42"));
  EXPECT_EQ(5, s.wantsol());
}

TEST(SolverTest, FormatDescription) {
  TestSolver s("testsolver");
  EXPECT_EQ(
    "\n"
    "      This is a very long option description that should be indented and\n"
    "      wrapped.\n",
    s.FormatDescription(
          "This is a very long option description "
          "that should be indented and wrapped."));
}

struct DummyOptionHandler {};

TEST(SolverTest, SolverWithDefaultOptionHandler) {
  Solver<DummyOptionHandler> s("testsolver");
  EXPECT_STREQ("testsolver", s.name());
  s.DisableOptionEcho();
  DummyOptionHandler handler;
  s.ParseOptions(Args("wantsol=3"), handler);
  EXPECT_EQ(3, s.wantsol());
}

struct OptSolver : Solver<OptSolver> {
  int intopt1;
  int intopt2;
  double dblopt1;
  double dblopt2;
  std::string stropt1;
  std::string stropt2;

  enum Tag { A, B, C, D };

  void SetIntOption(const char *name, int value) {
    EXPECT_STREQ("intopt1", name);
    intopt1 = value;
  }

  void SetIntOptionWithTag(const char *name, int value, Tag tag) {
    EXPECT_STREQ("intopt2", name);
    intopt2 = value;
    EXPECT_EQ(B, tag);
  }

  void SetDblOption(const char *name, double value) {
    EXPECT_STREQ("dblopt1", name);
    dblopt1 = value;
  }

  void SetDblOptionWithTag(const char *name, double value, Tag tag) {
    EXPECT_STREQ("dblopt2", name);
    dblopt2 = value;
    EXPECT_EQ(C, tag);
  }

  void SetStrOption(const char *name, const char *value) {
    EXPECT_STREQ("stropt1", name);
    stropt1 = value;
  }

  void SetStrOptionWithTag(const char *name, const char *value, Tag tag) {
    EXPECT_STREQ("stropt2", name);
    stropt2 = value;
    EXPECT_EQ(D, tag);
  }

  OptSolver()
  : Solver<OptSolver>("test"), intopt1(0), intopt2(0), dblopt1(0), dblopt2(0) {
    AddIntOption("intopt1", "Integer option 1", &OptSolver::SetIntOption);
    AddIntOption("intopt2", "Integer option 2",
        &OptSolver::SetIntOptionWithTag, B);
    AddDblOption("dblopt1", "Double option 1", &OptSolver::SetDblOption);
    AddDblOption("dblopt2", "Double option 2",
        &OptSolver::SetDblOptionWithTag, C);
    AddStrOption("stropt1", "Double option 1", &OptSolver::SetStrOption);
    AddStrOption("stropt2", "Double option 2",
        &OptSolver::SetStrOptionWithTag, D);
  }
};

TEST(SolverTest, SolverOptions) {
  OptSolver s;
  s.DisableOptionEcho();
  s.ParseOptions(Args("intopt1=3", "intopt2=7"), s);
  EXPECT_EQ(3, s.intopt1);
  EXPECT_EQ(7, s.intopt2);
  s.ParseOptions(Args("dblopt2=1.3", "dblopt1=5.4"), s);
  EXPECT_EQ(5.4, s.dblopt1);
  EXPECT_EQ(1.3, s.dblopt2);
  s.ParseOptions(Args("stropt1=abc", "stropt2=def"), s);
  EXPECT_EQ("abc", s.stropt1);
  EXPECT_EQ("def", s.stropt2);
}

// TODO: test parse errors and separate handler
