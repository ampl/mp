/*
 Basic solver tests.

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

#include "gtest/gtest.h"
#include "solvers/util/solver.h"
#include "tests/args.h"
#include "tests/config.h"
#include "tests/solution_handler.h"
#include "tests/util.h"

#ifdef WIN32
# define putenv _putenv
#endif

using ampl::BasicSolver;
using ampl::Problem;
using ampl::Solver;

namespace {

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

struct TestSolver : BasicSolver {
  TestSolver(const char *name, const char *long_name = 0, long date = 0)
  : BasicSolver(name, long_name, date) {}

  void set_long_name(const char *name) {
    BasicSolver::set_long_name(name);
  }

  void set_version(const char *version) {
    BasicSolver::set_version(version);
  }

  bool ParseOptions(char **argv, unsigned flags = BasicSolver::NO_OPTION_ECHO) {
    return BasicSolver::ParseOptions(argv, flags);
  }

  void AddKeyword() {
    BasicSolver::AddKeyword("testopt", "A Test Option", TestKeywordFunc, INFO);
  }

  void Solve(Problem &) {}
};
}

TEST(SolverTest, ObjPrec) {
  double value = 12.3456789123456789;
  char buffer[64];
  sprintf(buffer, "%.*g", obj_prec(), value);
  EXPECT_EQ(buffer, str(fmt::Format("{}") << ampl::ObjPrec(value)));
}

TEST(SolverTest, Format) {
  EXPECT_EQ(
    "     This is a very long option description that should be indented and\n"
    "     wrapped.\n",
    ampl::internal::Format(
          "This is a very long option description "
          "that should be indented and wrapped.", 5));
}

TEST(SolverTest, BasicSolverCtor) {
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

TEST(SolverTest, BasicSolverVirtualDtor) {
  bool destroyed = false;
  class DtorTestSolver : public BasicSolver {
   private:
    bool &destroyed_;

   public:
    DtorTestSolver(bool &destroyed)
    : BasicSolver("test", 0, 0), destroyed_(destroyed) {}
    ~DtorTestSolver() { destroyed_ = true; }
    void Solve(Problem &) {}
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
    s.ProcessArgs(args);
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
  EXPECT_EXIT({
    FILE *f = freopen("out", "w", stdout);
    s.ProcessArgs(args);
    fclose(f);
  }, ::testing::ExitedWithCode(0), "");
  fmt::Formatter format;
  format("Test Solver ({}), ASL({})\n") << sysdetails_ASL << ASLdate_ASL;
  EXPECT_EQ(format.str(), ReadFile("out"));
}

TEST(SolverTest, VersionWithDate) {
  TestSolver s("testsolver", "Test Solver", 20121227);
  Args args("program-name", "-v");
  EXPECT_EXIT({
    FILE *f = freopen("out", "w", stdout);
    s.ProcessArgs(args);
    fclose(f);
  }, ::testing::ExitedWithCode(0), "");
  fmt::Formatter format;
  format("Test Solver ({}), driver(20121227), ASL({})\n")
    << sysdetails_ASL << ASLdate_ASL;
  EXPECT_EQ(format.str(), ReadFile("out"));
}

TEST(SolverTest, SetVersion) {
  TestSolver s("testsolver", "Test Solver");
  const char *VERSION = "Solver Version 3.0";
  s.set_version(VERSION);
  EXPECT_STREQ(VERSION, s.version());
  Args args("program-name", "-v");
  EXPECT_EXIT({
    FILE *f = freopen("out", "w", stdout);
    s.ProcessArgs(args);
    fclose(f);
  }, ::testing::ExitedWithCode(0), "");
  fmt::Formatter format;
  format("{} ({}), ASL({})\n") << VERSION << sysdetails_ASL << ASLdate_ASL;
  EXPECT_EQ(format.str(), ReadFile("out"));
}

TEST(SolverTest, OptionsVar) {
  TestSolver s("testsolver");
  char options[] = "testsolver_options=wantsol=9";
  putenv(options);
  EXPECT_EQ(0, s.wantsol());
  EXPECT_TRUE(s.ParseOptions(Args(0)));
  EXPECT_EQ(9, s.wantsol());
}

TEST(SolverTest, ErrorHandler) {
  struct TestErrorHandler : ampl::ErrorHandler {
    std::string message;

    virtual ~TestErrorHandler() {}
    void HandleError(fmt::StringRef message) {
      this->message = message;
    }
  };

  TestErrorHandler eh;
  TestSolver s("test");
  s.set_error_handler(&eh);
  EXPECT_TRUE(&eh == s.error_handler());
  s.ReportError("test message");
  EXPECT_EQ("test message", eh.message);
}

TEST(SolverTest, SolutionHandler) {
  TestSolutionHandler sh;
  TestSolver s("test");
  s.set_solution_handler(&sh);
  EXPECT_TRUE(&sh == s.solution_handler());
  double primal = 0, dual = 0, obj = 42;
  s.HandleSolution("test message", &primal, &dual, obj);
  EXPECT_EQ(&s, sh.solver());
  EXPECT_EQ("test message", sh.message());
  EXPECT_EQ(&primal, sh.primal());
  EXPECT_EQ(&dual, sh.dual());
  EXPECT_EQ(42.0, sh.obj_value());
}

TEST(SolverTest, ReadProblem) {
  Args args("testprogram", "data/objconst.nl");
  TestSolver s("test");
  EXPECT_EQ(0, s.problem().num_vars());
  EXPECT_TRUE(s.ProcessArgs(args));
  EXPECT_EQ(1, s.problem().num_vars());
}

TEST(SolverTest, ReadProblemNoStub) {
  StderrRedirect redirect("out");
  Args args("testprogram");
  TestSolver s("test");
  EXPECT_EQ(0, s.problem().num_vars());
  EXPECT_FALSE(s.ProcessArgs(args));
  EXPECT_EQ(0, s.problem().num_vars());
}

TEST(SolverTest, ReadProblemError) {
  Args args("testprogram", "nonexistent");
  TestSolver s("test");
  EXPECT_EXIT({
    Stderr = stderr;
    s.ProcessArgs(args);
  }, ::testing::ExitedWithCode(1), "testprogram: can't open nonexistent.nl");
}

TEST(SolverTest, ReadingMinOrMaxWithZeroArgsFails) {
  const char *names[] = {"min", "max"};
  for (size_t i = 0, n = sizeof(names) / sizeof(*names); i < n; ++i) {
    std::string stub = str(fmt::Format("data/{}-with-zero-args") << names[i]);
    EXPECT_EXIT({
      Stderr = stderr;
      Problem p;
      p.Read(stub);
    }, ::testing::ExitedWithCode(1),
        c_str(fmt::Format("bad line 13 of {}.nl: 0") << stub));
  }
}

TEST(SolverTest, ReportError) {
  TestSolver s("test");
  EXPECT_EXIT({
    s.ReportError("File not found: {}") << "somefile";
    exit(0);
  }, ::testing::ExitedWithCode(0), "File not found: somefile");
}

TEST(SolverTest, OptionEcho) {
  EXPECT_EXIT({
    TestSolver s("test");
    FILE *f = freopen("out", "w", stdout);
    s.ParseOptions(Args("wantsol=3"));
    s.ParseOptions(Args("wantsol=5"), 0);
    s.ParseOptions(Args("wantsol=9"));
    fclose(f);
    exit(0);
  }, ::testing::ExitedWithCode(0), "");
  EXPECT_EQ("wantsol=5\n", ReadFile("out"));
}

TEST(SolverTest, ParseOptions) {
  TestSolver s("testsolver");
  s.AddKeyword();
  s.ParseOptions(Args("wantsol=5", "testopt=42"));
  EXPECT_EQ(5, s.wantsol());
}

TEST(SolverTest, SolverWithDefaultOptionHandler) {
  struct TestSolver : Solver<TestSolver> {
    TestSolver() : Solver<TestSolver>("testsolver") {}
    void Solve(Problem &) {}
  };
  TestSolver s;
  EXPECT_STREQ("testsolver", s.name());
  s.ParseOptions(Args("wantsol=3"), BasicSolver::NO_OPTION_ECHO);
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

  void Throw(const char *, int value) {
    if (value == 1)
      throw 1;
    throw std::runtime_error("Test exception in handler");
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
    AddIntOption("throw", "", &OptSolver::Throw);
  }

  bool ParseOptions(char **argv, unsigned flags = BasicSolver::NO_OPTION_ECHO) {
    return BasicSolver::ParseOptions(argv, flags);
  }

  void Solve(Problem &) {}
};

TEST(SolverTest, SolverOptions) {
  OptSolver s;
  EXPECT_TRUE(s.ParseOptions(Args("intopt1=3", "intopt2=7")));
  EXPECT_EQ(3, s.intopt1);
  EXPECT_EQ(7, s.intopt2);
  EXPECT_TRUE(s.ParseOptions(Args("dblopt2=1.3", "dblopt1=5.4")));
  EXPECT_EQ(5.4, s.dblopt1);
  EXPECT_EQ(1.3, s.dblopt2);
  EXPECT_TRUE(s.ParseOptions(Args("stropt1=abc", "stropt2=def")));
  EXPECT_EQ("abc", s.stropt1);
  EXPECT_EQ("def", s.stropt2);
}

struct TestSolver2 : Solver<TestSolver2> {
  int answer;

  void SetAnswer(const char *name, int value) {
    EXPECT_STREQ("answer", name);
    answer = value;
  }

  TestSolver2() : Solver<TestSolver2>("test") {
    AddIntOption("answer", "The answer to life the universe and everything",
        &TestSolver2::SetAnswer);
  }

  void Solve(Problem &) {}
};

TEST(SolverTest, SeparateOptionHandler) {
  TestSolver2 s;
  EXPECT_TRUE(s.ParseOptions(
      Args("answer=42"), BasicSolver::NO_OPTION_ECHO));
  EXPECT_EQ(42, s.answer);
}

TEST(SolverTest, OptionParseError) {
  OptSolver s;
  EXPECT_FALSE(s.ParseOptions(Args("badopt=3")));
  EXPECT_EXIT({
    FILE *f = freopen("out", "w", stdout);
    s.ParseOptions(Args("badopt=3"));
    fclose(f);
    exit(0);
  }, ::testing::ExitedWithCode(0), "");
  EXPECT_EQ("Unknown keyword \"badopt\"\n", ReadFile("out"));
}

struct TestErrorHandler : ampl::ErrorHandler {
  std::vector<std::string> errors_;

  virtual ~TestErrorHandler() {}
  void HandleError(fmt::StringRef message) {
    errors_.push_back(message);
  }
};

TEST(SolverTest, ExceptionInOptionHandler) {
  OptSolver s;
  TestErrorHandler handler;
  s.set_error_handler(&handler);
  EXPECT_FALSE(s.ParseOptions(Args("throw=1")));
  EXPECT_EQ(1u, handler.errors_.size());
  EXPECT_EQ("Unknown exception in option handler", handler.errors_[0]);
  EXPECT_FALSE(s.ParseOptions(Args("throw=2")));
  EXPECT_EQ(2u, handler.errors_.size());
  EXPECT_EQ("Test exception in handler", handler.errors_[1]);
}

TEST(SolverTest, ProcessArgsReadsProblem) {
  Args args("testprogram", "data/objconst.nl");
  OptSolver s;
  EXPECT_EQ(0, s.problem().num_vars());
  EXPECT_TRUE(s.ProcessArgs(args));
  EXPECT_EQ(1, s.problem().num_vars());
}

TEST(SolverTest, ProcessArgsWithouStub) {
  StderrRedirect redirect("out");
  Args args("testprogram");
  OptSolver s;
  EXPECT_EQ(0, s.problem().num_vars());
  EXPECT_FALSE(s.ProcessArgs(args));
  EXPECT_EQ(0, s.problem().num_vars());
}

TEST(SolverTest, ProcessArgsError) {
  Args args("testprogram", "nonexistent");
  OptSolver s;
  EXPECT_EXIT({
    Stderr = stderr;
    s.ProcessArgs(args);
  }, ::testing::ExitedWithCode(1), "testprogram: can't open nonexistent.nl");
}

TEST(SolverTest, ProcessArgsParsesSolverOptions) {
  OptSolver s;
  Args args("testprogram", "data/objconst.nl", "intopt1=3");
  EXPECT_TRUE(s.ProcessArgs(args, BasicSolver::NO_OPTION_ECHO));
  EXPECT_EQ(3, s.intopt1);
}

TEST(SolverTest, SignalHandler) {
  std::signal(SIGINT, SIG_DFL);
  TestSolver s("testsolver");
  EXPECT_EXIT({
    FILE *f = freopen("out", "w", stdout);
    ampl::SignalHandler sh(s);
    fmt::Print("{}") << ampl::SignalHandler::stop();
    std::fflush(stdout);
    std::raise(SIGINT);
    fmt::Print("{}") << ampl::SignalHandler::stop();
    fclose(f);
    exit(0);
  }, ::testing::ExitedWithCode(0), "");
  EXPECT_EQ("0\n<BREAK> (testsolver)\n1", ReadFile("out"));
}

TEST(SolverTest, SignalHandlerExitOnTwoSIGINTs) {
  std::signal(SIGINT, SIG_DFL);
  TestSolver s("testsolver");
  EXPECT_EXIT({
    ampl::SignalHandler sh(s);
    FILE *f = freopen("out", "w", stdout);
    std::raise(SIGINT);
    std::raise(SIGINT);
    std::fclose(f); // Unreachable, but silences a warning.
  }, ::testing::ExitedWithCode(1), "");
  EXPECT_EQ("\n<BREAK> (testsolver)\n\n<BREAK> (testsolver)\n",
      ReadFile("out"));
}
