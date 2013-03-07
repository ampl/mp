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

#include <csignal>
#include <fstream>

#include "gtest/gtest.h"
#include "solvers/util/solver.h"
#include "tests/args.h"
#include "tests/config.h"
#include "tests/solution_handler.h"
#include "tests/util.h"

#ifdef WIN32
# define putenv _putenv
#endif

using ampl::LinearObjExpr;
using ampl::LinearConExpr;
using ampl::Problem;
using ampl::ProblemChanges;
using ampl::BasicSolver;
using ampl::Solution;
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
    return DoParseOptions(argv, flags);
  }

  void AddKeyword() {
    BasicSolver::AddKeyword("testopt", "A Test Option", TestKeywordFunc, INFO);
  }

  std::string FormatDescription(const char *description) {
    return BasicSolver::FormatDescription(description);
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

class CBCPath {
 private:
  std::string path_;

 public:
  CBCPath() : path_("../solvers/cbc/bin/cbc") {
#ifdef WIN32
    std::replace(path_.begin(), path_.end(), '/', '\\');
#endif
  }

  operator const char *() const { return path_.c_str(); }
};

const CBCPath CBC_PATH;
}

TEST(SolutionTest, DefaultCtor) {
  Solution s;
  EXPECT_EQ(Solution::UNKNOWN, s.status());
  EXPECT_EQ(-1, s.solve_code());
  EXPECT_EQ(0, s.num_vars());
  EXPECT_EQ(0, s.num_cons());
  EXPECT_EQ(0, s.values());
  EXPECT_EQ(0, s.dual_values());
}

TEST(SolutionTest, Read) {
  WriteFile("test.sol", "test\n\n1\n3\n5\n7\n11\n");
  Solution s;
  s.Read("test", 3, 2);
  EXPECT_EQ(Solution::UNKNOWN, s.status());
  EXPECT_EQ(-1, s.solve_code());
  EXPECT_EQ(3, s.num_vars());
  EXPECT_EQ(2, s.num_cons());
  const double values[] = {5, 7, 11};
  for (int i = 0; i < 3; ++i) {
    EXPECT_EQ(values[i], s.value(i));
    EXPECT_EQ(values[i], s.values()[i]);
  }
  const double dual_values[] = {1, 3};
  for (int i = 0; i < 2; ++i) {
    EXPECT_EQ(dual_values[i], s.dual_value(i));
    EXPECT_EQ(dual_values[i], s.dual_values()[i]);
  }
}

TEST(SolutionTest, ReadError) {
  Solution s;
  StderrRedirect redirect("out");
  EXPECT_THROW(s.Read("nonexistent", 0, 0), ampl::Error);
}

TEST(SolutionTest, ReadEmpty) {
  WriteFile("test.sol", "test\n\n");
  Solution s;
  s.Read("test", 0, 0);
  EXPECT_EQ(0, s.num_vars());
  EXPECT_EQ(0, s.num_cons());
  EXPECT_EQ(0, s.solve_code());
}

TEST(SolutionTest, DoubleRead) {
  WriteFile("test.sol", "test\n\n1\n3\n5\n7\n11\n");
  Solution s;
  s.Read("test", 3, 2);
  WriteFile("test.sol", "test\n\n44\n22\n33\n");
  s.Read("test", 2, 1);
  EXPECT_EQ(2, s.num_vars());
  EXPECT_EQ(1, s.num_cons());
  EXPECT_EQ(22, s.value(0));
  EXPECT_EQ(33, s.value(1));
  EXPECT_EQ(44, s.dual_value(0));
}

TEST(SolutionTest, SolveCodes) {
  const Solution::Status STATES[] = {
      Solution::SOLVED,
      Solution::SOLVED_MAYBE,
      Solution::INFEASIBLE,
      Solution::UNBOUNDED,
      Solution::LIMIT,
      Solution::FAILURE
  };
  for (std::size_t i = 0; i < sizeof(STATES) / sizeof(*STATES); ++i) {
    {
      int solve_code = i * 100;
      WriteFile("test.sol",
          c_str(fmt::Format("test\n\n2\n2\nobjno 0 {}\n") << solve_code));
      Solution s;
      s.Read("test", 1, 1);
      EXPECT_EQ(STATES[i], s.status());
      EXPECT_EQ(solve_code, s.solve_code());
    }
    {
      int solve_code = i * 100 + 99;
      WriteFile("test.sol",
          c_str(fmt::Format("test\n\n2\n2\nobjno 0 {}\n") << solve_code));
      Solution s;
      s.Read("test", 1, 1);
      EXPECT_EQ(STATES[i], s.status());
      EXPECT_EQ(solve_code, s.solve_code());
    }
  }
  const double CODES[] = {-5, -1, 600, 1000};
  for (std::size_t i = 0; i < sizeof(CODES) / sizeof(*CODES); ++i) {
    WriteFile("test.sol",
        c_str(fmt::Format("test\n\n2\n2\nobjno 0 {}\n") << CODES[i]));
    Solution s;
    s.Read("test", 1, 1);
    EXPECT_EQ(Solution::UNKNOWN, s.status());
    EXPECT_EQ(CODES[i], s.solve_code());
  }
}

#ifndef NDEBUG
TEST(SolutionTest, BoundChecks) {
  WriteFile("test.sol", "test\n\n1\n3\n5\n7\n11\n");
  Solution s;
  s.Read("test", 3, 2);
  EXPECT_DEATH(s.value(-1), "Assertion");
  EXPECT_DEATH(s.value(3), "Assertion");
  EXPECT_DEATH(s.dual_value(-1), "Assertion");
  EXPECT_DEATH(s.dual_value(2), "Assertion");
}
#endif

TEST(SolutionTest, Swap) {
  WriteFile("test.sol", "test\n\n1\n3\n5\n7\n11\nobjno 0 10\n");
  Solution s1;
  s1.Read("test", 3, 2);
  WriteFile("test.sol", "test\n\n44\n22\n33\nobjno 0 20");
  Solution s2;
  s2.Read("test", 2, 1);
  s1.Swap(s2);

  EXPECT_EQ(20, s1.solve_code());
  EXPECT_EQ(2, s1.num_vars());
  EXPECT_EQ(1, s1.num_cons());
  EXPECT_EQ(22, s1.value(0));
  EXPECT_EQ(33, s1.value(1));
  EXPECT_EQ(44, s1.dual_value(0));

  EXPECT_EQ(10, s2.solve_code());
  EXPECT_EQ(3, s2.num_vars());
  EXPECT_EQ(2, s2.num_cons());
  EXPECT_EQ(5, s2.value(0));
  EXPECT_EQ(7, s2.value(1));
  EXPECT_EQ(11, s2.value(2));
  EXPECT_EQ(1, s2.dual_value(0));
  EXPECT_EQ(3, s2.dual_value(1));
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
    s.ReadProblem(args);
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
    s.ReadProblem(args);
    fclose(f);
  }, ::testing::ExitedWithCode(0), "");
  fmt::Formatter format;
  format("Test Solver ({0}), ASL({1})\n") << sysdetails_ASL << ASLdate_ASL;
  EXPECT_EQ(format.str(), ReadFile("out"));
}

TEST(SolverTest, VersionWithDate) {
  TestSolver s("testsolver", "Test Solver", 20121227);
  Args args("program-name", "-v");
  EXPECT_EXIT({
    FILE *f = freopen("out", "w", stdout);
    s.ReadProblem(args);
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
  EXPECT_EXIT({
    FILE *f = freopen("out", "w", stdout);
    s.ReadProblem(args);
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
  EXPECT_EQ(0, s.wantsol());
  EXPECT_TRUE(s.ParseOptions(Args(0)));
  EXPECT_EQ(9, s.wantsol());
}

TEST(SolverTest, ErrorHandler) {
  struct TestErrorHandler : ampl::ErrorHandler {
    std::string message;

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
  EXPECT_TRUE(s.ReadProblem(args));
  EXPECT_EQ(1, s.problem().num_vars());
}

TEST(SolverTest, ReadProblemNoStub) {
  StderrRedirect redirect("out");
  Args args("testprogram");
  TestSolver s("test");
  EXPECT_EQ(0, s.problem().num_vars());
  EXPECT_FALSE(s.ReadProblem(args));
  EXPECT_EQ(0, s.problem().num_vars());
}

TEST(SolverTest, ReadProblemError) {
  Args args("testprogram", "nonexistent");
  TestSolver s("test");
  EXPECT_EXIT({
    Stderr = stderr;
    s.ReadProblem(args);
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
  DummyOptionHandler handler;
  s.ParseOptions(Args("wantsol=3"), handler, BasicSolver::NO_OPTION_ECHO);
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

  bool ParseOptions(char **argv) {
    return Solver<OptSolver>::ParseOptions(
        argv, *this, BasicSolver::NO_OPTION_ECHO);
  }
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

struct TestOptionHandler {
  int answer;

  void SetAnswer(const char *name, int value) {
    EXPECT_STREQ("answer", name);
    answer = value;
  }
};

struct TestSolver2 : Solver<TestOptionHandler> {
  TestSolver2() : Solver<TestOptionHandler>("test") {
    AddIntOption("answer", "The answer to life the universe and everything",
        &TestOptionHandler::SetAnswer);
  }
};

TEST(SolverTest, SeparateOptionHandler) {
  TestSolver2 s;
  TestOptionHandler handler;
  EXPECT_TRUE(s.ParseOptions(
      Args("answer=42"), handler, BasicSolver::NO_OPTION_ECHO));
  EXPECT_EQ(42, handler.answer);
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
  EXPECT_TRUE(s.ProcessArgs(args, s));
  EXPECT_EQ(1, s.problem().num_vars());
}

TEST(SolverTest, ProcessArgsWithouStub) {
  StderrRedirect redirect("out");
  Args args("testprogram");
  OptSolver s;
  EXPECT_EQ(0, s.problem().num_vars());
  EXPECT_FALSE(s.ProcessArgs(args, s));
  EXPECT_EQ(0, s.problem().num_vars());
}

TEST(SolverTest, ProcessArgsError) {
  Args args("testprogram", "nonexistent");
  OptSolver s;
  EXPECT_EXIT({
    Stderr = stderr;
    s.ProcessArgs(args, s);
  }, ::testing::ExitedWithCode(1), "testprogram: can't open nonexistent.nl");
}

TEST(SolverTest, ProcessArgsParsesSolverOptions) {
  OptSolver s;
  Args args("testprogram", "data/objconst.nl", "intopt1=3");
  EXPECT_TRUE(s.ProcessArgs(args, s, BasicSolver::NO_OPTION_ECHO));
  EXPECT_EQ(3, s.intopt1);
}

TEST(SolverTest, ObjPrec) {
  double value = 12.3456789123456789;
  char buffer[64];
  sprintf(buffer, "%.*g", obj_prec(), value);
  EXPECT_EQ(buffer, str(fmt::Format("{}") << ampl::ObjPrec(value)));
}

TEST(SolverTest, EmptyProblem) {
  Problem p;
  EXPECT_EQ(0, p.num_vars());
  EXPECT_EQ(0, p.num_objs());
  EXPECT_EQ(0, p.num_cons());
  EXPECT_EQ(0, p.num_integer_vars());
  EXPECT_EQ(0, p.num_continuous_vars());
  EXPECT_EQ(0, p.num_nonlinear_objs());
  EXPECT_EQ(0, p.num_nonlinear_cons());
  EXPECT_EQ(0, p.num_logical_cons());
  EXPECT_EQ(-1, p.solve_code());
}

TEST(SolverTest, ProblemAccessors) {
  Args args("", "data/test.nl");
  TestSolver s("");
  EXPECT_TRUE(s.ReadProblem(args));
  Problem &p = s.problem();
  EXPECT_EQ(5, p.num_vars());
  EXPECT_EQ(19, p.num_objs());
  EXPECT_EQ(13, p.num_cons());
  EXPECT_EQ(2, p.num_integer_vars());
  EXPECT_EQ(3, p.num_continuous_vars());
  EXPECT_EQ(17, p.num_nonlinear_objs());
  EXPECT_EQ(11, p.num_nonlinear_cons());
  EXPECT_EQ(7, p.num_logical_cons());

  EXPECT_EQ(11, p.var_lb(0));
  EXPECT_EQ(15, p.var_lb(p.num_vars() - 1));
  EXPECT_EQ(21, p.var_ub(0));
  EXPECT_EQ(25, p.var_ub(p.num_vars() - 1));

  EXPECT_EQ(101, p.con_lb(0));
  EXPECT_EQ(113, p.con_lb(p.num_cons() - 1));
  EXPECT_EQ(201, p.con_ub(0));
  EXPECT_EQ(213, p.con_ub(p.num_cons() - 1));

  EXPECT_EQ(ampl::MIN, p.obj_type(0));
  EXPECT_EQ(ampl::MAX, p.obj_type(p.num_objs() - 1));

  {
    LinearObjExpr expr = p.linear_obj_expr(0);
    EXPECT_EQ(31, expr.begin()->coef());
    EXPECT_EQ(0, expr.begin()->var_index());
    expr = p.linear_obj_expr(p.num_objs() - 1);
    EXPECT_EQ(52, expr.begin()->coef());
    EXPECT_EQ(3, expr.begin()->var_index());
  }

  {
    LinearConExpr expr = p.linear_con_expr(0);
    EXPECT_EQ(61, expr.begin()->coef());
    EXPECT_EQ(0, expr.begin()->var_index());
    expr = p.linear_con_expr(p.num_cons() - 1);
    EXPECT_EQ(82, expr.begin()->coef());
    EXPECT_EQ(2, expr.begin()->var_index());
  }

  EXPECT_EQ(OP_sin, p.nonlinear_obj_expr(0).opcode());
  EXPECT_EQ(OP_cos, p.nonlinear_obj_expr(p.num_nonlinear_objs() - 1).opcode());

  EXPECT_EQ(OP_log, p.nonlinear_con_expr(0).opcode());
  EXPECT_EQ(OP_exp, p.nonlinear_con_expr(p.num_nonlinear_cons() - 1).opcode());

  EXPECT_EQ(NE, p.logical_con_expr(0).opcode());
  EXPECT_EQ(OPAND, p.logical_con_expr(p.num_logical_cons() - 1).opcode());

  EXPECT_EQ(-1, p.solve_code());
  p.set_solve_code(42);
  EXPECT_EQ(42, p.solve_code());
}

#ifndef NDEBUG
TEST(SolverTest, ProblemBoundChecks) {
  Args args("", "data/test.nl");
  TestSolver s("");
  EXPECT_TRUE(s.ReadProblem(args));
  Problem &p = s.problem();

  EXPECT_DEATH(p.var_lb(-1), "Assertion");
  EXPECT_DEATH(p.var_lb(p.num_vars()), "Assertion");
  EXPECT_DEATH(p.var_ub(-1), "Assertion");
  EXPECT_DEATH(p.var_ub(p.num_vars()), "Assertion");

  EXPECT_DEATH(p.con_lb(-1), "Assertion");
  EXPECT_DEATH(p.con_lb(p.num_cons()), "Assertion");
  EXPECT_DEATH(p.con_ub(-1), "Assertion");
  EXPECT_DEATH(p.con_ub(p.num_cons()), "Assertion");

  EXPECT_DEATH(p.obj_type(-1), "Assertion");
  EXPECT_DEATH(p.obj_type(p.num_objs()), "Assertion");

  EXPECT_DEATH(p.linear_obj_expr(-1), "Assertion");
  EXPECT_DEATH(p.linear_obj_expr(p.num_objs()), "Assertion");

  EXPECT_DEATH(p.linear_con_expr(-1), "Assertion");
  EXPECT_DEATH(p.linear_con_expr(p.num_cons()), "Assertion");

  EXPECT_DEATH(p.nonlinear_obj_expr(-1), "Assertion");
  EXPECT_DEATH(p.nonlinear_obj_expr(p.num_objs()), "Assertion");

  EXPECT_DEATH(p.nonlinear_con_expr(-1), "Assertion");
  EXPECT_DEATH(p.nonlinear_con_expr(p.num_cons()), "Assertion");

  EXPECT_DEATH(p.logical_con_expr(-1), "Assertion");
  EXPECT_DEATH(p.logical_con_expr(p.num_logical_cons()), "Assertion");
}
#endif

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

#ifdef HAVE_CBC
TEST(SolverTest, Solve) {
  Args args("", "data/simple.nl");
  TestSolver solver("testsolver");
  solver.ReadProblem(args);
  Solution s;
  solver.problem().Solve(CBC_PATH, s);
  EXPECT_EQ(2, s.num_vars());
  EXPECT_EQ(1, s.num_cons());
  EXPECT_EQ(2, s.value(0));
  EXPECT_NEAR(0, s.value(1), 1e-5);
  EXPECT_EQ(0.5, s.dual_value(0));
}

TEST(SolverTest, AddVarAndSolve) {
  Args args("", "data/simple.nl");
  TestSolver solver("testsolver");
  solver.ReadProblem(args);
  Solution s;
  ProblemChanges changes(solver.problem());
  EXPECT_EQ(0, changes.num_vars());
  EXPECT_EQ(0, changes.num_cons());
  EXPECT_EQ(0, changes.num_objs());
  changes.AddVar(42, 42);
  EXPECT_EQ(1, changes.num_vars());
  EXPECT_EQ(0, changes.num_cons());
  EXPECT_EQ(0, changes.num_objs());
  solver.problem().Solve(CBC_PATH, s, &changes);
  EXPECT_EQ(3, s.num_vars());
  EXPECT_EQ(1, s.num_cons());
  EXPECT_EQ(2, s.value(0));
  EXPECT_NEAR(0, s.value(1), 1e-5);
  EXPECT_EQ(42, s.value(2));
  EXPECT_EQ(0.5, s.dual_value(0));
}

TEST(SolverTest, AddConAndSolve) {
  Args args("", "data/simple.nl");
  TestSolver solver("testsolver");
  solver.ReadProblem(args);
  Solution s;
  ProblemChanges changes(solver.problem());
  const double coefs[] = {1, 0};
  EXPECT_EQ(0, changes.num_vars());
  EXPECT_EQ(0, changes.num_cons());
  EXPECT_EQ(0, changes.num_objs());
  changes.AddCon(coefs, -Infinity, 1);
  EXPECT_EQ(0, changes.num_vars());
  EXPECT_EQ(1, changes.num_cons());
  EXPECT_EQ(0, changes.num_objs());
  solver.problem().Solve(CBC_PATH, s, &changes);
  EXPECT_EQ(2, s.num_vars());
  EXPECT_EQ(2, s.num_cons());
  EXPECT_EQ(1, s.value(0));
  EXPECT_EQ(0.5, s.value(1));
  EXPECT_EQ(0.5, s.dual_value(0));
  EXPECT_EQ(0.5, s.dual_value(1));
}

TEST(SolverTest, AddObjAndSolve) {
  Args args("", "data/noobj.nl");
  TestSolver solver("testsolver");
  solver.ReadProblem(args);
  Solution s;
  ProblemChanges changes(solver.problem());
  double coef = -1;
  int var = 0;
  EXPECT_EQ(0, changes.num_vars());
  EXPECT_EQ(0, changes.num_cons());
  EXPECT_EQ(0, changes.num_objs());
  changes.AddObj(ampl::MAX, 1, &coef, &var);
  EXPECT_EQ(0, changes.num_vars());
  EXPECT_EQ(0, changes.num_cons());
  EXPECT_EQ(1, changes.num_objs());
  solver.problem().Solve(CBC_PATH, s, &changes);
  EXPECT_EQ(1, s.num_vars());
  EXPECT_EQ(1, s.num_cons());
  EXPECT_EQ(0, s.value(0));
  EXPECT_EQ(-1, s.dual_value(0));
}

TEST(SolverTest, SolveIgnoreFunctions) {
  Args args("", "data/ssd.nl");
  TestSolver solver("testsolver");
  char amplfunc[] = "AMPLFUNC=../solvers/ssdsolver/ssd.dll";
  putenv(amplfunc);
  solver.ReadProblem(args);
  Solution s;
  solver.problem().Solve(CBC_PATH, s, 0, Problem::IGNORE_FUNCTIONS);
  EXPECT_EQ(42, s.value(0));
}
#endif

TEST(SolverTest, SolveWithUnknownSolver) {
  Args args("", "data/simple.nl");
  TestSolver solver("testsolver");
  solver.ReadProblem(args);
  Solution s;
  EXPECT_THROW(solver.problem().Solve("unknownsolver", s), ampl::Error);
}

TEST(SolverTest, WriteProblem) {
  Args args("", "data/simple.nl");
  TestSolver solver("testsolver");
  solver.ReadProblem(args);
  fmt::Writer writer;
  writer << solver.problem();
  EXPECT_EQ(
      "var x1 <= 2;\n"
      "var x2;\n"
      "maximize o: x1 + x2;\n"
      "s.t. c1: x1 + 2 * x2 <= 2;\n", writer.str());
}
