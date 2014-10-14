/*
 Problem tests.

 Copyright (C) 2013 AMPL Optimization Inc

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

#include <gmock/gmock.h>
#include "asl/aslbuilder.h"
#include "asl/problem.h"
#include "mp/nl.h"
#include "../util.h"
#include "stderr-redirect.h"

using mp::var::CONTINUOUS;
using mp::var::INTEGER;
using mp::LinearConExpr;
using mp::LinearObjExpr;
using mp::Problem;
using mp::ProblemChanges;
using mp::Solution;

namespace obj = mp::obj;
namespace suf = mp::suf;

#ifdef _WIN32
# define putenv _putenv
#endif

#ifndef MP_TEST_DATA_DIR
# define MP_TEST_DATA_DIR "../data"
#endif

TEST(SolutionTest, DefaultCtor) {
  Solution s;
  EXPECT_EQ(mp::sol::UNKNOWN, s.status());
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
  EXPECT_EQ(mp::sol::UNKNOWN, s.status());
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
  EXPECT_THROW(s.Read("nonexistent", 0, 0), mp::Error);
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
  namespace sol = mp::sol;
  const sol::Status STATES[] = {
      sol::SOLVED,
      sol::UNSOLVED,
      sol::INFEASIBLE,
      sol::UNBOUNDED,
      sol::LIMIT,
      sol::FAILURE
  };
  for (int i = 0,
         n = static_cast<int>(sizeof(STATES) / sizeof(*STATES)); i < n; ++i) {
    {
      int solve_code = i * 100;
      WriteFile("test.sol",
          fmt::format("test\n\n2\n2\nobjno 0 {}\n", solve_code));
      Solution s;
      s.Read("test", 1, 1);
      EXPECT_EQ(STATES[i], s.status());
      EXPECT_EQ(solve_code, s.solve_code());
    }
    {
      int solve_code = i * 100 + 99;
      WriteFile("test.sol",
          fmt::format("test\n\n2\n2\nobjno 0 {}\n", solve_code));
      Solution s;
      s.Read("test", 1, 1);
      EXPECT_EQ(STATES[i], s.status());
      EXPECT_EQ(solve_code, s.solve_code());
    }
  }
  const double CODES[] = {-5, -1, 600, 1000};
  for (std::size_t i = 0; i < sizeof(CODES) / sizeof(*CODES); ++i) {
    WriteFile("test.sol",
        fmt::format("test\n\n2\n2\nobjno 0 {}\n", CODES[i]));
    Solution s;
    s.Read("test", 1, 1);
    EXPECT_EQ(sol::UNKNOWN, s.status());
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

TEST(ProblemTest, EmptyProblem) {
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

TEST(ProblemTest, ProblemAccessors) {
  Problem p;
  p.Read(MP_TEST_DATA_DIR "/test");
  EXPECT_EQ(5, p.num_vars());
  EXPECT_EQ(19, p.num_objs());
  EXPECT_EQ(13, p.num_cons());
  EXPECT_EQ(2, p.num_integer_vars());
  EXPECT_EQ(3, p.num_continuous_vars());
  EXPECT_EQ(17, p.num_nonlinear_objs());
  EXPECT_EQ(11, p.num_nonlinear_cons());
  EXPECT_EQ(7, p.num_logical_cons());

  EXPECT_EQ(CONTINUOUS, p.var_type(0));
  EXPECT_EQ(CONTINUOUS, p.var_type(1));
  EXPECT_EQ(CONTINUOUS, p.var_type(2));
  EXPECT_EQ(INTEGER, p.var_type(3));
  EXPECT_EQ(INTEGER, p.var_type(4));

  EXPECT_EQ(11, p.var_lb(0));
  EXPECT_EQ(15, p.var_lb(p.num_vars() - 1));
  EXPECT_EQ(21, p.var_ub(0));
  EXPECT_EQ(25, p.var_ub(p.num_vars() - 1));

  EXPECT_EQ(101, p.con_lb(0));
  EXPECT_EQ(113, p.con_lb(p.num_cons() - 1));
  EXPECT_EQ(201, p.con_ub(0));
  EXPECT_EQ(213, p.con_ub(p.num_cons() - 1));

  EXPECT_EQ(obj::MIN, p.obj_type(0));
  EXPECT_EQ(obj::MAX, p.obj_type(p.num_objs() - 1));

  {
    LinearObjExpr expr = p.linear_obj_expr(0);
    EXPECT_EQ(31, expr.begin()->coef());
    EXPECT_EQ(0, expr.begin()->var_index());
    EXPECT_EQ(5, std::distance(expr.begin(), expr.end()));
    expr = p.linear_obj_expr(p.num_objs() - 1);
    EXPECT_EQ(52, expr.begin()->coef());
    EXPECT_EQ(3, expr.begin()->var_index());
  }

  {
    LinearConExpr expr = p.linear_con_expr(0);
    EXPECT_EQ(61, expr.begin()->coef());
    EXPECT_EQ(0, expr.begin()->var_index());
    EXPECT_EQ(5, std::distance(expr.begin(), expr.end()));
    expr = p.linear_con_expr(p.num_cons() - 1);
    EXPECT_EQ(82, expr.begin()->coef());
    EXPECT_EQ(2, expr.begin()->var_index());
  }

  EXPECT_EQ(mp::expr::SIN, p.nonlinear_obj_expr(0).kind());
  EXPECT_EQ(mp::expr::COS,
            p.nonlinear_obj_expr(p.num_nonlinear_objs() - 1).kind());

  EXPECT_EQ(mp::expr::LOG, p.nonlinear_con_expr(0).kind());
  EXPECT_EQ(mp::expr::EXP,
            p.nonlinear_con_expr(p.num_nonlinear_cons() - 1).kind());

  EXPECT_EQ(mp::expr::NE, p.logical_con_expr(0).kind());
  EXPECT_EQ(mp::expr::AND, p.logical_con_expr(p.num_logical_cons() - 1).kind());

  EXPECT_EQ(-1, p.solve_code());
  p.set_solve_code(42);
  EXPECT_EQ(42, p.solve_code());
}

TEST(ProblemTest, VarType) {
  Problem p;
  p.AddVar(0, 0, CONTINUOUS);
  p.AddVar(0, 0, INTEGER);
  p.AddVar(0, 0, INTEGER);
  EXPECT_EQ(CONTINUOUS, p.var_type(0));
  EXPECT_EQ(INTEGER, p.var_type(1));
  EXPECT_EQ(INTEGER, p.var_type(2));
}

#ifndef NDEBUG
TEST(ProblemTest, BoundChecks) {
  Problem p;
  p.Read(MP_TEST_DATA_DIR "/test");

  EXPECT_DEATH(p.var_type(-1), "Assertion");
  EXPECT_DEATH(p.var_type(p.num_vars()), "Assertion");

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

#ifdef HAVE_ILOGCP
static const std::string SOLVER_PATH = GetExecutableDir() + "/ilogcp";

TEST(ProblemTest, Solve) {
  Problem p;
  p.Read(MP_TEST_DATA_DIR "/simple");
  Solution s;
  p.Solve(SOLVER_PATH, s);
  EXPECT_EQ(2, s.num_vars());
  EXPECT_EQ(1, s.num_cons());
  EXPECT_EQ(2, s.value(0));
  EXPECT_NEAR(0, s.value(1), 1e-5);
  EXPECT_EQ(1, s.dual_value(0));
}

TEST(ProblemChangesTest, AddVarAndSolve) {
  Problem p;
  p.Read(MP_TEST_DATA_DIR "/simple");
  Solution s;
  ProblemChanges changes(p);
  EXPECT_EQ(0, changes.num_vars());
  EXPECT_EQ(0, changes.num_cons());
  EXPECT_EQ(0, changes.num_objs());
  changes.AddVar(42, 42);
  EXPECT_EQ(1, changes.num_vars());
  EXPECT_EQ(0, changes.num_cons());
  EXPECT_EQ(0, changes.num_objs());
  p.Solve(SOLVER_PATH, s, &changes);
  EXPECT_EQ(3, s.num_vars());
  EXPECT_EQ(1, s.num_cons());
  EXPECT_EQ(2, s.value(0));
  EXPECT_NEAR(0, s.value(1), 1e-5);
  EXPECT_EQ(42, s.value(2));
  EXPECT_EQ(1, s.dual_value(0));
}

TEST(ProblemChangesTest, AddConAndSolve) {
  Problem p;
  p.Read(MP_TEST_DATA_DIR "/simple");
  Solution s;
  ProblemChanges changes(p);
  const double coefs[] = {1, 0};
  EXPECT_EQ(0, changes.num_vars());
  EXPECT_EQ(0, changes.num_cons());
  EXPECT_EQ(0, changes.num_objs());
  changes.AddCon(coefs, -Infinity, 1);
  EXPECT_EQ(0, changes.num_vars());
  EXPECT_EQ(1, changes.num_cons());
  EXPECT_EQ(0, changes.num_objs());
  p.Solve(SOLVER_PATH, s, &changes);
  EXPECT_EQ(2, s.num_vars());
  EXPECT_EQ(2, s.num_cons());
  EXPECT_EQ(1, s.value(0));
  EXPECT_EQ(0.5, s.value(1));
  EXPECT_EQ(0.5, s.dual_value(0));
  EXPECT_EQ(0.5, s.dual_value(1));
}

TEST(ProblemChangesTest, AddObjAndSolve) {
  Problem p;
  p.Read(MP_TEST_DATA_DIR "/noobj");
  Solution s;
  ProblemChanges changes(p);
  double coef = -1;
  int var = 0;
  EXPECT_EQ(0, changes.num_vars());
  EXPECT_EQ(0, changes.num_cons());
  EXPECT_EQ(0, changes.num_objs());
  changes.AddObj(obj::MAX, 1, &coef, &var);
  EXPECT_EQ(0, changes.num_vars());
  EXPECT_EQ(0, changes.num_cons());
  EXPECT_EQ(1, changes.num_objs());
  p.Solve(SOLVER_PATH, s, &changes);
  EXPECT_EQ(1, s.num_vars());
  EXPECT_EQ(1, s.num_cons());
  EXPECT_EQ(0, s.value(0));
  EXPECT_EQ(-1, s.dual_value(0));
}

TEST(ProblemChangesTest, CopyConstructorCon) {
  Problem p;
  p.Read(MP_TEST_DATA_DIR "/simple");

  ProblemChanges changes(p);
  for (int i = 0; i < 10; ++i) {
    ProblemChanges next(changes);
    EXPECT_EQ(next.num_cons(), changes.num_cons());
    const double coefs[] = {1, 0};
    next.AddCon(coefs, -Infinity, Infinity);
    changes = next;
    EXPECT_EQ(next.num_cons(), changes.num_cons());
  }
  EXPECT_EQ(10, changes.num_cons());
  Solution s;
  p.Solve(SOLVER_PATH, s, &changes);
  EXPECT_EQ(2, s.num_vars());
  EXPECT_EQ(11, s.num_cons());
  EXPECT_EQ(2, s.value(0));
  EXPECT_NEAR(0, s.value(1), 1e-5);
  EXPECT_EQ(1, s.dual_value(0));
}

TEST(ProblemTest, SolveIgnoreFunctions) {
  char amplfunc[] = "AMPLFUNC=../../solvers/ssdsolver/ssd.dll";
  putenv(amplfunc);
  Problem p;
  p.Read(MP_TEST_DATA_DIR "/ssd");
  Solution s;
  p.Solve(SOLVER_PATH, s, 0, Problem::IGNORE_FUNCTIONS);
  EXPECT_EQ(42, s.value(0));
}
#endif

TEST(ProblemTest, SolveWithUnknownSolver) {
  Problem p;
  p.Read(MP_TEST_DATA_DIR "/simple");
  Solution s;
  EXPECT_THROW(p.Solve("unknownsolver", s), mp::Error);
}

TEST(ProblemTest, Write) {
  Problem p;
  p.Read(MP_TEST_DATA_DIR "/simple");
  fmt::MemoryWriter writer;
  writer << p;
  EXPECT_EQ(
      "var x1 >= 0;\n"
      "var x2 >= 0;\n"
      "maximize o: x1 + x2;\n"
      "s.t. c1: x1 + 2 * x2 <= 2;\n", writer.str());
}

TEST(ProblemTest, WriteVarBounds) {
  Problem p;
  p.AddVar(42, 42);
  fmt::MemoryWriter writer;
  writer << p;
  EXPECT_EQ("var x1 = 42;\n", writer.str());
}

TEST(ProblemTest, AddVar) {
  Problem p;
  EXPECT_EQ(0, p.num_vars());
  p.AddVar(111, 222);
  EXPECT_EQ(1, p.num_vars());
  EXPECT_EQ(0, p.num_integer_vars());
  EXPECT_EQ(1, p.num_continuous_vars());
  EXPECT_EQ(CONTINUOUS, p.var_type(0));
  EXPECT_EQ(111, p.var_lb(0));
  EXPECT_EQ(222, p.var_ub(0));

  p.AddVar(333, 444, INTEGER);
  EXPECT_EQ(2, p.num_vars());
  EXPECT_EQ(1, p.num_integer_vars());
  EXPECT_EQ(1, p.num_continuous_vars());
  EXPECT_EQ(INTEGER, p.var_type(1));
  EXPECT_EQ(333, p.var_lb(1));
  EXPECT_EQ(444, p.var_ub(1));

  p.Read(MP_TEST_DATA_DIR "/simple");
  EXPECT_THROW(p.AddVar(0, 0), mp::Error);
}

class TestASLBuilder : public mp::internal::ASLBuilder {
 private:
  static mp::ProblemInfo MakeProblemInfo() {
    mp::ProblemInfo info = mp::ProblemInfo();
    info.num_vars = info.num_objs = 1;
    return info;
  }

 public:
  explicit TestASLBuilder(mp::ProblemInfo info = MakeProblemInfo()) {
    set_flags(mp::internal::ASL_STANDARD_OPCODES);
    SetInfo(info);
  }
};

TEST(ProblemTest, AddCon) {
  Problem p;
  p.AddVar(0, 0);
  EXPECT_EQ(0, p.num_logical_cons());
  TestASLBuilder builder;
  mp::LogicalExpr expr = builder.MakeRelational(
        mp::expr::EQ, builder.MakeVariable(0), builder.MakeNumericConstant(0));
  p.AddCon(expr);
  EXPECT_EQ(0, p.num_cons());
  EXPECT_EQ(1, p.num_logical_cons());
  EXPECT_EQ(expr, p.logical_con_expr(0));

  p.Read(MP_TEST_DATA_DIR "/test");
  EXPECT_THROW(p.AddCon(expr), mp::Error);
}

TEST(ProblemTest, AddObj) {
  Problem p;
  p.AddVar(0, 0);
  EXPECT_EQ(0, p.num_objs());
  TestASLBuilder builder;
  mp::NumericExpr expr = builder.MakeBinary(
        mp::expr::ADD, builder.MakeVariable(0), builder.MakeNumericConstant(1));
  p.AddObj(obj::MAX, expr);
  EXPECT_EQ(1, p.num_objs());
  EXPECT_EQ(obj::MAX, p.obj_type(0));
  EXPECT_EQ(expr, p.nonlinear_obj_expr(0));

  p.Read(MP_TEST_DATA_DIR "/simple");
  EXPECT_THROW(p.AddObj(obj::MAX, expr), mp::Error);
}

TEST(ProblemTest, ReadFunctionWithoutLibrary) {
  Problem p;
  // It shouldn't be an error to have a function without an implementation
  // (provided by a function library) because some functions are not evaluated
  // but translated into the solver representation.
  p.Read(MP_TEST_DATA_DIR "/element");
  EXPECT_EQ(1, p.num_objs());
}

TEST(ProblemTest, Proxy) {
  mp::ProblemInfo info = mp::ProblemInfo();
  info.num_vars = 42;
  info.num_objs = 1;
  TestASLBuilder builder(info);
  builder.MakeVariable(0);
  Problem p(builder.GetProblem());
  EXPECT_EQ(info.num_vars, p.num_vars());
}

typedef ::testing::TestWithParam<int> SuffixTest;

TEST_P(SuffixTest, FindSuffix) {
  TestASLBuilder builder;
  builder.set_flags(ASL_keep_all_suffixes);
  int kind = GetParam();
  builder.AddSuffix(kind, 1, "foo");
  builder.AddSuffix(kind, 2, "bar");
  Problem p(builder.GetProblem());
  mp::ASLSuffixPtr suffix = p.suffixes(kind).Find("foo");
  EXPECT_TRUE(suffix);
  EXPECT_STREQ("foo", suffix->name());
  EXPECT_EQ(kind, suffix->kind() & suf::MASK);
  EXPECT_STREQ("bar", p.suffixes(kind).Find("bar")->name());
  int kinds[] = {suf::VAR, suf::CON, suf::OBJ, suf::PROBLEM};
  for (std::size_t i = 0, n = sizeof(kinds) / sizeof(*kinds); i != n; ++i) {
    if (kinds[i] != kind)
      EXPECT_FALSE(p.suffixes(kinds[i]).Find("foo"));
  }
}

TEST_P(SuffixTest, SuffixView) {
  TestASLBuilder builder;
  builder.set_flags(ASL_keep_all_suffixes);
  int kind = GetParam();
  builder.AddSuffix(kind, 1, "foo");
  builder.AddSuffix(kind, 2, "bar");
  Problem p(builder.GetProblem());
  mp::SuffixView view = p.suffixes(kind);
  mp::SuffixView::iterator it = view.begin();
  ASSERT_NE(it, view.end());
  // Suffixes are returned in the reverse order of insertion.
  EXPECT_STREQ("bar", (*it).name());
  EXPECT_STREQ("foo", (*++it).name());
  EXPECT_STREQ("foo", it->name());
  EXPECT_STREQ("foo", (*it++).name());
  EXPECT_EQ(it, view.end());
}

TEST_P(SuffixTest, EmptySuffixView) {
  Problem p;
  mp::SuffixView view = p.suffixes(GetParam());
  EXPECT_EQ(view.begin(), view.end());
}

struct MockValueVisitor {
  MOCK_METHOD2(Visit, void (int index, int value));
};

TEST_P(SuffixTest, VisitValues) {
  mp::ProblemInfo info = mp::ProblemInfo();
  info.num_vars = info.num_objs = info.num_algebraic_cons = 3;
  TestASLBuilder builder(info);
  builder.set_flags(ASL_keep_all_suffixes);
  int kind = GetParam();
  mp::internal::ASLBuilder::SuffixHandler handler =
      builder.AddSuffix(kind, 3, "foo");
  handler.SetValue(0, 11);
  if (kind != mp::suf::PROBLEM) {
    handler.SetValue(1, 22);
    handler.SetValue(2, 33);
  }
  Problem p(builder.GetProblem());
  testing::StrictMock<MockValueVisitor> visitor;
  mp::ASLSuffixPtr suffix = p.suffixes(kind).Find("foo");
  testing::InSequence dummy;
  EXPECT_CALL(visitor, Visit(0, 11));
  if (kind != mp::suf::PROBLEM) {
    EXPECT_CALL(visitor, Visit(1, 22));
    EXPECT_CALL(visitor, Visit(2, 33));
  }
  suffix->VisitValues(visitor);
}

INSTANTIATE_TEST_CASE_P(, SuffixTest, ::testing::Values(
                          int(suf::VAR), int(suf::CON), int(suf::OBJ),
                          int(suf::PROBLEM)));
