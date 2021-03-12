/*
 Mock problem builder

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

#ifndef MP_MOCK_PROBLEM_BUILDER_H_
#define MP_MOCK_PROBLEM_BUILDER_H_

#include "gmock/gmock.h"
#include "mp/problem-builder.h"

enum IDType { NULL_ID = 0, ID = 42, ID2, ID3, ID4 };

// Defines a field of type IDType, a construtor to set it and operator==
// that compares objects by their ids. The constructor takes an argument of
// type IDType, a distinct type used to make sure that the constructor is not
// called with a wrong argument in the tested code.
#define DEFINE_ID(Class) \
 public: \
  explicit Class(IDType id = NULL_ID) : id_(id) {} \
  friend bool operator==(Class lhs, Class rhs) { return lhs.id_ == rhs.id_; } \
 private: \
  IDType id_

class TestFunction {
 private:
  // Safe bool type.
  typedef void (TestFunction::*SafeBool)() const;

  // A member function representing the true value of SafeBool.
  void True() const {}

 public:
  // Returns a value convertible to bool that can be used in conditions but
  // not in comparisons and evaluates to "true" if this function is not null
  // and "false" otherwise.
  // Example:
  //   if (f) {
  //     // Do something if f is not null.
  //   }
  operator SafeBool() const { return &TestFunction::True; }

  DEFINE_ID(TestFunction);
};

template <int I>
struct BasicTestExpr {
  DEFINE_ID(BasicTestExpr);
};

typedef BasicTestExpr<0> TestExpr;
typedef BasicTestExpr<1> TestNumericExpr;
typedef BasicTestExpr<2> TestLogicalExpr;

struct TestReference : TestNumericExpr { DEFINE_ID(TestReference); };
struct TestCountExpr : TestNumericExpr { DEFINE_ID(TestCountExpr); };

template <int I>
struct BasicTestLinearExprBuilder {
  void AddTerm(int, double) {}
  DEFINE_ID(BasicTestLinearExprBuilder);
};

typedef BasicTestLinearExprBuilder<0> TestLinearObjBuilder;
typedef BasicTestLinearExprBuilder<1> TestLinearConBuilder;
typedef BasicTestLinearExprBuilder<1> TestLinearExprBuilder;

struct TestColumnSizeHandler {
  void Add(int) {}
  DEFINE_ID(TestColumnSizeHandler);
};

template <int I>
struct TestSuffixHandler {
  void SetValue(int, int) {}
  void SetValue(int, double) {}
  DEFINE_ID(TestSuffixHandler);
};

template <int I, typename Arg = TestNumericExpr>
struct TestExprBuilder {
  void AddArg(Arg) {}
  DEFINE_ID(TestExprBuilder);
};

typedef TestExprBuilder<0> TestNumericExprBuilder;
typedef TestExprBuilder<1> TestIteratedExprBuilder;
typedef TestExprBuilder<2> TestNumberOfExprBuilder;
typedef TestExprBuilder<0, TestExpr> TestSymbolicNumberOfExprBuilder;
typedef TestExprBuilder<3, TestLogicalExpr> TestCountExprBuilder;
typedef TestExprBuilder<4, TestLogicalExpr> TestIteratedLogicalExprBuilder;
typedef TestExprBuilder<5> TestPairwiseExprBuilder;

struct TestCallExprBuilder {
  void AddArg(TestNumericExpr) {}
  void AddArg(TestExpr) {}
  DEFINE_ID(TestCallExprBuilder);
};

struct TestPLTermBuilder {
  void AddSlope(double) {}
  void AddBreakpoint(double) {}
  DEFINE_ID(TestPLTermBuilder);
};

// A mock problem builder.
class MockProblemBuilder {
 public:
  typedef TestFunction Function;
  typedef TestExpr Expr;
  typedef TestNumericExpr NumericExpr;
  typedef TestLogicalExpr LogicalExpr;
  typedef TestCountExpr CountExpr;
  typedef TestReference Reference;

  MockProblemBuilder() {}

  // Constructs a MockProblemBuilder object and stores a pointer to it
  // in ``solver``.
  template <typename Solver>
  explicit MockProblemBuilder(Solver &s, fmt::StringRef) {
    s.builder = this;
  }

  typedef mp::MutSuffix Suffix;
  typedef mp::MutIntSuffix IntSuffix;
  typedef mp::SuffixSet SuffixSet;

  MOCK_METHOD1(suffixes, SuffixSet &(int));

  MOCK_METHOD0(num_vars, int ());
  MOCK_METHOD0(num_algebraic_cons, int ());

  MOCK_METHOD1(SetInfo, void (const mp::ProblemInfo &info));
  MOCK_METHOD0(EndBuild, void ());

  MOCK_METHOD3(AddVar, void (double lb, double ub, mp::var::Type));
  MOCK_METHOD2(AddVars, void (int num_vars, mp::var::Type));

  typedef TestLinearObjBuilder LinearObjBuilder;

  MOCK_METHOD1(AddObj, void (mp::obj::Type type));
  MOCK_METHOD1(AddObjs, void (int num_objs));

  typedef TestLinearConBuilder LinearConBuilder;

  MOCK_METHOD2(AddCon, void (double lb, double ub));
  MOCK_METHOD1(AddCon, void (LogicalExpr expr));
  MOCK_METHOD1(AddAlgebraicCons, void (int num_cons));
  MOCK_METHOD1(AddLogicalCons, void (int num_cons));

  typedef TestLinearExprBuilder LinearExprBuilder;

  class MutCommonExpr {
   public:
    MutCommonExpr() {}
    MutCommonExpr(const MutCommonExpr &) {}

    MOCK_CONST_METHOD1(set_linear_expr,
                       LinearExprBuilder (int num_linear_terms));
    MOCK_CONST_METHOD1(set_nonlinear_expr, void (NumericExpr expr));
    MOCK_CONST_METHOD1(set_position, void (int position));
  };

  MOCK_METHOD1(common_expr, MutCommonExpr &(int var_index));

  MOCK_METHOD1(AddCommonExpr, MutCommonExpr &(NumericExpr expr));
  MOCK_METHOD1(AddCommonExprs, void (int num_exprs));

  MOCK_METHOD3(SetComplementarity,
               void (int con_index, int var_index, mp::ComplInfo info));

  class MutVariable {
   public:
    MutVariable() {}
    MutVariable(const MutVariable &) {}
    MOCK_CONST_METHOD1(set_lb, void (double lb));
    MOCK_CONST_METHOD1(set_ub, void (double ub));
    MOCK_CONST_METHOD1(set_value, void (double value));
  };

  MOCK_METHOD1(var, MutVariable &(int var_index));

  class Objective {
   public:
    MOCK_CONST_METHOD1(set_type, void (mp::obj::Type type));
    MOCK_CONST_METHOD1(set_linear_expr,
                       LinearObjBuilder (int num_linear_terms));
    MOCK_CONST_METHOD1(set_nonlinear_expr, void (NumericExpr expr));
  };

  MOCK_METHOD1(obj, Objective &(int obj_index));

  class AlgebraicCon {
   public:
    MOCK_CONST_METHOD1(set_lb, void (double value));
    MOCK_CONST_METHOD1(set_ub, void (double value));
    MOCK_CONST_METHOD1(set_dual, void (double value));
    MOCK_CONST_METHOD1(set_linear_expr,
                       LinearConBuilder (int num_linear_terms));
    MOCK_CONST_METHOD1(set_nonlinear_expr, void (NumericExpr expr));
  };

  MOCK_METHOD1(algebraic_con, AlgebraicCon &(int index));

  class LogicalCon {
   public:
    MOCK_CONST_METHOD1(set_expr, void (LogicalExpr expr));
  };

  MOCK_METHOD1(logical_con, LogicalCon &(int index));

  typedef TestColumnSizeHandler ColumnSizeHandler;

  MOCK_METHOD0(GetColumnSizeHandler, ColumnSizeHandler ());

  MOCK_METHOD1(AddFunctions, void (int num_funcs));

  MOCK_METHOD3(AddFunction, Function (fmt::StringRef name,
                                      int num_args, mp::func::Type type));

  MOCK_METHOD4(DefineFunction, Function (int index, fmt::StringRef name,
                                         int num_args, mp::func::Type type));

  MOCK_METHOD1(function, Function (int index));

  typedef TestSuffixHandler<0> IntSuffixHandler;

  MOCK_METHOD3(AddIntSuffix,
               IntSuffixHandler (fmt::StringRef name,
                                 int kind, int num_values));

  typedef TestSuffixHandler<1> DblSuffixHandler;

  MOCK_METHOD3(AddDblSuffix,
               DblSuffixHandler (fmt::StringRef name,
                                 mp::suf::Kind kind, int num_values));


  typedef TestNumericExprBuilder NumericExprBuilder;

  MOCK_METHOD1(MakeNumericConstant, NumericExpr (double value));
  MOCK_METHOD1(MakeVariable, Reference (int var_index));
  MOCK_METHOD1(MakeCommonExpr, Reference (int index));
  MOCK_METHOD2(MakeUnary, NumericExpr (mp::expr::Kind kind, NumericExpr arg));
  MOCK_METHOD3(MakeBinary,
               NumericExpr (mp::expr::Kind kind,
                            NumericExpr lhs, NumericExpr rhs));

  MOCK_METHOD3(MakeIf,
               NumericExpr (LogicalExpr condition,
                            NumericExpr then_expr, NumericExpr else_expr));

  typedef TestPLTermBuilder PLTermBuilder;

  MOCK_METHOD1(BeginPLTerm, PLTermBuilder (int num_breakpoints));
  MOCK_METHOD2(EndPLTerm, NumericExpr (PLTermBuilder builder, NumericExpr arg));

  typedef TestCallExprBuilder CallExprBuilder;

  MOCK_METHOD2(BeginCall, CallExprBuilder (Function func, int num_args));
  MOCK_METHOD1(EndCall, NumericExpr (CallExprBuilder builder));

  typedef TestIteratedExprBuilder IteratedExprBuilder;

  MOCK_METHOD2(BeginIterated,
               IteratedExprBuilder (mp::expr::Kind kind, int num_args));
  MOCK_METHOD1(EndIterated, NumericExpr (IteratedExprBuilder builder));

  MOCK_METHOD1(BeginSum, IteratedExprBuilder (int num_args));
  MOCK_METHOD1(EndSum, NumericExpr (IteratedExprBuilder builder));

  typedef TestCountExprBuilder CountExprBuilder;

  MOCK_METHOD1(BeginCount, CountExprBuilder (int num_args));
  MOCK_METHOD1(EndCount, CountExpr (CountExprBuilder builder));

  typedef TestNumberOfExprBuilder NumberOfExprBuilder;

  MOCK_METHOD2(BeginNumberOf,
               NumberOfExprBuilder (int num_args, NumericExpr arg0));
  MOCK_METHOD1(EndNumberOf, NumericExpr (NumberOfExprBuilder builder));

  typedef TestSymbolicNumberOfExprBuilder SymbolicNumberOfExprBuilder;

  MOCK_METHOD2(BeginSymbolicNumberOf,
               SymbolicNumberOfExprBuilder (int num_args, Expr arg0));
  MOCK_METHOD1(EndSymbolicNumberOf,
               NumericExpr (SymbolicNumberOfExprBuilder builder));

  MOCK_METHOD1(MakeLogicalConstant, LogicalExpr (bool value));
  MOCK_METHOD1(MakeNot, LogicalExpr (LogicalExpr arg));
  MOCK_METHOD3(MakeBinaryLogical,
               LogicalExpr (mp::expr::Kind kind,
                            LogicalExpr lhs, LogicalExpr rhs));
  MOCK_METHOD3(MakeRelational,
               LogicalExpr (mp::expr::Kind kind,
                            NumericExpr lhs, NumericExpr rhs));
  MOCK_METHOD3(MakeLogicalCount,
               LogicalExpr (mp::expr::Kind kind,
                            NumericExpr lhs, CountExpr rhs));
  MOCK_METHOD3(MakeImplication,
               LogicalExpr (LogicalExpr condition,
                            LogicalExpr then_expr, LogicalExpr else_expr));

  typedef TestIteratedLogicalExprBuilder IteratedLogicalExprBuilder;

  MOCK_METHOD2(BeginIteratedLogical,
               IteratedLogicalExprBuilder (mp::expr::Kind kind, int num_args));
  MOCK_METHOD1(EndIteratedLogical,
               LogicalExpr (IteratedLogicalExprBuilder builder));

  typedef TestPairwiseExprBuilder PairwiseExprBuilder;

  MOCK_METHOD2(BeginPairwise,
               PairwiseExprBuilder (mp::expr::Kind kind, int num_args));
  MOCK_METHOD1(EndPairwise, LogicalExpr (PairwiseExprBuilder builder));

  MOCK_METHOD1(MakeStringLiteral, Expr (fmt::StringRef value));

  MOCK_METHOD3(MakeSymbolicIf,
               Expr (LogicalExpr condition, Expr then_expr, Expr else_expr));

  MOCK_METHOD0(problem, MockProblemBuilder &());
};

#endif  // MP_MOCK_PROBLEM_BUILDER_H_
