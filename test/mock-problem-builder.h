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

template <int I>
struct BasicTestExpr {
  DEFINE_ID(BasicTestExpr);
};

typedef BasicTestExpr<0> TestExpr;
typedef BasicTestExpr<1> TestNumericExpr;
typedef BasicTestExpr<2> TestLogicalExpr;

struct TestVariable : TestNumericExpr { DEFINE_ID(TestVariable); };
struct TestCountExpr : TestNumericExpr { DEFINE_ID(TestCountExpr); };

template <int I>
struct TestLinearExprBuilder {
  void AddTerm(int, double) {}
  DEFINE_ID(TestLinearExprBuilder);
};

typedef TestLinearExprBuilder<0> TestLinearObjBuilder;
typedef TestLinearExprBuilder<1> TestLinearConBuilder;
typedef TestLinearExprBuilder<2> TestLinearVarBuilder;

struct TestColumnSizeHandler {
  void Add(int) {}
  DEFINE_ID(TestColumnSizeHandler);
};

struct TestSuffixHandler {
  void SetValue(int, int) {}
  void SetValue(int, double) {}
  DEFINE_ID(TestSuffixHandler);
};

template <int I, typename Arg = TestNumericExpr>
struct TestArgHandler {
  void AddArg(Arg) {}
  DEFINE_ID(TestArgHandler);
};

typedef TestArgHandler<0> TestNumericArgHandler;
typedef TestArgHandler<1, TestLogicalExpr> TestLogicalArgHandler;
typedef TestArgHandler<2> TestAllDiffArgHandler;

struct TestCallArgHandler {
  void AddArg(TestNumericExpr) {}
  void AddArg(TestExpr) {}
  DEFINE_ID(TestCallArgHandler);
};

struct TestPLTermHandler {
  void AddSlope(double) {}
  void AddBreakpoint(double) {}
  DEFINE_ID(TestPLTermHandler);
};

// A mock problem builder.
class MockProblemBuilder {
 public:
  typedef TestExpr Expr;
  typedef TestNumericExpr NumericExpr;
  typedef TestLogicalExpr LogicalExpr;
  typedef TestCountExpr CountExpr;
  typedef TestVariable Variable;

  MockProblemBuilder() {}

  // Constructs a MockProblemBuilder object and stores a pointer to it
  // in ``builder``.
  explicit MockProblemBuilder(MockProblemBuilder *&builder) { builder = this; }

  typedef mp::Suffix *SuffixPtr;
  typedef mp::SuffixSet SuffixSet;

  MOCK_METHOD1(suffixes, SuffixSet &(int));

  MOCK_METHOD0(num_vars, int ());
  MOCK_METHOD0(num_cons, int ());

  MOCK_METHOD1(SetInfo, void (const mp::ProblemInfo &info));
  MOCK_METHOD0(EndBuild, void ());

  MOCK_METHOD3(SetObj, void (int index, mp::obj::Type type, NumericExpr expr));
  MOCK_METHOD2(SetCon, void (int index, NumericExpr expr));
  MOCK_METHOD2(SetLogicalCon, void (int index, LogicalExpr expr));
  MOCK_METHOD3(SetCommonExpr, void (int index, NumericExpr expr, int position));
  MOCK_METHOD3(SetComplement, void (int con_index, int var_index, int flags));

  typedef TestLinearObjBuilder LinearObjBuilder;

  MOCK_METHOD2(GetLinearObjBuilder,
               LinearObjBuilder (int obj_index, int num_linear_terms));

  typedef TestLinearConBuilder LinearConBuilder;

  MOCK_METHOD2(GetLinearConBuilder,
               LinearConBuilder (int con_index, int num_linear_terms));

  typedef TestLinearVarBuilder LinearVarBuilder;

  MOCK_METHOD2(GetLinearVarBuilder,
               LinearVarBuilder (int var_index, int num_linear_terms));

  MOCK_METHOD3(SetVarBounds, void (int index, double lb, double ub));
  MOCK_METHOD3(SetConBounds, void (int index, double lb, double ub));

  MOCK_METHOD2(SetInitialValue, void (int var_index, double value));
  MOCK_METHOD2(SetInitialDualValue, void (int con_index, double value));

  typedef TestColumnSizeHandler ColumnSizeHandler;

  MOCK_METHOD0(GetColumnSizeHandler, ColumnSizeHandler ());

  MOCK_METHOD4(SetFunction, void (int index, fmt::StringRef name,
                                  int num_args, mp::func::Type type));

  typedef TestSuffixHandler SuffixHandler;

  MOCK_METHOD3(AddSuffix,
               SuffixHandler (int kind, int num_values, fmt::StringRef name));

  typedef TestNumericArgHandler NumericArgHandler;
  typedef TestLogicalArgHandler LogicalArgHandler;
  typedef TestCallArgHandler CallArgHandler;

  MOCK_METHOD1(MakeNumericConstant, NumericExpr (double value));
  MOCK_METHOD1(MakeVariable, Variable (int var_index));
  MOCK_METHOD2(MakeUnary, NumericExpr (mp::expr::Kind kind, NumericExpr arg));
  MOCK_METHOD3(MakeBinary,
               NumericExpr (mp::expr::Kind kind,
                            NumericExpr lhs, NumericExpr rhs));

  MOCK_METHOD3(MakeIf,
               NumericExpr (LogicalExpr condition,
                            NumericExpr true_expr, NumericExpr false_expr));

  typedef TestPLTermHandler PLTermHandler;

  MOCK_METHOD1(BeginPLTerm, PLTermHandler (int num_breakpoints));
  MOCK_METHOD2(EndPLTerm, NumericExpr (PLTermHandler handler, NumericExpr arg));

  MOCK_METHOD2(BeginCall, CallArgHandler (int func_index, int num_args));
  MOCK_METHOD1(EndCall, NumericExpr (CallArgHandler handler));

  MOCK_METHOD2(BeginVarArg,
               NumericArgHandler (mp::expr::Kind kind, int num_args));
  MOCK_METHOD1(EndVarArg, NumericExpr (NumericArgHandler handler));

  MOCK_METHOD1(BeginSum, NumericArgHandler (int num_args));
  MOCK_METHOD1(EndSum, NumericExpr (NumericArgHandler handler));

  MOCK_METHOD1(BeginCount, LogicalArgHandler (int num_args));
  MOCK_METHOD1(EndCount, CountExpr (LogicalArgHandler handler));

  MOCK_METHOD2(BeginNumberOf,
               NumericArgHandler (int num_args, NumericExpr value));
  MOCK_METHOD1(EndNumberOf, NumericExpr (NumericArgHandler handler));

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
                            LogicalExpr true_expr, LogicalExpr false_expr));

  MOCK_METHOD2(BeginIteratedLogical,
               LogicalArgHandler (mp::expr::Kind kind, int num_args));
  MOCK_METHOD1(EndIteratedLogical, LogicalExpr (LogicalArgHandler handler));

  typedef TestAllDiffArgHandler AllDiffArgHandler;

  MOCK_METHOD1(BeginAllDiff, AllDiffArgHandler (int num_args));
  MOCK_METHOD1(EndAllDiff, LogicalExpr (AllDiffArgHandler handler));

  MOCK_METHOD1(MakeStringLiteral, Expr (fmt::StringRef value));
};

#endif  // MP_MOCK_PROBLEM_BUILDER_H_
