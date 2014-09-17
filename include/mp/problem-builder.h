/*
 A minimal implementation of the ProblemBuilder concept.

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

#ifndef MP_PROBLEM_BUILDER_H_
#define MP_PROBLEM_BUILDER_H_

#include <memory>

#include "mp/error.h"
#include "mp/problem-base.h"

#define MP_UNUSED(x) (void)(x)

namespace mp {

// A minimal implementation of the ProblemBuilder concept.
template <typename Impl, typename ExprT>
class ProblemBuilder {
 private:
  struct LinearExprBuilder {
    void AddTerm(int var_index, double coef) {
      MP_UNUSED(var_index); MP_UNUSED(coef);
    }
  };

  struct ArgHandler {
    void AddArg(ExprT arg) { MP_UNUSED(arg); }
  };

 public:
  typedef ExprT Expr;
  typedef Expr NumericExpr;
  typedef Expr LogicalExpr;
  typedef Expr CountExpr;
  typedef Expr Variable;

  void ReportUnhandledConstruct(fmt::StringRef name) {
    throw Error("unsupported: {}", name);
  }

  void SetInfo(const ProblemInfo &) {}

  // Sets an objective type and expression.
  // index: Index of an objective; 0 <= index < num_objs.
  void SetObj(int index, obj::Type type, NumericExpr expr) {
    MP_UNUSED(index); MP_UNUSED(type); MP_UNUSED(expr);
    MP_DISPATCH(ReportUnhandledConstruct("objective"));
  }

  // Sets an algebraic constraint expression.
  // index: Index of an algebraic contraint; 0 <= index < num_algebraic_cons.
  void SetCon(int index, NumericExpr expr) {
    MP_UNUSED(index); MP_UNUSED(expr);
    MP_DISPATCH(ReportUnhandledConstruct("nonlinear constraint"));
  }

  // Sets a logical constraint expression.
  // index: Index of a logical contraint; 0 <= index < num_logical_cons.
  void SetLogicalCon(int index, LogicalExpr expr) {
    MP_UNUSED(index); MP_UNUSED(expr);
    MP_DISPATCH(ReportUnhandledConstruct("logical constraint"));
  }

  // Sets a defined variable expression.
  // index: Index of a defined variable;
  // num_vars <= index < num_vars + num_defined_vars.
  void SetVar(int index, NumericExpr expr, int position) {
    MP_UNUSED(index); MP_UNUSED(expr); MP_UNUSED(position);
    MP_DISPATCH(ReportUnhandledConstruct("nonlinear defined variable"));
  }

  // Sets a complementarity relation.
  void SetComplement(int con_index, int var_index, int flags) {
    MP_UNUSED(con_index); MP_UNUSED(var_index); MP_UNUSED(flags);
    MP_DISPATCH(ReportUnhandledConstruct("complementarity constraint"));
  }

  typedef LinearExprBuilder LinearObjBuilder;

  // Returns a handler for receiving linear terms in an objective.
  LinearObjBuilder GetLinearObjBuilder(int obj_index, int num_linear_terms) {
    MP_UNUSED(obj_index); MP_UNUSED(num_linear_terms);
    MP_DISPATCH(ReportUnhandledConstruct("linear objective"));
    return LinearObjBuilder();
  }

  typedef LinearExprBuilder LinearConBuilder;

  // Returns a handler for receiving linear terms in a constraint.
  LinearConBuilder GetLinearConBuilder(int con_index, int num_linear_terms) {
    MP_UNUSED(con_index); MP_UNUSED(num_linear_terms);
    MP_DISPATCH(ReportUnhandledConstruct("linear constraint"));
    return LinearConBuilder();
  }

  typedef LinearExprBuilder LinearVarBuilder;

  // Returns a handler for receiving linear terms in a defined variable
  // expression.
  LinearVarBuilder GetLinearVarBuilder(int var_index, int num_linear_terms) {
    MP_UNUSED(var_index); MP_UNUSED(num_linear_terms);
    MP_DISPATCH(ReportUnhandledConstruct("linear defined variable"));
    return LinearVarBuilder();
  }

  void SetVarBounds(int index, double lb, double ub) {
    MP_UNUSED(index); MP_UNUSED(lb); MP_UNUSED(ub);
    MP_DISPATCH(ReportUnhandledConstruct("variable bound"));
  }
  void SetConBounds(int index, double lb, double ub) {
    MP_UNUSED(index); MP_UNUSED(lb); MP_UNUSED(ub);
    MP_DISPATCH(ReportUnhandledConstruct("constraint bound"));
  }

  void SetInitialValue(int var_index, double value) {
    MP_UNUSED(var_index); MP_UNUSED(value);
    MP_DISPATCH(ReportUnhandledConstruct("initial value"));
  }
  void SetInitialDualValue(int con_index, double value) {
    MP_UNUSED(con_index); MP_UNUSED(value);
    MP_DISPATCH(ReportUnhandledConstruct("initial dual value"));
  }

  struct ColumnSizeHandler {
    void Add(int size) { MP_UNUSED(size); }
  };

  // Returns a handler that receives column sizes in Jacobian.
  ColumnSizeHandler GetColumnSizeHandler() {
    MP_DISPATCH(ReportUnhandledConstruct("Jacobian column size"));
    return ColumnSizeHandler();
  }

  // Sets a function at the given index.
  void SetFunction(int index, fmt::StringRef name,
                   int num_args, func::Type type) {
    MP_UNUSED(index); MP_UNUSED(name); MP_UNUSED(num_args); MP_UNUSED(type);
    MP_DISPATCH(ReportUnhandledConstruct("function"));
  }

  struct SuffixHandler {
    void SetValue(int index, int value) {
      MP_UNUSED(index); MP_UNUSED(value);
    }
    void SetValue(int index, double value) {
      MP_UNUSED(index); MP_UNUSED(value);
    }
  };

  // Adds a suffix.
  SuffixHandler AddSuffix(int kind, int num_values, fmt::StringRef name) {
    MP_UNUSED(kind); MP_UNUSED(num_values); MP_UNUSED(name);
    MP_DISPATCH(ReportUnhandledConstruct("suffix"));
    return SuffixHandler();
  }

  typedef ArgHandler NumericArgHandler;
  typedef ArgHandler LogicalArgHandler;
  typedef ArgHandler CallArgHandler;

  NumericExpr MakeNumericConstant(double value) {
    MP_UNUSED(value);
    MP_DISPATCH(ReportUnhandledConstruct(
                  "numeric constant in nonlinear expression"));
    return NumericExpr();
  }

  Variable MakeVariable(int var_index) {
    MP_UNUSED(var_index);
    MP_DISPATCH(ReportUnhandledConstruct("nonlinear variable"));
    return Variable();
  }

  NumericExpr MakeUnary(expr::Kind kind, NumericExpr arg) {
    MP_UNUSED(kind); MP_UNUSED(arg);
    MP_DISPATCH(ReportUnhandledConstruct("unary expression"));
    return NumericExpr();
  }

  NumericExpr MakeBinary(expr::Kind kind, NumericExpr lhs, NumericExpr rhs) {
    MP_UNUSED(kind); MP_UNUSED(lhs); MP_UNUSED(rhs);
    MP_DISPATCH(ReportUnhandledConstruct("binary expression"));
    return NumericExpr();
  }

  NumericExpr MakeIf(LogicalExpr condition,
      NumericExpr true_expr, NumericExpr false_expr) {
    MP_UNUSED(condition); MP_UNUSED(true_expr); MP_UNUSED(false_expr);
    MP_DISPATCH(ReportUnhandledConstruct("if expression"));
    return NumericExpr();
  }

  struct PLTermHandler {
    void AddSlope(double slope) { MP_UNUSED(slope); }
    void AddBreakpoint(double breakpoint) { MP_UNUSED(breakpoint); }
  };

  PLTermHandler BeginPLTerm(int num_breakpoints) {
    MP_UNUSED(num_breakpoints);
    MP_DISPATCH(ReportUnhandledConstruct("piecewise-linear term"));
    return PLTermHandler();
  }
  NumericExpr EndPLTerm(PLTermHandler handler, Variable var) {
    MP_UNUSED(handler); MP_UNUSED(var);
    return NumericExpr();
  }

  CallArgHandler BeginCall(int func_index, int num_args) {
    MP_UNUSED(func_index); MP_UNUSED(num_args);
    MP_DISPATCH(ReportUnhandledConstruct("function call"));
    return CallArgHandler();
  }
  NumericExpr EndCall(CallArgHandler handler) {
    MP_UNUSED(handler);
    return NumericExpr();
  }

  NumericArgHandler BeginVarArg(expr::Kind kind, int num_args) {
    MP_UNUSED(kind); MP_UNUSED(num_args);
    MP_DISPATCH(ReportUnhandledConstruct("variable-argument expression"));
    return NumericArgHandler();
  }
  NumericExpr EndVarArg(NumericArgHandler handler) {
    MP_UNUSED(handler);
    return NumericExpr();
  }

  NumericArgHandler BeginSum(int num_args) {
    MP_UNUSED(num_args);
    MP_DISPATCH(ReportUnhandledConstruct("sum"));
    return NumericArgHandler();
  }
  NumericExpr EndSum(NumericArgHandler handler) {
    MP_UNUSED(handler);
    return NumericExpr();
  }

  LogicalArgHandler BeginCount(int num_args) {
    MP_UNUSED(num_args);
    MP_DISPATCH(ReportUnhandledConstruct("count expression"));
    return LogicalArgHandler();
  }
  NumericExpr EndCount(LogicalArgHandler handler) {
    MP_UNUSED(handler);
    return NumericExpr();
  }

  NumericArgHandler BeginNumberOf(int num_args) {
    MP_UNUSED(num_args);
    MP_DISPATCH(ReportUnhandledConstruct("numberof expression"));
    return NumericArgHandler();
  }
  NumericExpr EndNumberOf(NumericArgHandler handler) {
    MP_UNUSED(handler);
    return NumericExpr();
  }

  LogicalExpr MakeLogicalConstant(bool value) {
    MP_UNUSED(value);
    MP_DISPATCH(ReportUnhandledConstruct("logical constant"));
    return LogicalExpr();
  }

  LogicalExpr MakeNot(LogicalExpr arg) {
    MP_UNUSED(arg);
    MP_DISPATCH(ReportUnhandledConstruct("logical not"));
    return LogicalExpr();
  }

  LogicalExpr MakeBinaryLogical(
      expr::Kind kind, LogicalExpr lhs, LogicalExpr rhs) {
    MP_UNUSED(kind); MP_UNUSED(lhs); MP_UNUSED(rhs);
    MP_DISPATCH(ReportUnhandledConstruct("binary logical expression"));
    return LogicalExpr();
  }

  LogicalExpr MakeRelational(
      expr::Kind kind, NumericExpr lhs, NumericExpr rhs) {
    MP_UNUSED(kind); MP_UNUSED(lhs); MP_UNUSED(rhs);
    MP_DISPATCH(ReportUnhandledConstruct("relational expression"));
    return LogicalExpr();
  }

  LogicalExpr MakeLogicalCount(
      expr::Kind kind, NumericExpr lhs, CountExpr rhs) {
    MP_UNUSED(kind); MP_UNUSED(lhs); MP_UNUSED(rhs);
    MP_DISPATCH(ReportUnhandledConstruct("logical count expression"));
    return LogicalExpr();
  }

  LogicalExpr MakeImplication(
      LogicalExpr condition, LogicalExpr true_expr, LogicalExpr false_expr) {
    MP_UNUSED(condition); MP_UNUSED(true_expr); MP_UNUSED(false_expr);
    MP_DISPATCH(ReportUnhandledConstruct("implication expression"));
    return LogicalExpr();
  }

  LogicalArgHandler BeginIteratedLogical(expr::Kind kind, int num_args) {
    MP_UNUSED(kind); MP_UNUSED(num_args);
    MP_DISPATCH(ReportUnhandledConstruct("iterated logical expression"));
    return LogicalArgHandler();
  }
  LogicalExpr EndIteratedLogical(LogicalArgHandler handler) {
    MP_UNUSED(handler);
    return LogicalExpr();
  }

  typedef ArgHandler AllDiffArgHandler;

  AllDiffArgHandler BeginAllDiff(int num_args) {
    MP_UNUSED(num_args);
    MP_DISPATCH(ReportUnhandledConstruct("alldiff expression"));
    return AllDiffArgHandler();
  }
  LogicalExpr EndAllDiff(AllDiffArgHandler handler) {
    MP_UNUSED(handler);
    return LogicalExpr();
  }

  // Constructs a StringLiteral object.
  // value: string value which may not be null-terminated.
  Expr MakeStringLiteral(fmt::StringRef value) {
    MP_UNUSED(value);
    MP_DISPATCH(ReportUnhandledConstruct("string literal"));
    return Expr();
  }
};
}  // namespace mp

#endif  // MP_PROBLEM_BUILDER_H_
