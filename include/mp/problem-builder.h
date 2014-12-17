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

#include <cassert>
#include <cstring>
#include <memory>
#include <set>
#include <vector>

#include "mp/common.h"
#include "mp/error.h"
#include "mp/nl.h"
#include "mp/suffix.h"

namespace mp {

// A minimal implementation of the ProblemBuilder concept.
template <typename Impl, typename ExprType>
class ProblemBuilder : public SuffixManager {
 private:
  struct ExprBuilder {
    void AddArg(ExprType arg) { MP_UNUSED(arg); }
  };

  template <typename T>
  struct SuffixHandler {
    void SetValue(int index, T value) { MP_UNUSED2(index, value); }
  };

 public:
  typedef ExprType Expr;
  typedef Expr NumericExpr;
  typedef Expr LogicalExpr;
  typedef Expr CountExpr;
  typedef Expr Reference;

  static void ReportUnhandledConstruct(fmt::StringRef name) {
    throw MakeUnsupportedError(name);
  }

  void SetInfo(const ProblemInfo &) {}
  void EndBuild() {}

  // Adds a variable.
  void AddVar(double lb, double ub, var::Type type) {
    MP_UNUSED3(lb, ub, type);
    MP_DISPATCH(ReportUnhandledConstruct("variable"));
  }

  struct LinearExprBuilder {
    void AddTerm(int var_index, double coef) { MP_UNUSED2(var_index, coef); }
  };

  typedef LinearExprBuilder LinearObjBuilder;

  // Adds an objective.
  // Returns a builder for the linear part of the objective expression.
  LinearObjBuilder AddObj(
      obj::Type type, NumericExpr expr, int num_linear_terms) {
    MP_UNUSED3(type, expr, num_linear_terms);
    MP_DISPATCH(ReportUnhandledConstruct("objective"));
    return LinearObjBuilder();
  }

  typedef LinearExprBuilder LinearConBuilder;

  // Adds an algebraic constraint.
  // Returns a builder for the linear part of the constraint expression.
  LinearConBuilder AddCon(double lb, double ub, NumericExpr expr,
                          int num_linear_terms) {
    MP_UNUSED3(expr, lb, ub); MP_UNUSED(num_linear_terms);
    MP_DISPATCH(ReportUnhandledConstruct("algebraic constraint"));
    return LinearConBuilder();
  }

  // Adds a logical constraint.
  void AddCon(LogicalExpr expr) {
    MP_UNUSED(expr);
    MP_DISPATCH(ReportUnhandledConstruct("logical constraint"));
  }

  // Adds a common expression (defined variable).
  // Returns a builder for the linear part of the common expression.
  LinearExprBuilder BeginCommonExpr(int num_linear_terms) {
    MP_UNUSED(num_linear_terms);
    MP_DISPATCH(ReportUnhandledConstruct("common expression"));
    return LinearExprBuilder();
  }

  void EndCommonExpr(LinearExprBuilder builder,
                     NumericExpr expr, int position) {
    MP_UNUSED3(builder, expr, position);
  }

  // Sets a complementarity relation.
  void SetComplement(int con_index, int var_index, int flags) {
    MP_UNUSED3(con_index, var_index, flags);
    MP_DISPATCH(ReportUnhandledConstruct("complementarity constraint"));
  }

  void SetInitialValue(int var_index, double value) {
    MP_UNUSED2(var_index, value);
    // Initial values are ignored by default.
  }
  void SetInitialDualValue(int con_index, double value) {
    MP_UNUSED2(con_index, value);
    // Initial dual values are ignored by default.
  }

  struct ColumnSizeHandler {
    void Add(int size) { MP_UNUSED(size); }
  };

  // Returns a handler that receives column sizes in Jacobian.
  ColumnSizeHandler GetColumnSizeHandler() {
    MP_DISPATCH(ReportUnhandledConstruct("Jacobian column size"));
    return ColumnSizeHandler();
  }

  struct Function {};

  // Adds a function.
  Function AddFunction(fmt::StringRef name, int num_args, func::Type type) {
    MP_UNUSED3(name, num_args, type);
    MP_DISPATCH(ReportUnhandledConstruct("function"));
    return Function();
  }

  typedef SuffixHandler<int> IntSuffixHandler;

  // Adds a suffix.
  IntSuffixHandler AddIntSuffix(fmt::StringRef name, int kind, int num_values) {
    MP_UNUSED3(kind, num_values, name);
    MP_DISPATCH(ReportUnhandledConstruct("integer suffix"));
    return IntSuffixHandler();
  }

  typedef SuffixHandler<double> DblSuffixHandler;

  // Adds a suffix.
  DblSuffixHandler AddDblSuffix(fmt::StringRef name, int kind, int num_values) {
    MP_UNUSED3(kind, num_values, name);
    MP_DISPATCH(ReportUnhandledConstruct("double suffix"));
    return DblSuffixHandler();
  }

  typedef ExprBuilder NumericExprBuilder;
  typedef ExprBuilder VarArgExprBuilder;
  typedef ExprBuilder CallExprBuilder;
  typedef ExprBuilder NumberOfExprBuilder;
  typedef ExprBuilder CountExprBuilder;
  typedef ExprBuilder IteratedLogicalExprBuilder;

  NumericExpr MakeNumericConstant(double value) {
    MP_UNUSED(value);
    MP_DISPATCH(ReportUnhandledConstruct(
                  "numeric constant in nonlinear expression"));
    return NumericExpr();
  }

  Reference MakeVariable(int var_index) {
    MP_UNUSED(var_index);
    MP_DISPATCH(ReportUnhandledConstruct("variable in nonlinear expression"));
    return Reference();
  }

  NumericExpr MakeUnary(expr::Kind kind, NumericExpr arg) {
    MP_UNUSED2(kind, arg);
    MP_DISPATCH(ReportUnhandledConstruct(str(kind)));
    return NumericExpr();
  }

  NumericExpr MakeBinary(expr::Kind kind, NumericExpr lhs, NumericExpr rhs) {
    MP_UNUSED3(kind, lhs, rhs);
    MP_DISPATCH(ReportUnhandledConstruct(str(kind)));
    return NumericExpr();
  }

  NumericExpr MakeIf(LogicalExpr condition,
      NumericExpr true_expr, NumericExpr false_expr) {
    MP_UNUSED3(condition, true_expr, false_expr);
    MP_DISPATCH(ReportUnhandledConstruct("if expression"));
    return NumericExpr();
  }

  struct PLTermBuilder {
    void AddSlope(double slope) { MP_UNUSED(slope); }
    void AddBreakpoint(double breakpoint) { MP_UNUSED(breakpoint); }
  };

  PLTermBuilder BeginPLTerm(int num_breakpoints) {
    MP_UNUSED(num_breakpoints);
    MP_DISPATCH(ReportUnhandledConstruct("piecewise-linear term"));
    return PLTermBuilder();
  }
  NumericExpr EndPLTerm(PLTermBuilder builder, Reference arg) {
    MP_UNUSED2(builder, arg);
    MP_DISPATCH(ReportUnhandledConstruct("piecewise-linear term"));
    return NumericExpr();
  }

  CallExprBuilder BeginCall(Function func, int num_args) {
    MP_UNUSED2(func, num_args);
    MP_DISPATCH(ReportUnhandledConstruct("function call"));
    return CallExprBuilder();
  }
  NumericExpr EndCall(CallExprBuilder builder) {
    MP_UNUSED(builder);
    MP_DISPATCH(ReportUnhandledConstruct("function call"));
    return NumericExpr();
  }

  VarArgExprBuilder BeginVarArg(expr::Kind kind, int num_args) {
    MP_UNUSED2(kind, num_args);
    MP_DISPATCH(ReportUnhandledConstruct(str(kind)));
    return VarArgExprBuilder();
  }
  NumericExpr EndVarArg(VarArgExprBuilder builder) {
    MP_UNUSED(builder);
    MP_DISPATCH(ReportUnhandledConstruct("vararg expression"));
    return NumericExpr();
  }

  NumericExprBuilder BeginSum(int num_args) {
    MP_UNUSED(num_args);
    MP_DISPATCH(ReportUnhandledConstruct("sum"));
    return NumericExprBuilder();
  }
  NumericExpr EndSum(NumericExprBuilder builder) {
    MP_UNUSED(builder);
    MP_DISPATCH(ReportUnhandledConstruct("sum"));
    return NumericExpr();
  }

  NumberOfExprBuilder BeginNumberOf(int num_args, NumericExpr arg0) {
    MP_UNUSED2(num_args, arg0);
    MP_DISPATCH(ReportUnhandledConstruct("numberof expression"));
    return NumberOfExprBuilder();
  }
  NumericExpr EndNumberOf(NumberOfExprBuilder builder) {
    MP_UNUSED(builder);
    MP_DISPATCH(ReportUnhandledConstruct("numberof expression"));
    return NumericExpr();
  }

  typedef ExprBuilder SymbolicNumberOfExprBuilder;

  SymbolicNumberOfExprBuilder BeginSymbolicNumberOf(int num_args, Expr arg0) {
    MP_UNUSED2(num_args, arg0);
    MP_DISPATCH(ReportUnhandledConstruct("symbolic numberof expression"));
    return SymbolicNumberOfExprBuilder();
  }
  NumericExpr EndSymbolicNumberOf(SymbolicNumberOfExprBuilder builder) {
    MP_UNUSED(builder);
    MP_DISPATCH(ReportUnhandledConstruct("symbolic numberof expression"));
    return NumericExpr();
  }

  CountExprBuilder BeginCount(int num_args) {
    MP_UNUSED(num_args);
    MP_DISPATCH(ReportUnhandledConstruct("count expression"));
    return CountExprBuilder();
  }
  NumericExpr EndCount(CountExprBuilder builder) {
    MP_UNUSED(builder);
    MP_DISPATCH(ReportUnhandledConstruct("count expression"));
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
    MP_UNUSED3(kind, lhs, rhs);
    MP_DISPATCH(ReportUnhandledConstruct(str(kind)));
    return LogicalExpr();
  }

  LogicalExpr MakeRelational(
      expr::Kind kind, NumericExpr lhs, NumericExpr rhs) {
    MP_UNUSED3(kind, lhs, rhs);
    MP_DISPATCH(ReportUnhandledConstruct(str(kind)));
    return LogicalExpr();
  }

  LogicalExpr MakeLogicalCount(
      expr::Kind kind, NumericExpr lhs, CountExpr rhs) {
    MP_UNUSED3(kind, lhs, rhs);
    MP_DISPATCH(ReportUnhandledConstruct(str(kind)));
    return LogicalExpr();
  }

  LogicalExpr MakeImplication(
      LogicalExpr condition, LogicalExpr true_expr, LogicalExpr false_expr) {
    MP_UNUSED3(condition, true_expr, false_expr);
    MP_DISPATCH(ReportUnhandledConstruct("implication expression"));
    return LogicalExpr();
  }

  IteratedLogicalExprBuilder BeginIteratedLogical(
      expr::Kind kind, int num_args) {
    MP_UNUSED2(kind, num_args);
    MP_DISPATCH(ReportUnhandledConstruct(str(kind)));
    return IteratedLogicalExprBuilder();
  }
  LogicalExpr EndIteratedLogical(IteratedLogicalExprBuilder builder) {
    MP_UNUSED(builder);
    MP_DISPATCH(ReportUnhandledConstruct("iterated logical expression"));
    return LogicalExpr();
  }

  typedef ExprBuilder PairwiseExprBuilder;

  PairwiseExprBuilder BeginPairwise(expr::Kind kind, int num_args) {
    MP_UNUSED2(kind, num_args);
    MP_DISPATCH(ReportUnhandledConstruct("alldiff expression"));
    return PairwiseExprBuilder();
  }
  LogicalExpr EndPairwise(PairwiseExprBuilder builder) {
    MP_UNUSED(builder);
    MP_DISPATCH(ReportUnhandledConstruct("alldiff expression"));
    return LogicalExpr();
  }

  // Constructs a StringLiteral object.
  // value: string value which may not be null-terminated.
  Expr MakeStringLiteral(fmt::StringRef value) {
    MP_UNUSED(value);
    MP_DISPATCH(ReportUnhandledConstruct("string literal"));
    return Expr();
  }

  Expr MakeSymbolicIf(LogicalExpr condition, Expr true_expr, Expr false_expr) {
    MP_UNUSED3(condition, true_expr, false_expr);
    MP_DISPATCH(ReportUnhandledConstruct("symbolic if expression"));
    return Expr();
  }
};

// Adapts ProblemBuilder for use as an .nl handler.
// It doesn't call ProblemBuilder::EndBuild to allow for modification
// of a problem after it has been read.
template <typename ProblemBuilder>
class ProblemBuilderToNLAdapter {
 public:
  typedef typename ProblemBuilder::Function Function;
  typedef typename ProblemBuilder::Expr Expr;
  typedef typename ProblemBuilder::NumericExpr NumericExpr;
  typedef typename ProblemBuilder::LogicalExpr LogicalExpr;
  typedef typename ProblemBuilder::CountExpr CountExpr;
  typedef typename ProblemBuilder::Reference Reference;

 private:
  ProblemBuilder &builder_;
  int num_continuous_vars_;

  struct ObjInfo {
    obj::Type type;
    NumericExpr expr;
    ObjInfo() : type(obj::MIN), expr() {}
  };
  std::vector<ObjInfo> objs_;

  int obj_index_;

  // Algebraic constraints
  struct ConInfo {
    NumericExpr expr;
    double lb, ub;
    ConInfo() : expr(), lb(0), ub(0) {}
  };
  std::vector<ConInfo> cons_;

  std::vector<Function> funcs_;

 protected:
  void set_obj_index(int index) { obj_index_ = index; }

 public:
  // Possible values for the objective index.
  enum {
    // Skip all objectives.
    SKIP_ALL_OBJS = -1,

    // Pass all objectives to the builder.
    NEED_ALL_OBJS = -2
  };

  // TODO: check for permuted indices in segments

  // Index of the objective to pass to the builder, SKIP_ALL_OBJS to
  // skip all objectives, NEED_ALL_OBJS to pass all objectives.
  int obj_index() const { return obj_index_; }

  explicit ProblemBuilderToNLAdapter(ProblemBuilder &builder, int obj_index = 0)
    : builder_(builder), num_continuous_vars_(0), obj_index_(obj_index) {}

  ProblemBuilder &builder() { return builder_; }

  // Receives notification of an .nl header.
  void OnHeader(const NLHeader &h) {
    num_continuous_vars_ = h.num_continuous_vars();
    objs_.resize(h.num_objs);
    cons_.resize(h.num_algebraic_cons);
    funcs_.resize(h.num_funcs);

    // Update the number of objectives if necessary.
    int num_objs = 0;
    if (obj_index_ >= 0)
      num_objs = std::min(h.num_objs, 1);
    else if (obj_index_ == NEED_ALL_OBJS)
      num_objs = h.num_objs;
    if (num_objs == h.num_objs) {
      builder_.SetInfo(h);
    } else {
      ProblemInfo info(h);
      info.num_objs = num_objs;
      builder_.SetInfo(info);
    }
  }

  // Returns true if objective should be handled.
  bool NeedObj(int obj_index) const {
    if (obj_index == obj_index_)
      return true;
    return obj_index_ == NEED_ALL_OBJS;
  }

  // Receives notification of an objective type and the nonlinear part of
  // an objective expression.
  void OnObj(int index, obj::Type type, NumericExpr expr) {
    assert(0 <= index && static_cast<unsigned>(index) < objs_.size());
    if (!NeedObj(index))
      return;  // Ignore inactive objective.
    ObjInfo &obj = objs_[index];
    obj.type = type;
    obj.expr = expr;
  }

  // Receives notification of the nonlinear part of an algebraic constraint
  // expression.
  void OnAlgebraicCon(int index, NumericExpr expr) {
    assert(0 <= index && static_cast<unsigned>(index) < cons_.size());
    cons_[index].expr = expr;
  }

  // Receives notification of a logical constraint.
  void OnLogicalCon(int, LogicalExpr expr) {
    builder_.AddCon(expr);
  }

  // Receives notification of a complementarity relation.
  void OnComplement(int con_index, int var_index, int flags) {
    builder_.SetComplement(con_index, var_index, flags);
  }

  typedef typename ProblemBuilder::LinearObjBuilder LinearObjHandler;

  // Receives notification of the linear part of an objective expression.
  LinearObjHandler OnLinearObjExpr(int obj_index, int num_linear_terms) {
    assert(0 <= obj_index && static_cast<unsigned>(obj_index) < objs_.size());
    const ObjInfo &obj_info = objs_[obj_index];
    return builder_.AddObj(obj_info.type, obj_info.expr, num_linear_terms);
  }

  typedef typename ProblemBuilder::LinearConBuilder LinearConHandler;

  // Receives notification of the linear part of a constraint expression.
  LinearConHandler OnLinearConExpr(int con_index, int num_linear_terms) {
    const ConInfo &con = cons_[con_index];
    return builder_.AddCon(con.lb, con.ub, con.expr, num_linear_terms);
  }

  typedef typename ProblemBuilder::LinearExprBuilder LinearExprHandler;

  // Receives notification of a commmon expression (defined variable).
  LinearExprHandler BeginCommonExpr(int index, int num_linear_terms) {
    MP_UNUSED(index);
    return builder_.BeginCommonExpr(num_linear_terms);
  }
  void EndCommonExpr(LinearExprHandler handler,
                     NumericExpr expr, int position) {
    builder_.EndCommonExpr(handler, expr, position);
  }

  // Receives notification of variable bounds.
  void OnVarBounds(int index, double lb, double ub) {
    var::Type type =
        index < num_continuous_vars_ ? var::CONTINUOUS : var::INTEGER;
    builder_.AddVar(lb, ub, type);
  }

  // Receives notification of constraint bounds (ranges).
  void OnConBounds(int index, double lb, double ub) {
    ConInfo &con = cons_[index];
    con.lb = lb;
    con.ub = ub;
  }

  // Receives notification of the initial value for a variable.
  void OnInitialValue(int var_index, double value) {
    builder_.SetInitialValue(var_index, value);
  }

  // Receives notification of the initial value for a dual variable.
  void OnInitialDualValue(int con_index, double value) {
    builder_.SetInitialDualValue(con_index, value);
  }

  typedef typename ProblemBuilder::ColumnSizeHandler ColumnSizeHandler;

  // Receives notification of Jacobian column sizes.
  ColumnSizeHandler OnColumnSizes() {
    return builder_.GetColumnSizeHandler();
  }

  // Receives notification of a function.
  void OnFunction(int index, fmt::StringRef name,
                  int num_args, func::Type type) {
    funcs_[index] = builder_.AddFunction(name, num_args, type);
  }

  typedef typename ProblemBuilder::IntSuffixHandler IntSuffixHandler;

  // Receives notification of an integer suffix.
  IntSuffixHandler OnIntSuffix(fmt::StringRef name, int kind, int num_values) {
    return builder_.AddIntSuffix(name, kind, num_values);
  }

  typedef typename ProblemBuilder::DblSuffixHandler DblSuffixHandler;

  // Receives notification of a double suffix.
  DblSuffixHandler OnDblSuffix(fmt::StringRef name, int kind, int num_values) {
    return builder_.AddDblSuffix(name, kind, num_values);
  }

  typedef typename ProblemBuilder::NumericExprBuilder NumericArgHandler;

  // Receives notification of a numeric constant in a nonlinear expression.
  NumericExpr OnNumericConstant(double value) {
    return builder_.MakeNumericConstant(value);
  }

  // Receives notification of a variable in a nonlinear expression.
  Reference OnVariableRef(int var_index) {
    return builder_.MakeVariable(var_index);
  }

  // Receives notification of a common expression (defined variable) reference.
  Reference OnCommonExprRef(int expr_index) {
    return builder_.MakeCommonExpr(expr_index);
  }

  // Receives notification of a unary expression.
  NumericExpr OnUnary(expr::Kind kind, NumericExpr arg) {
    return builder_.MakeUnary(kind, arg);
  }

  // Receives notification of a binary expression.
  NumericExpr OnBinary(expr::Kind kind, NumericExpr lhs, NumericExpr rhs) {
    return builder_.MakeBinary(kind, lhs, rhs);
  }

  // Receives notification of an if expression.
  NumericExpr OnIf(LogicalExpr condition,
      NumericExpr true_expr, NumericExpr false_expr) {
    return builder_.MakeIf(condition, true_expr, false_expr);
  }

  typedef typename ProblemBuilder::PLTermBuilder PLTermHandler;

  // Receives notification of the beginning of a piecewise-linear term.
  PLTermHandler BeginPLTerm(int num_breakpoints) {
    return builder_.BeginPLTerm(num_breakpoints);
  }
  // Receives notification of the end of a piecewise-linear term.
  NumericExpr EndPLTerm(PLTermHandler handler, Reference arg) {
    return builder_.EndPLTerm(handler, arg);
  }

  typedef typename ProblemBuilder::CallExprBuilder CallArgHandler;

  // Receives notification of the beginning of a call expression.
  CallArgHandler BeginCall(int func_index, int num_args) {
    return builder_.BeginCall(funcs_[func_index], num_args);
  }
  // Receives notification of the end of a call expression.
  NumericExpr EndCall(CallArgHandler handler) {
    return builder_.EndCall(handler);
  }

  typedef typename ProblemBuilder::VarArgExprBuilder VarArgHandler;

  // Receives notification of the beginning of a vararg expression (min or max).
  VarArgHandler BeginVarArg(expr::Kind kind, int num_args) {
    return builder_.BeginVarArg(kind, num_args);
  }
  // Receives notification of the end of a vararg expression (min or max).
  NumericExpr EndVarArg(VarArgHandler handler) {
    return builder_.EndVarArg(handler);
  }

  // Receives notification of the beginning of a sum expression.
  NumericArgHandler BeginSum(int num_args) {
    return builder_.BeginSum(num_args);
  }
  // Receives notification of the end of a sum expression.
  NumericExpr EndSum(NumericArgHandler handler) {
    return builder_.EndSum(handler);
  }

  typedef typename ProblemBuilder::NumberOfExprBuilder NumberOfArgHandler;

  // Receives notification of the beginning of a numberof expression.
  NumberOfArgHandler BeginNumberOf(int num_args, NumericExpr arg0) {
    return builder_.BeginNumberOf(num_args, arg0);
  }
  // Receives notification of the end of a numberof expression.
  NumericExpr EndNumberOf(NumberOfArgHandler handler) {
    return builder_.EndNumberOf(handler);
  }

  typedef typename ProblemBuilder::SymbolicNumberOfExprBuilder
          SymbolicArgHandler;

  // Receives notification of the beginning of a symbolic numberof expression.
  SymbolicArgHandler BeginSymbolicNumberOf(int num_args, Expr arg0) {
    return builder_.BeginSymbolicNumberOf(num_args, arg0);
  }
  // Receives notification of the end of a symbolic numberof expression.
  NumericExpr EndSymbolicNumberOf(SymbolicArgHandler handler) {
    return builder_.EndSymbolicNumberOf(handler);
  }

  typedef typename ProblemBuilder::CountExprBuilder CountArgHandler;

  // Receives notification of the beginning of a count expression.
  CountArgHandler BeginCount(int num_args) {
    return builder_.BeginCount(num_args);
  }
  // Receives notification of the end of a count expression.
  CountExpr EndCount(CountArgHandler handler) {
    return builder_.EndCount(handler);
  }

  // Receives notification of a logical constant.
  LogicalExpr OnLogicalConstant(bool value) {
    return builder_.MakeLogicalConstant(value);
  }

  // Receives notification of a logical not expression.
  LogicalExpr OnNot(LogicalExpr arg) {
    return builder_.MakeNot(arg);
  }

  // Receives notification of a binary logical expression.
  LogicalExpr OnBinaryLogical(
      expr::Kind kind, LogicalExpr lhs, LogicalExpr rhs) {
    return builder_.MakeBinaryLogical(kind, lhs, rhs);
  }

  // Receives notification of a relational expression.
  LogicalExpr OnRelational(
      expr::Kind kind, NumericExpr lhs, NumericExpr rhs) {
    return builder_.MakeRelational(kind, lhs, rhs);
  }

  // Receives notification of a logical count expression.
  LogicalExpr OnLogicalCount(
      expr::Kind kind, NumericExpr lhs, CountExpr rhs) {
    return builder_.MakeLogicalCount(kind, lhs, rhs);
  }

  // Receives notification of an implication expression.
  LogicalExpr OnImplication(
      LogicalExpr condition, LogicalExpr true_expr, LogicalExpr false_expr) {
    return builder_.MakeImplication(condition, true_expr, false_expr);
  }

  typedef typename ProblemBuilder::IteratedLogicalExprBuilder LogicalArgHandler;

  // Receives notification of the beginning of an iterated logical expression.
  LogicalArgHandler BeginIteratedLogical(expr::Kind kind, int num_args) {
    return builder_.BeginIteratedLogical(kind, num_args);
  }
  // Receives notification of the end of an iterated logical expression.
  LogicalExpr EndIteratedLogical(LogicalArgHandler handler) {
    return builder_.EndIteratedLogical(handler);
  }

  typedef typename ProblemBuilder::PairwiseExprBuilder PairwiseArgHandler;

  // Receives notification of the beginning of a pairwise expression.
  PairwiseArgHandler BeginPairwise(expr::Kind kind, int num_args) {
    return builder_.BeginPairwise(kind, num_args);
  }
  // Receives notification of the end of a pairwise expression.
  LogicalExpr EndPairwise(PairwiseArgHandler handler) {
    return builder_.EndPairwise(handler);
  }

  // Receives notification of a string literal.
  Expr OnStringLiteral(fmt::StringRef value) {
    return builder_.MakeStringLiteral(value);
  }

  // Receives notification of a symbolic if expression.
  Expr OnSymbolicIf(LogicalExpr condition, Expr true_expr, Expr false_expr) {
    return builder_.MakeSymbolicIf(condition, true_expr, false_expr);
  }
};
}  // namespace mp

#endif  // MP_PROBLEM_BUILDER_H_
