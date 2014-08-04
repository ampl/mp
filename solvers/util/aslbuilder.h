/*
 An ASL problem builder.

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

#ifndef SOLVERS_UTIL_ASLBUILDER_H_
#define SOLVERS_UTIL_ASLBUILDER_H_

#include "solvers/util/format.h"
#include "solvers/util/problem-base.h"
#include "solvers/util/safeint.h"
#include "solvers/util/expr.h"

struct Static;

namespace ampl {

struct NLHeader;
class ExprFactory;

namespace internal {

// An exception representing an ASL error.
class ASLError: public ampl::Error {
private:
  int error_code_;

public:
  ASLError(int error_code, fmt::StringRef message) :
      Error(message), error_code_(error_code) {}

  int error_code() const { return error_code_; }
};

// Use standard opcodes when building expressions.
enum { ASL_STANDARD_OPCODES = 0x1000000 };

// Provides methods for building an ASL problem object.
class ASLBuilder {
 private:
  ASL *asl_;
  bool own_asl_;
  efunc **r_ops_;
  efunc *standard_opcodes_[N_OPS];
  int flags_;  // Flags passed to BeginBuild.
  int nz_;
  int nderp_;
  static const double DVALUE[];

  // "Static" data for the functions in fg_read.
  Static *static_;

  FMT_DISALLOW_COPY_AND_ASSIGN(ASLBuilder);

  static void CheckOpCode(int opcode, expr::Kind kind, const char *expr_name) {
    if (expr::kind(opcode) != kind)
      throw Error("invalid {} expression code {}", expr_name, opcode);
  }

  template <typename T>
  T *Allocate(safeint::SafeInt<int> size = safeint::SafeInt<int>(sizeof(T)));

  // Sets objective or constraint expression; adapted from co_read.
  void SetObjOrCon(
      int index, cde *d, int *cexp1_end, NumericExpr expr, int **z);

  ::expr *MakeConstant(double value);
  ::expr *DoMakeUnary(int opcode, Expr lhs);
  ::expr *MakeBinary(int opcode, expr::Kind kind, Expr lhs, Expr rhs);
  ::expr *MakeIf(int opcode,
      LogicalExpr condition, Expr true_expr, Expr false_expr);
  ::expr *MakeIterated(int opcode, ArrayRef<Expr> args);

  template <expr::Kind K>
  BasicIteratedExpr<K> MakeIterated(int opcode, ArrayRef<Expr> args) {
    return Expr::Create< BasicIteratedExpr<K> >(MakeIterated(opcode, args));
  }

 public:
  explicit ASLBuilder(ASL *asl = 0);
  ~ASLBuilder();

  // Initializes the ASL object in a similar way to jac0dim, but
  // doesn't read the .nl file as it is the responsibility of NLReader.
  // Instead it uses the information provided in NLHeader.
  void InitASL(const char *stub, const NLHeader &h);

  // Begins building the ASL object.
  // flags: reader flags, see ASL_reader_flag_bits.
  // Throws ASLError on error.
  void BeginBuild(const char *stub, const NLHeader &h,
                  int flags = ASL_STANDARD_OPCODES);

  // Ends building the ASL object.
  void EndBuild();

  // Sets objective type and expression.
  void SetObj(int index, obj::Type type, NumericExpr expr);

  // Sets constraint expression.
  void SetCon(int index, NumericExpr expr);

  Function AddFunction(const char *name, ufunc f, int num_args,
                       func::Type type = func::NUMERIC, void *info = 0);

  // Sets a function at the given index.
  // If the function with the specified name doesn't exist and the flag
  // ASL_allow_missing_funcs is not set, SetFunction throws ASLError.
  Function SetFunction(int index, const char *name, int num_args,
                       func::Type type = func::NUMERIC);

  // The Make* methods construct expression objects. These objects are
  // local to the currently built ASL problem and shouldn't be used with
  // other problems. The expression objects are not accessible via the
  // problem API until they are added as a part of objective or constraint
  // expression. For this reason the methods are called Make* rather than Add*.

  NumericConstant MakeNumericConstant(double value) {
    return Expr::Create<NumericConstant>(MakeConstant(value));
  }

  Variable MakeVariable(int var_index);

  UnaryExpr MakeUnary(int opcode, NumericExpr arg);

  BinaryExpr MakeBinary(int opcode, NumericExpr lhs, NumericExpr rhs) {
    return Expr::Create<BinaryExpr>(MakeBinary(opcode, expr::BINARY, lhs, rhs));
  }

  IfExpr MakeIf(LogicalExpr condition,
      NumericExpr true_expr, NumericExpr false_expr) {
    return Expr::Create<IfExpr>(
        MakeIf(OPIFnl, condition, true_expr, false_expr));
  }

  PiecewiseLinearExpr MakePiecewiseLinear(int num_breakpoints,
      const double *breakpoints, const double *slopes, Variable var);

  CallExpr MakeCall(Function f, ArrayRef<Expr> args);

  CallExpr MakeCall(int func_index, ArrayRef<Expr> args) {
    Function f(asl_->i.funcs_[func_index]);
    if (!f)
      throw Error("function {} undefined", func_index);
    return MakeCall(f, args);
  }

  VarArgExpr MakeVarArg(int opcode, ArrayRef<NumericExpr> args);

  SumExpr MakeSum(ArrayRef<NumericExpr> args) {
    return MakeIterated<expr::SUM>(OPSUMLIST, args);
  }

  CountExpr MakeCount(ArrayRef<LogicalExpr> args) {
    return MakeIterated<expr::COUNT>(OPCOUNT, args);
  }

  NumberOfExpr MakeNumberOf(ArrayRef<NumericExpr> args) {
    assert(args.size() >= 1);
    return MakeIterated<expr::NUMBEROF>(OPNUMBEROF, args);
  }

  LogicalConstant MakeLogicalConstant(bool value) {
    return Expr::Create<LogicalConstant>(MakeConstant(value));
  }

  NotExpr MakeNot(LogicalExpr arg) {
    return Expr::Create<NotExpr>(DoMakeUnary(OPNOT, arg));
  }

  BinaryLogicalExpr MakeBinaryLogical(
      int opcode, LogicalExpr lhs, LogicalExpr rhs) {
    return Expr::Create<BinaryLogicalExpr>(
        MakeBinary(opcode, expr::BINARY_LOGICAL, lhs, rhs));
  }

  RelationalExpr MakeRelational(int opcode, NumericExpr lhs, NumericExpr rhs) {
    return Expr::Create<RelationalExpr>(
        MakeBinary(opcode, expr::RELATIONAL, lhs, rhs));
  }

  LogicalCountExpr MakeLogicalCount(
      int opcode, NumericExpr lhs, CountExpr rhs) {
    return Expr::Create<LogicalCountExpr>(
        MakeBinary(opcode, expr::LOGICAL_COUNT, lhs, rhs));
  }

  ImplicationExpr MakeImplication(
      LogicalExpr condition, LogicalExpr true_expr, LogicalExpr false_expr) {
    return Expr::Create<ImplicationExpr>(
        MakeIf(OPIMPELSE, condition, true_expr, false_expr));
  }

  IteratedLogicalExpr MakeIteratedLogical(
      int opcode, ArrayRef<LogicalExpr> args);

  AllDiffExpr MakeAllDiff(ArrayRef<NumericExpr> args) {
    return MakeIterated<expr::ALLDIFF>(OPALLDIFF, args);
  }

  // Constructs a StringLiteral object.
  // value: string value which may not be null-terminated.
  StringLiteral MakeStringLiteral(fmt::StringRef value);
};
}
}  // namespace ampl

#endif  // SOLVERS_UTIL_ASLBUILDER_H_
