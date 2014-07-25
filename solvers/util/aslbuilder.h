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

#include "solvers/util/expr.h"
#include "solvers/util/noncopyable.h"

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
class ASLBuilder : Noncopyable {
 private:
  ASL *asl_;
  bool own_asl_;
  efunc **r_ops_;
  efunc *standard_opcodes_[N_OPS];
  int nv1_;
  int nz_;
  int nderp_;
  static const double DVALUE[];

  static void CheckOpCode(int opcode, Expr::Kind kind, const char *expr_name) {
    if (Expr::INFO[opcode].kind != kind)
      throw Error("invalid {} expression code {}", expr_name, opcode);
  }

  // Creates a binary or unary expression.
  expr *MakeExpr(int opcode, Expr lhs, Expr rhs = Expr());

  template <typename T>
  T *Allocate(unsigned size = sizeof(T)) {
    assert(size >= sizeof(T));
    return reinterpret_cast<T*>(mem_ASL(asl_, size));
  }

  expr *MakeConstant(double value);

  expr *MakeBinary(int opcode, Expr::Kind kind, Expr lhs, Expr rhs);

  expr *MakeIf(int opcode,
      LogicalExpr condition, Expr true_expr, Expr false_expr);

  // Makes an iterated expression.
  expr *MakeIterated(int opcode, int num_args, const Expr *args);

  // Makes an iterated expression.
  template <Expr::Kind K>
  BasicIteratedExpr<K> MakeIterated(
      int opcode, int num_args, const Expr *args) {
    return Expr::Create< BasicIteratedExpr<K> >(
        MakeIterated(opcode, num_args, args));
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
  void BeginBuild(const char *stub, const NLHeader &h, int flags);

  // Ends building the ASL object.
  void EndBuild();

  // Adds an objective.
  void AddObj(int obj_index, bool maximize, NumericExpr expr);

  Function AddFunction(int index, const char *name, int num_args, int type = 0);

  NumericConstant MakeNumericConstant(double value) {
    return Expr::Create<NumericConstant>(MakeConstant(value));
  }

  Variable MakeVariable(int var_index);

  // The Make* methods construct expression objects. These objects are
  // local to the currently built ASL problem and shouldn't be used with
  // other problems. The expression objects are not accessible via the
  // ASL API until they are added to the problem as a part of objective
  // or constraint expression. For this reason the methods below use a
  // different naming convention from the Add* methods.

  UnaryExpr MakeUnary(int opcode, NumericExpr arg);

  BinaryExpr MakeBinary(int opcode, NumericExpr lhs, NumericExpr rhs) {
    return Expr::Create<BinaryExpr>(MakeBinary(opcode, Expr::BINARY, lhs, rhs));
  }

  IfExpr MakeIf(LogicalExpr condition,
      NumericExpr true_expr, NumericExpr false_expr) {
    return Expr::Create<IfExpr>(
        MakeIf(OPIFnl, condition, true_expr, false_expr));
  }

  PiecewiseLinearExpr MakePiecewiseLinear(int num_breakpoints,
      const double *breakpoints, const double *slopes, Variable var);

  CallExpr MakeCall(Function f, int num_args, const Expr *args);

  VarArgExpr MakeVarArg(int opcode, int num_args, NumericExpr *args);

  SumExpr MakeSum(int num_args, const NumericExpr *args) {
    return MakeIterated<Expr::SUM>(OPSUMLIST, num_args, args);
  }

  CountExpr MakeCount(int num_args, const LogicalExpr *args) {
    return MakeIterated<Expr::COUNT>(OPCOUNT, num_args, args);
  }

  NumberOfExpr MakeNumberOf(int num_args, const NumericExpr *args) {
    assert(num_args >= 1);
    return MakeIterated<Expr::NUMBEROF>(OPNUMBEROF, num_args, args);
  }

  LogicalConstant MakeLogicalConstant(bool value) {
    return Expr::Create<LogicalConstant>(MakeConstant(value));
  }

  NotExpr MakeNot(LogicalExpr arg);

  BinaryLogicalExpr MakeBinaryLogical(
      int opcode, LogicalExpr lhs, LogicalExpr rhs) {
    return Expr::Create<BinaryLogicalExpr>(
        MakeBinary(opcode, Expr::BINARY_LOGICAL, lhs, rhs));
  }

  RelationalExpr MakeRelational(int opcode, NumericExpr lhs, NumericExpr rhs) {
    return Expr::Create<RelationalExpr>(
        MakeBinary(opcode, Expr::RELATIONAL, lhs, rhs));
  }

  LogicalCountExpr MakeLogicalCount(
      int opcode, NumericExpr lhs, CountExpr rhs) {
    return Expr::Create<LogicalCountExpr>(
        MakeBinary(opcode, Expr::LOGICAL_COUNT, lhs, rhs));
  }

  ImplicationExpr MakeImplication(
      LogicalExpr condition, LogicalExpr true_expr, LogicalExpr false_expr) {
    return Expr::Create<ImplicationExpr>(
        MakeIf(OPIMPELSE, condition, true_expr, false_expr));
  }

  IteratedLogicalExpr MakeIteratedLogical(
      int opcode, int num_args, const LogicalExpr *args);

  AllDiffExpr MakeAllDiff(int num_args, const NumericExpr *args) {
    return MakeIterated<Expr::ALLDIFF>(OPALLDIFF, num_args, args);
  }

  StringLiteral MakeStringLiteral(int size, const char *value);
};
}
}

#endif  // SOLVERS_UTIL_ASLBUILDER_H_
