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

#ifndef MP_ASLBUILDER_H_
#define MP_ASLBUILDER_H_

#include "mp/format.h"
#include "mp/problem-base.h"
#include "mp/safeint.h"
#include "expr.h"

struct Static;

namespace mp {

struct NLHeader;
class ExprFactory;

namespace internal {

// An exception representing an ASL error.
class ASLError: public mp::Error {
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

  void SetBounds(double *lbs, double *&ubs, int index, double lb, double ub) {
    if (!ubs)
      ubs = lbs;
    if (lbs != ubs) {
      lbs[index] = lb;
      ubs[index] = ub;
    } else {
      lbs[2 * index] = lb;
      ubs[2 * index + 1] = ub;
    }
  }

  static void CheckOpCode(int opcode, expr::Kind kind, const char *expr_name) {
    if (expr::kind(opcode) != kind)
      throw Error("invalid {} expression code {}", expr_name, opcode);
  }

  template <typename T>
  T *Allocate(SafeInt<int> size = SafeInt<int>(sizeof(T)));

  template <typename T>
  T *ZapAllocate(std::size_t size);

  template <typename T>
  T *AllocateSuffixValues(T *&values, int num_values, int nx, int nx1);

  // Sets objective or constraint expression; adapted from co_read.
  void SetObjOrCon(int index, cde *d, int *cexp1_end, ::expr *e, int **z);

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
  typedef mp::Expr Expr;
  typedef mp::NumericExpr NumericExpr;
  typedef mp::LogicalExpr LogicalExpr;
  typedef mp::Variable Variable;
  typedef mp::CountExpr CountExpr;

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

  void SetVarBounds(int index, double lb, double ub) {
    SetBounds(asl_->i.LUv_, asl_->i.Uvx_, index, lb, ub);
  }
  void SetConBounds(int index, double lb, double ub) {
    SetBounds(asl_->i.LUrhs_, asl_->i.Urhsx_, index, lb, ub);
  }

  void SetComplement(int, int, int) {
    // TODO
  }

  void SetVar(int, NumericExpr, int) {
    // TODO
  }

  // Sets objective type and expression.
  // index: Index of an objective; 0 <= index < num_objs.
  void SetObj(int index, obj::Type type, NumericExpr expr);

  // Sets an algebraic constraint expression.
  // index: Index of an algebraic contraint; 0 <= index < num_algebraic_cons.
  void SetCon(int index, NumericExpr expr);

  // Sets a logical constraint expression.
  // index: Index of a logical contraint; 0 <= index < num_logical_cons.
  void SetLogicalCon(int index, LogicalExpr expr);

  class LinearExprHandler {
   public:
    void AddTerm(int, double) {
      // TODO
    }
  };
  LinearExprHandler GetLinearVarHandler(int, int) {
    // TODO
    return LinearExprHandler();
  }
  LinearExprHandler GetLinearObjHandler(int, int) {
    // TODO
    return LinearExprHandler();
  }
  LinearExprHandler GetLinearConHandler(int, int) {
    // TODO
    return LinearExprHandler();
  }

  struct ColumnSizeHandler {
    void Add(int) {
      // TODO
    }
  };
  ColumnSizeHandler GetColumnSizeHandler() {
    // TODO
    return ColumnSizeHandler();
  }

  void SetInitialValue(int, double) {}
  void SetInitialDualValue(int, double) {}

  Function AddFunction(const char *name, ufunc f, int num_args,
                       func::Type type = func::NUMERIC, void *info = 0);

  // Sets a function at the given index.
  // If the function with the specified name doesn't exist and the flag
  // ASL_allow_missing_funcs is not set, SetFunction throws ASLError.
  Function SetFunction(int index, fmt::StringRef name, int num_args,
                       func::Type type = func::NUMERIC);

  class SuffixHandler {
   private:
    int *int_values_;
    double *dbl_values_;

   public:
    explicit SuffixHandler(int *values = 0)
      : int_values_(values), dbl_values_(0) {}
    explicit SuffixHandler(double *values)
      : int_values_(0), dbl_values_(values) {}

    void SetValue(int index, int value) {
      if (int_values_)
        int_values_[index] = value;
      else if (dbl_values_)
        dbl_values_[index] = value;
    }

    void SetValue(int index, double value) {
      if (int_values_)
        int_values_[index] = static_cast<int>(value + 0.5);
      else if (dbl_values_)
        dbl_values_[index] = value;
    }
  };

  SuffixHandler AddSuffix(int kind, int num_values, fmt::StringRef name);

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
}  // namespace mp

#endif  // MP_ASLBUILDER_H_
