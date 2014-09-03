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
  efunc *standard_opcodes_[expr::MAX_OPCODE + 1];
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

  template <typename ExprT>
  static void CheckKind(expr::Kind kind, const char *expr_name) {
    if (!internal::Is<ExprT>(kind))
      throw Error("invalid {} expression kind {}", expr_name, kind);
  }

  template <typename T>
  T *Allocate(SafeInt<int> size = SafeInt<int>(sizeof(T))) {
    if (size.value() > std::numeric_limits<int>::max())
      throw std::bad_alloc();
    return reinterpret_cast<T*>(mem_ASL(asl_, size.value()));
  }

  template <typename T>
  T *ZapAllocate(std::size_t size);

  template <typename T>
  T *AllocateSuffixValues(T *&values, int num_values, int nx, int nx1);

  // Sets objective or constraint expression; adapted from co_read.
  void SetObjOrCon(int index, cde *d, int *cexp1_end, ::expr *e, int **z);

  ::expr *MakeConstant(double value);

  ::expr *DoMakeUnary(expr::Kind kind, Expr arg);

  template <typename ExprT>
  ExprT MakeUnary(expr::Kind kind, Expr arg) {
    return Expr::Create<ExprT>(DoMakeUnary(kind, arg));
  }

  ::expr *DoMakeBinary(expr::Kind kind, Expr lhs, Expr rhs);

  template <typename ExprT>
  ExprT MakeBinary(expr::Kind kind, Expr lhs, Expr rhs, const char *name) {
    CheckKind<ExprT>(kind, name);
    return Expr::Create<ExprT>(DoMakeBinary(kind, lhs, rhs));
  }

  ::expr *MakeIf(expr::Kind kind,
      LogicalExpr condition, Expr true_expr, Expr false_expr);
  ::expr *MakeIterated(expr::Kind kind, ArrayRef<Expr> args);

  template <typename IteratedExpr>
  IteratedExpr MakeIterated(expr::Kind kind, ArrayRef<Expr> args) {
    return Expr::Create<IteratedExpr>(MakeIterated(kind, args));
  }

 public:
  typedef mp::Expr Expr;
  typedef mp::NumericExpr NumericExpr;
  typedef mp::LogicalExpr LogicalExpr;
  typedef mp::Variable Variable;
  typedef mp::CountExpr CountExpr;

  explicit ASLBuilder(ASL *asl = 0);
  ~ASLBuilder();

  void set_flags(int flags) { flags_ = flags; }

  // Initializes the ASL object in a similar way to jac0dim, but
  // doesn't read the .nl file as it is the responsibility of NLReader.
  // Instead it uses the information provided in NLHeader.
  void InitASL(const char *stub, const NLHeader &h);

  // Begins building the ASL object.
  // flags: reader flags, see ASL_reader_flag_bits.
  // Throws ASLError on error.
  void BeginBuild(const char *stub, const NLHeader &h);

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

  template <typename Grad>
  class LinearExprHandler {
   private:
    ASLBuilder *builder_;
    Grad **term_;

   public:
    LinearExprHandler(ASLBuilder *b, Grad **term) : builder_(b), term_(term) {}
    ~LinearExprHandler() { *term_ = 0; }

    void AddTerm(int var_index, double coef) {
      Grad *og = builder_->Allocate<Grad>();
      *term_ = og;
      term_ = &og->next;
      og->varno = var_index;
      og->coef = coef;
    }
  };

  typedef LinearExprHandler<ograd> LinearObjHandler;
  typedef LinearExprHandler<cgrad> LinearConHandler;

  LinearConHandler GetLinearVarHandler(int, int) {
    // TODO
    return LinearConHandler(0, 0);
  }
  LinearObjHandler GetLinearObjHandler(int index, int) {
    return LinearObjHandler(this, asl_->i.Ograd_ + index);
  }
  LinearConHandler GetLinearConHandler(int index, int) {
    return LinearConHandler(this, asl_->i.Cgrad_ + index);
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

  UnaryExpr MakeUnary(expr::Kind kind, NumericExpr arg);

  BinaryExpr MakeBinary(expr::Kind kind, NumericExpr lhs, NumericExpr rhs) {
    return MakeBinary<BinaryExpr>(kind, lhs, rhs, "binary");
  }

  IfExpr MakeIf(LogicalExpr condition,
      NumericExpr true_expr, NumericExpr false_expr) {
    return Expr::Create<IfExpr>(
        MakeIf(expr::IF, condition, true_expr, false_expr));
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

  VarArgExpr MakeVarArg(expr::Kind kind, ArrayRef<NumericExpr> args);

  SumExpr MakeSum(ArrayRef<NumericExpr> args) {
    return MakeIterated<SumExpr>(expr::SUM, args);
  }

  CountExpr MakeCount(ArrayRef<LogicalExpr> args) {
    return MakeIterated<CountExpr>(expr::COUNT, args);
  }

  NumberOfExpr MakeNumberOf(ArrayRef<NumericExpr> args) {
    assert(args.size() >= 1);
    return MakeIterated<NumberOfExpr>(expr::NUMBEROF, args);
  }

  LogicalConstant MakeLogicalConstant(bool value) {
    return Expr::Create<LogicalConstant>(MakeConstant(value));
  }

  NotExpr MakeNot(LogicalExpr arg) { return MakeUnary<NotExpr>(expr::NOT, arg); }

  BinaryLogicalExpr MakeBinaryLogical(
      expr::Kind kind, LogicalExpr lhs, LogicalExpr rhs) {
    return MakeBinary<BinaryLogicalExpr>(kind, lhs, rhs, "binary logical");
  }

  RelationalExpr MakeRelational(
      expr::Kind kind, NumericExpr lhs, NumericExpr rhs) {
    return MakeBinary<RelationalExpr>(kind, lhs, rhs, "relational");
  }

  LogicalCountExpr MakeLogicalCount(
      expr::Kind kind, NumericExpr lhs, CountExpr rhs) {
    return MakeBinary<LogicalCountExpr>(kind, lhs, rhs, "logical count");
  }

  ImplicationExpr MakeImplication(
      LogicalExpr condition, LogicalExpr true_expr, LogicalExpr false_expr) {
    return Expr::Create<ImplicationExpr>(
        MakeIf(expr::IMPLICATION, condition, true_expr, false_expr));
  }

  IteratedLogicalExpr MakeIteratedLogical(
      expr::Kind kind, ArrayRef<LogicalExpr> args);

  AllDiffExpr MakeAllDiff(ArrayRef<NumericExpr> args) {
    return MakeIterated<AllDiffExpr>(expr::ALLDIFF, args);
  }

  // Constructs a StringLiteral object.
  // value: string value which may not be null-terminated.
  StringLiteral MakeStringLiteral(fmt::StringRef value);
};
}
}  // namespace mp

#endif  // MP_ASLBUILDER_H_
