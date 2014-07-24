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

  template <typename ExprT>
  ExprT MakeExpr(int opcode, Expr lhs, Expr rhs);

  template <typename T>
  T *Allocate(unsigned size = sizeof(T)) {
    assert(size >= sizeof(T));
    return reinterpret_cast<T*>(mem_ASL(asl_, size));
  }

  template <typename ExprT>
  ExprT MakeConstant(double value);

  template <Expr::Kind KIND, typename Arg>
  BasicBinaryExpr<KIND, Arg> MakeBinary(int opcode, Arg lhs, Arg rhs) {
    CheckOpCode(opcode, KIND, "binary");
    typedef BasicBinaryExpr<KIND, Arg> BinaryExpr;
    BinaryExpr expr = MakeExpr<BinaryExpr>(opcode, lhs, rhs);
    expr.expr_->dL = 1;
    expr.expr_->dR = DVALUE[opcode];  // for PLUS, MINUS, REM
    return expr;
  }

  // Makes an iterated expression.
  template <typename ExprT>
  ExprT MakeIterated(int opcode, int num_args);

  // Makes an iterated expression.
  template <Expr::Kind K>
  BasicIteratedExpr<K> MakeIterated(int opcode,
      int num_args, typename IteratedExprInfo<K>::Arg *args);

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

  // The Make* methods construct expression objects. These objects are
  // local to the currently built ASL problem and shouldn't be used with
  // other problems. The expression objects are not accessible via the
  // ASL API until they are added to the problem as a part of objective
  // or constraint expression. For this reason the methods below use a
  // different naming convention from the Add* methods.

  UnaryExpr MakeUnary(int opcode, NumericExpr arg);

  BinaryExpr MakeBinary(int opcode, NumericExpr lhs, NumericExpr rhs) {
    return MakeBinary<Expr::BINARY, NumericExpr>(opcode, lhs, rhs);
  }

  VarArgExpr MakeVarArg(int opcode, int num_args, NumericExpr *args);

  SumExpr MakeSum(int num_args, NumericExpr *args) {
    return MakeIterated<Expr::SUM>(OPSUMLIST, num_args, args);
  }

  CountExpr MakeCount(int num_args, LogicalExpr *args) {
    return MakeIterated<Expr::COUNT>(OPCOUNT, num_args, args);
  }

  IfExpr MakeIf(LogicalExpr condition,
      NumericExpr true_expr, NumericExpr false_expr);

  PiecewiseLinearExpr MakePiecewiseLinear(int num_breakpoints,
      const double *breakpoints, const double *slopes, Variable var);

  Variable MakeVariable(int var_index);

  NumberOfExpr MakeNumberOf(
      NumericExpr value, int num_args, const NumericExpr *args);

  CallExpr MakeCall(Function f, int num_args, const Expr *args);

  NumericConstant MakeNumericConstant(double value) {
    return MakeConstant<NumericConstant>(value);
  }

  LogicalConstant MakeLogicalConstant(bool value) {
    return MakeConstant<LogicalConstant>(value);
  }

  BinaryLogicalExpr MakeBinaryLogical(
      int opcode, LogicalExpr lhs, LogicalExpr rhs) {
    return MakeBinary<Expr::BINARY_LOGICAL, LogicalExpr>(opcode, lhs, rhs);
  }

  StringLiteral MakeStringLiteral(int size, const char *value);
};

template <typename ExprT>
ExprT ASLBuilder::MakeExpr(int opcode, Expr lhs, Expr rhs) {
  expr *e = Allocate<expr>();
  e->op = reinterpret_cast<efunc*>(opcode);
  e->L.e = lhs.expr_;
  e->R.e = rhs.expr_;
  e->a = asl_->i.n_var_ + asl_->i.nsufext[ASL_Sufkind_var];
  e->dL = DVALUE[opcode];  // for UMINUS, FLOOR, CEIL
  return Expr::Create<ExprT>(e);
}

template <typename ExprT>
ExprT ASLBuilder::MakeConstant(double value) {
  expr_n *result = Allocate<expr_n>(asl_->i.size_expr_n_);
  result->op = reinterpret_cast<efunc_n*>(OPNUM);
  result->v = value;
  return Expr::Create<ExprT>(reinterpret_cast<expr*>(result));
}

template <typename ExprT>
ExprT ASLBuilder::MakeIterated(int opcode, int num_args) {
  assert(num_args >= 0);
  expr *result = Allocate<expr>(
      sizeof(expr) - sizeof(double) + (num_args + 1) * sizeof(expr*));
  result->op = reinterpret_cast<efunc*>(opcode);
  result->L.ep = reinterpret_cast<expr**>(&result->dR);
  result->R.ep = result->L.ep + num_args;
  return Expr::Create<ExprT>(result);
}

template <Expr::Kind K>
BasicIteratedExpr<K> ASLBuilder::MakeIterated(
    int opcode, int num_args, typename IteratedExprInfo<K>::Arg *args) {
  BasicIteratedExpr<K> result =
      MakeIterated< BasicIteratedExpr<K> >(opcode, num_args);
  expr **arg_ptrs = result.expr_->L.ep;
  for (int i = 0; i < num_args; ++i)
    arg_ptrs[i] = args[i].expr_;
  return result;
}
}
}

#endif  // SOLVERS_UTIL_ASLBUILDER_H_
