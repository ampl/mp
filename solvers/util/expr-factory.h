/*
 An AMPL expression factory.

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

#ifndef SOLVERS_UTIL_EXPR_FACTORY_H_
#define SOLVERS_UTIL_EXPR_FACTORY_H_

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
      Error(message), error_code_(error_code) {
  }

  int error_code() const { return error_code_; }
};

class ASLBuilder {
 private:
  ASL &asl_;
  efunc **r_ops_;
  int nv1_;
  int nz_;
  int nderp_;

 public:
  ASLBuilder(ASL &asl, const char *stub, const NLHeader &h);

  // Begin building the ASL.
  // flags: reader flags, see ASL_reader_flag_bits.
  // Throws ASLError on error.
  void BeginBuild(int flags);

  // End building the ASL.
  void EndBuild();
};
}

class ExprFactory : Noncopyable {
 private:
  ASL *asl_;
  efunc *r_ops_[N_OPS];

  static void CheckOpCode(int opcode, Expr::Kind kind, const char *expr_name) {
    if (Expr::INFO[opcode].kind != kind)
      throw Error("invalid {} expression code {}", expr_name, opcode);
  }

  template <typename ExprT>
  ExprT MakeExpr(int opcode, NumericExpr lhs, NumericExpr rhs);

  template <typename T>
  T *Allocate(unsigned size = sizeof(T)) {
    assert(size >= sizeof(T));
    return reinterpret_cast<T*>(mem_ASL(asl_, size));
  }

  template <typename ExprT>
  ExprT MakeConstant(double value);

  // Make sum or count expression.
  template <typename Arg>
  BasicSumExpr<Arg> MakeSum(int opcode, int num_args, Arg *args);

 public:
  // Constructs an ExprFactory object.
  // flags: reader flags, see ASL_reader_flag_bits.
  ExprFactory(const NLHeader &h, const char *stub, int flags = 0);
  ~ExprFactory();

  UnaryExpr MakeUnary(int opcode, NumericExpr arg);
  BinaryExpr MakeBinary(int opcode, NumericExpr lhs, NumericExpr rhs);
  VarArgExpr MakeVarArg(int opcode, int num_args, NumericExpr *args);

  SumExpr MakeSum(int num_args, NumericExpr *args) {
    return MakeSum<NumericExpr>(OPSUMLIST, num_args, args);
  }

  CountExpr MakeCount(int num_args, LogicalExpr *args) {
    return MakeSum<LogicalExpr>(OPCOUNT, num_args, args);
  }

  IfExpr MakeIf(LogicalExpr condition,
      NumericExpr true_expr, NumericExpr false_expr);

  PiecewiseLinearExpr MakePiecewiseLinear(int num_breakpoints,
      const double *breakpoints, const double *slopes, Variable var);

  NumericConstant MakeNumericConstant(double value) {
    return MakeConstant<NumericConstant>(value);
  }

  Variable MakeVariable(int var_index);

  LogicalConstant MakeLogicalConstant(bool value) {
    return MakeConstant<LogicalConstant>(value);
  }
};

template <typename ExprT>
ExprT ExprFactory::MakeConstant(double value) {
  expr_n *result = Allocate<expr_n>(asl_->i.size_expr_n_);
  result->op = reinterpret_cast<efunc_n*>(r_ops_[OPNUM]);
  result->v = value;
  return Expr::Create<ExprT>(reinterpret_cast<expr*>(result));
}

template <typename Arg>
BasicSumExpr<Arg> ExprFactory::MakeSum(int opcode, int num_args, Arg *args) {
  assert(num_args >= 0);
  expr *result = Allocate<expr>(
      sizeof(expr) - sizeof(double) + num_args * sizeof(expr*));
  result->op = r_ops_[opcode];
  expr **arg_ptrs = result->L.ep = reinterpret_cast<expr**>(&result->dR);
  for (int i = 0; i < num_args; ++i)
    arg_ptrs[i] = args[i].expr_;
  result->R.ep = arg_ptrs + num_args;
  return Expr::Create< BasicSumExpr<Arg> >(result);
}
}

#endif  // SOLVERS_UTIL_EXPR_FACTORY_H_
