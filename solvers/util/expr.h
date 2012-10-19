/*
 A C++ interface to AMPL expression trees.

 Copyright (C) 2012 AMPL Optimization LLC

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization LLC disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Victor Zverovich
 */

#ifndef SOLVERS_UTIL_EXPR_H_
#define SOLVERS_UTIL_EXPR_H_

#include <cassert>
#include <cstddef>
#include <iterator>

#include "solvers/util/util.h"

extern "C" {
#include "solvers/nlp.h"
#include "solvers/opcode.hd"
}

namespace ampl {

// An operation type.
// Operation types should be in sync with the codes in op_type.hd.
enum OpType {
  OPTYPE_UNARY    =  1, // Unary operation
  OPTYPE_BINARY   =  2, // Binary operation
  OPTYPE_VARARG   =  3, // Variable-argument function such as min or max
  OPTYPE_PLTERM   =  4, // Piecewise-linear term
  OPTYPE_IF       =  5, // The if-then-else expression
  OPTYPE_SUM      =  6, // The sum expression
  OPTYPE_FUNCALL  =  7, // Function call
  OPTYPE_STRING   =  8, // String
  OPTYPE_NUMBER   =  9, // Number
  OPTYPE_VARIABLE = 10, // Variable
  OPTYPE_COUNT    = 11  // The count expression
};

// An expression.
// An Expr object represents a handle (reference) to an expression so
// it is cheap to construct and pass by value. A type safe way to
// process expressions of different types is by using ExprVisitor.
class Expr {
 protected:
  const expr *expr_;

  // Returns true iff this expression is null or has type t.
  bool HasTypeOrNull(OpType t) const {
    return !expr_ || optype[opcode()] == t;
  }

  void True() const {}
  typedef void (Expr::*SafeBool)() const;

 public:
  // Constructs a Expr object representing a reference to the AMPL
  // expression e.
  explicit Expr(const expr *e = 0) : expr_(e) {}

  // TODO(viz): remove
  const expr *get() const { return expr_; }

  operator SafeBool() const { return expr_ != 0 ? &Expr::True : 0; }

  // Returns the operation code (opcode) of this expression.
  // The opcodes are defined in opcode.hd.
  unsigned opcode() const {
    return reinterpret_cast<std::size_t>(expr_->op);
  }
};

template <typename T>
T Cast(Expr e);

// An unary expression.
class UnaryExpr : public Expr {
 private:
  explicit UnaryExpr(Expr e) : Expr(e) {
    assert(HasTypeOrNull(OPTYPE_UNARY));
  }

  template <typename Impl, typename Result>
  friend class ExprVisitor;

 public:
  // Returns the argument of this expression.
  Expr arg() const { return Expr(expr_->L.e); }
};

// A binary expression.
class BinaryExpr : public Expr {
 private:
  explicit BinaryExpr(Expr e) : Expr(e) {
    assert(HasTypeOrNull(OPTYPE_BINARY));
  }

  template <typename Impl, typename Result>
  friend class ExprVisitor;

 public:
  // Returns the left-hand side (the first argument) of this expression.
  Expr lhs() const { return Expr(expr_->L.e); }

  // Returns the right-hand side (the second argument) of this expression.
  Expr rhs() const { return Expr(expr_->R.e); }
};

// A variable argument expression.
class VarArgExpr : public Expr {
 private:
  explicit VarArgExpr(Expr e) : Expr(e) {
    assert(HasTypeOrNull(OPTYPE_VARARG));
  }

  template <typename Impl, typename Result>
  friend class ExprVisitor;

  static const de END;

 public:
  // An argument iterator.
  class iterator : public std::iterator<std::forward_iterator_tag, Expr> {
   private:
    const de *de_;

    friend class VarArgExpr;

    explicit iterator(const de *d) : de_(d) {}

   public:
    iterator() : de_(&END) {}

    Expr operator*() const { return Expr(de_->e); }

    iterator &operator++() {
      ++de_;
      return *this;
    }

    iterator operator++(int) {
      iterator it(*this);
      ++it.de_;
      return it;
    }

    bool operator==(iterator other) const { return de_->e == other.de_->e; }
    bool operator!=(iterator other) const { return de_->e != other.de_->e; }
  };

  iterator begin() const {
    return iterator(reinterpret_cast<const expr_va*>(expr_)->L.d);
  }

  iterator end() const {
    return iterator();
  }
};

// An argument iterator.
class ArgIterator : public std::iterator<std::forward_iterator_tag, Expr> {
 private:
  const expr *const *ptr_;

  friend class SumExpr;

 public:
  explicit ArgIterator(const expr *const *p = 0) : ptr_(p) {}

  Expr operator*() const { return Expr(*ptr_); }

  ArgIterator &operator++() {
    ++ptr_;
    return *this;
  }

  ArgIterator operator++(int) {
    ArgIterator it(*this);
    ++it.ptr_;
    return it;
  }

  bool operator==(ArgIterator other) const { return ptr_ == other.ptr_; }
  bool operator!=(ArgIterator other) const { return ptr_ != other.ptr_; }
};

// A sum or count expression.
class SumExpr : public Expr {
 private:
  explicit SumExpr(Expr e) : Expr(e) {
    assert(HasTypeOrNull(OPTYPE_SUM) || HasTypeOrNull(OPTYPE_COUNT));
  }

  template <typename Impl, typename Result>
  friend class ExprVisitor;

 public:
  typedef ArgIterator iterator;

  iterator begin() const {
    return iterator(expr_->L.ep);
  }

  iterator end() const {
    return iterator(expr_->R.ep);
  }
};

// An if-then-else expression.
class IfExpr : public Expr {
 private:
  explicit IfExpr(Expr e) : Expr(e) {
    assert(HasTypeOrNull(OPTYPE_IF));
  }

  template <typename Impl, typename Result>
  friend class ExprVisitor;

 public:
  Expr condition() const {
    return Expr(reinterpret_cast<const expr_if*>(expr_)->e);
  }

  Expr true_expr() const {
    return Expr(reinterpret_cast<const expr_if*>(expr_)->T);
  }

  Expr false_expr() const {
    return Expr(reinterpret_cast<const expr_if*>(expr_)->F);
  }
};

// A piecewise-linear term.
class PiecewiseLinearTerm : public Expr {
 private:
  explicit PiecewiseLinearTerm(Expr e) : Expr(e) {
    assert(HasTypeOrNull(OPTYPE_PLTERM));
  }

  template <typename Impl, typename Result>
  friend class ExprVisitor;

 public:
  int num_slopes() const {
    assert(expr_->L.p->n >= 1);
    return expr_->L.p->n;
  }

  int num_breakpoints() const {
    return num_slopes() - 1;
  }

  double slope(int index) const {
    assert(index >= 0 && index < num_slopes());
    return expr_->L.p->bs[2 * index];
  }

  double breakpoint(int index) const {
    assert(index >= 0 && index < num_breakpoints());
    return expr_->L.p->bs[2 * index + 1];
  }

  int var_index() const {
    return reinterpret_cast<expr_v*>(expr_->R.e)->a;
  }
};

// A number.
class Number : public Expr {
 private:
  explicit Number(Expr e) : Expr(e) {
    assert(HasTypeOrNull(OPTYPE_NUMBER));
  }

  template <typename Impl, typename Result>
  friend class ExprVisitor;

  friend Number Cast<Number>(Expr);

 public:
  // Returns the value of this number.
  double value() const { return reinterpret_cast<const expr_n*>(expr_)->v; }
};

template <>
inline Number Cast<Number>(Expr e) {
  return Number(e.opcode() == OPNUM ? e : Expr());
}

// A variable.
class Variable : public Expr {
 private:
  explicit Variable(Expr e) : Expr(e) {
    assert(HasTypeOrNull(OPTYPE_VARIABLE));
  }

  template <typename Impl, typename Result>
  friend class ExprVisitor;

 public:
  // Returns the index of this variable.
  int index() const { return expr_->a; }
};

// A numberof expression.
class NumberOfExpr : public Expr {
 private:
  explicit NumberOfExpr(Expr e) : Expr(e) {
    assert(opcode() == OPNUMBEROF);
  }

  template <typename Impl, typename Result>
  friend class ExprVisitor;

 public:
  Expr target() const { return Expr(*expr_->L.ep); }

  typedef ArgIterator iterator;

  iterator begin() const {
    return iterator(expr_->L.ep + 1);
  }

  iterator end() const {
    return iterator(expr_->R.ep);
  }
};

// An expression visitor.
// To use ExprVisitor define a subclass that implements some or all of the
// Visit* methods with the same signatures as the methods in ExprVisitor,
// for example, VisitDiv(BinaryExpr).
// Specify the subclass name as the Impl template parameter. Then calling
// ExprVisitor::Visit for some expression will dispatch to a Visit* method
// specific to the expression type. For example, if the expression is
// a division then VisitDiv(BinaryExpr) method of a subclass will be called.
// If the subclass doesn't contain a method with this signature, then
// a corresponding method of ExprVisitor will be called.
//
// Example:
//  class MyExprVisitor : public ExprVisitor<MyExprVisitor, double> {
//   public:
//    double VisitPlus(BinaryExpr e) { return Visit(e.lhs()) + Visit(e.rhs()); }
//    double VisitNumber(Number n) { return n.value(); }
//  };
//
// ExprVisitor uses the curiously recurring template pattern:
// http://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
template <typename Impl, typename Result>
class ExprVisitor {
 public:
  Result Visit(Expr e);

  Result VisitUnhandledExpr(Expr e) {
    throw UnsupportedExprError(GetOpName(e.opcode()));
  }

  Result VisitPlus(BinaryExpr e) { return VisitUnhandledExpr(e); }
  Result VisitMinus(BinaryExpr e) { return VisitUnhandledExpr(e); }
  Result VisitMult(BinaryExpr e) { return VisitUnhandledExpr(e); }
  Result VisitDiv(BinaryExpr e) { return VisitUnhandledExpr(e); }
  Result VisitRem(BinaryExpr e) { return VisitUnhandledExpr(e); }
  Result VisitPow(BinaryExpr e) { return VisitUnhandledExpr(e); }
  Result VisitLess(BinaryExpr e) { return VisitUnhandledExpr(e); }
  Result VisitMin(VarArgExpr e) { return VisitUnhandledExpr(e); }
  Result VisitMax(VarArgExpr e) { return VisitUnhandledExpr(e); }
  Result VisitFloor(UnaryExpr e) { return VisitUnhandledExpr(e); }
  Result VisitCeil(UnaryExpr e) { return VisitUnhandledExpr(e); }
  Result VisitAbs(UnaryExpr e) { return VisitUnhandledExpr(e); }
  Result VisitMinus(UnaryExpr e) { return VisitUnhandledExpr(e); }
  Result VisitIf(IfExpr e) { return VisitUnhandledExpr(e); }
  Result VisitTanh(UnaryExpr e) { return VisitUnhandledExpr(e); }
  Result VisitTan(UnaryExpr e) { return VisitUnhandledExpr(e); }
  Result VisitSqrt(UnaryExpr e) { return VisitUnhandledExpr(e); }
  Result VisitSinh(UnaryExpr e) { return VisitUnhandledExpr(e); }
  Result VisitSin(UnaryExpr e) { return VisitUnhandledExpr(e); }
  Result VisitLog10(UnaryExpr e) { return VisitUnhandledExpr(e); }
  Result VisitLog(UnaryExpr e) { return VisitUnhandledExpr(e); }
  Result VisitExp(UnaryExpr e) { return VisitUnhandledExpr(e); }
  Result VisitCosh(UnaryExpr e) { return VisitUnhandledExpr(e); }
  Result VisitCos(UnaryExpr e) { return VisitUnhandledExpr(e); }
  Result VisitAtanh(UnaryExpr e) { return VisitUnhandledExpr(e); }
  Result VisitAtan2(BinaryExpr e) { return VisitUnhandledExpr(e); }
  Result VisitAtan(UnaryExpr e) { return VisitUnhandledExpr(e); }
  Result VisitAsinh(UnaryExpr e) { return VisitUnhandledExpr(e); }
  Result VisitAsin(UnaryExpr e) { return VisitUnhandledExpr(e); }
  Result VisitAcosh(UnaryExpr e) { return VisitUnhandledExpr(e); }
  Result VisitAcos(UnaryExpr e) { return VisitUnhandledExpr(e); }
  Result VisitSum(SumExpr e) { return VisitUnhandledExpr(e); }
  Result VisitIntDiv(BinaryExpr e) { return VisitUnhandledExpr(e); }
  Result VisitPrecision(BinaryExpr e) { return VisitUnhandledExpr(e); }
  Result VisitRound(BinaryExpr e) { return VisitUnhandledExpr(e); }
  Result VisitTrunc(BinaryExpr e) { return VisitUnhandledExpr(e); }
  Result VisitCount(SumExpr e) { return VisitUnhandledExpr(e); }
  Result VisitNumberOf(NumberOfExpr e) { return VisitUnhandledExpr(e); }
  Result VisitConstExpPow(BinaryExpr e) { return VisitUnhandledExpr(e); }
  Result VisitPow2(UnaryExpr e) { return VisitUnhandledExpr(e); }
  Result VisitConstBasePow(BinaryExpr e) { return VisitUnhandledExpr(e); }
  Result VisitPLTerm(PiecewiseLinearTerm t) { return VisitUnhandledExpr(t); }
  Result VisitNumber(Number n) { return VisitUnhandledExpr(n); }
  Result VisitVariable(Variable v) { return VisitUnhandledExpr(v); }
};

#define AMPL_DISPATCH(call) static_cast<Impl*>(this)->call

template <typename Impl, typename Result>
Result ExprVisitor<Impl, Result>::Visit(Expr e) {
  // TODO(viz): check that all opcodes are handled
  switch (e.opcode()) {
  case OPPLUS:
    return AMPL_DISPATCH(VisitPlus(BinaryExpr(e)));
  case OPMINUS:
    return AMPL_DISPATCH(VisitMinus(BinaryExpr(e)));
  case OPMULT:
    return AMPL_DISPATCH(VisitMult(BinaryExpr(e)));
  case OPDIV:
    return AMPL_DISPATCH(VisitDiv(BinaryExpr(e)));
  case OPREM:
    return AMPL_DISPATCH(VisitRem(BinaryExpr(e)));
  case OPPOW:
    return AMPL_DISPATCH(VisitPow(BinaryExpr(e)));
  case OPLESS:
    return AMPL_DISPATCH(VisitLess(BinaryExpr(e)));
  case MINLIST:
    return AMPL_DISPATCH(VisitMin(VarArgExpr(e)));
  case MAXLIST:
    return AMPL_DISPATCH(VisitMax(VarArgExpr(e)));
  case FLOOR:
    return AMPL_DISPATCH(VisitFloor(UnaryExpr(e)));
  case CEIL:
    return AMPL_DISPATCH(VisitCeil(UnaryExpr(e)));
  case ABS:
    return AMPL_DISPATCH(VisitAbs(UnaryExpr(e)));
  case OPUMINUS:
    return AMPL_DISPATCH(VisitMinus(UnaryExpr(e)));
  case OPIFnl:
    return AMPL_DISPATCH(VisitIf(IfExpr(e)));
  case OP_tanh:
    return AMPL_DISPATCH(VisitTanh(UnaryExpr(e)));
  case OP_tan:
    return AMPL_DISPATCH(VisitTan(UnaryExpr(e)));
  case OP_sqrt:
    return AMPL_DISPATCH(VisitSqrt(UnaryExpr(e)));
  case OP_sinh:
    return AMPL_DISPATCH(VisitSinh(UnaryExpr(e)));
  case OP_sin:
    return AMPL_DISPATCH(VisitSin(UnaryExpr(e)));
  case OP_log10:
    return AMPL_DISPATCH(VisitLog10(UnaryExpr(e)));
  case OP_log:
    return AMPL_DISPATCH(VisitLog(UnaryExpr(e)));
  case OP_exp:
    return AMPL_DISPATCH(VisitExp(UnaryExpr(e)));
  case OP_cosh:
    return AMPL_DISPATCH(VisitCosh(UnaryExpr(e)));
  case OP_cos:
    return AMPL_DISPATCH(VisitCos(UnaryExpr(e)));
  case OP_atanh:
    return AMPL_DISPATCH(VisitAtanh(UnaryExpr(e)));
  case OP_atan2:
    return AMPL_DISPATCH(VisitAtan2(BinaryExpr(e)));
  case OP_atan:
    return AMPL_DISPATCH(VisitAtan(UnaryExpr(e)));
  case OP_asinh:
    return AMPL_DISPATCH(VisitAsinh(UnaryExpr(e)));
  case OP_asin:
    return AMPL_DISPATCH(VisitAsin(UnaryExpr(e)));
  case OP_acosh:
    return AMPL_DISPATCH(VisitAcosh(UnaryExpr(e)));
  case OP_acos:
    return AMPL_DISPATCH(VisitAcos(UnaryExpr(e)));
  case OPSUMLIST:
    return AMPL_DISPATCH(VisitSum(SumExpr(e)));
  case OPintDIV:
    return AMPL_DISPATCH(VisitIntDiv(BinaryExpr(e)));
  case OPprecision:
    return AMPL_DISPATCH(VisitPrecision(BinaryExpr(e)));
  case OPround:
    return AMPL_DISPATCH(VisitRound(BinaryExpr(e)));
  case OPtrunc:
    return AMPL_DISPATCH(VisitTrunc(BinaryExpr(e)));
  case OP1POW:
    return AMPL_DISPATCH(VisitConstExpPow(BinaryExpr(e)));
  case OP2POW:
    return AMPL_DISPATCH(VisitPow2(UnaryExpr(e)));
  case OPCPOW:
    return AMPL_DISPATCH(VisitConstBasePow(BinaryExpr(e)));
  case OPPLTERM:
    return AMPL_DISPATCH(VisitPLTerm(PiecewiseLinearTerm(e)));
  case OPNUM:
    return AMPL_DISPATCH(VisitNumber(Number(e)));
  case OPVARVAL:
    return AMPL_DISPATCH(VisitVariable(Variable(e)));
  case OPCOUNT:
    return AMPL_DISPATCH(VisitCount(SumExpr(e)));
  case OPNUMBEROF:
    return AMPL_DISPATCH(VisitNumberOf(NumberOfExpr(e)));
  default:
    return VisitUnhandledExpr(e);
  }
}

#undef AMPL_DISPATCH
}

#endif  // SOLVERS_UTIL_EXPR_H_
