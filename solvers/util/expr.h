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
#include <stdexcept>

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

// A base class for all expression classes.
class ExprBase {
 protected:
  expr *expr_;

  // Returns true iff this expression is null or has type t.
  bool HasTypeOrNull(OpType t) const {
    return !expr_ || type() == t;
  }

  void True() const {}
  typedef void (ExprBase::*SafeBool)() const;

  // Constructs a ExprBase object representing a reference to the AMPL
  // expression e.
  explicit ExprBase(expr *e = 0) : expr_(e) {}

 public:
  // TODO(viz): remove
  expr *get() const { return expr_; }
  
  operator SafeBool() const { return expr_ != 0 ? &ExprBase::True : 0; }

  // Returns the operation code (opcode) of this expression.
  // The opcodes are defined in opcode.hd.
  unsigned opcode() const {
    return reinterpret_cast<std::size_t>(expr_->op);
  }
  
  // Returns the operation name of this expression.
  const char *opname() const;

  // Returns the type of this expression.
  unsigned type() const { return optype[opcode()]; }
};

// An arithmetic expression.
// An Expr object represents a handle (reference) to an expression so
// it is cheap to construct and pass by value. A type safe way to
// process expressions of different types is by using ExprVisitor.
class Expr : public ExprBase {
 public:
  // Constructs a Expr object representing a reference to the AMPL
  // expression e.
  explicit Expr(expr *e = 0) : ExprBase(e) {}
};

// Returns true if the expressions e1 and e2 are structurally equal.
bool Equal(Expr e1, Expr e2);

// A logical or constraint expression.
class LogicalExpr : public ExprBase {
 public:
  explicit LogicalExpr(expr *e = 0) : ExprBase(e) {}
};

template <typename T>
T Cast(Expr e);

template <typename ExprT>
class ExprProxy {
 private:
  ExprT expr_;

 public:
  explicit ExprProxy(expr *e) : expr_(e) {}

  ExprT *operator->() const { return &expr_; }
};

// An unary expression.
class UnaryExpr : public Expr {
 private:
  explicit UnaryExpr(Expr e) : Expr(e) {
    assert(HasTypeOrNull(OPTYPE_UNARY));
  }

  template <typename Impl, typename Result, typename LResult>
  friend class ExprVisitor;

 public:
  // Returns the argument of this expression.
  Expr arg() const { return Expr(expr_->L.e); }
};

// A binary expression.
class BinaryExpr : public Expr {
 private:
  explicit BinaryExpr(Expr e) : Expr(e) {
    assert(!expr_ || opcode() == OP1POW || opcode() == OPCPOW ||
           HasTypeOrNull(OPTYPE_BINARY));
  }

  template <typename Impl, typename Result, typename LResult>
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

  template <typename Impl, typename Result, typename LResult>
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
    ExprProxy<Expr> operator->() const { return ExprProxy<Expr>(de_->e); }

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
    return iterator(reinterpret_cast<expr_va*>(expr_)->L.d);
  }

  iterator end() const {
    return iterator();
  }
};

// An argument iterator.
template <typename ExprT>
class ArgIterator : public std::iterator<std::forward_iterator_tag, ExprT> {
 private:
  expr *const *ptr_;

 public:
  explicit ArgIterator(expr *const *p = 0) : ptr_(p) {}

  ExprT operator*() const { return ExprT(*ptr_); }
  ExprProxy<ExprT> operator->() const { return ExprProxy<ExprT>(*ptr_); }

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

// An sum expression.
class SumExpr : public Expr {
 private:
  explicit SumExpr(Expr e) : Expr(e) {
    assert(HasTypeOrNull(OPTYPE_SUM));
  }

  template <typename Impl, typename Result, typename LResult>
  friend class ExprVisitor;

 public:
  typedef ArgIterator<Expr> iterator;

  iterator begin() const {
    return iterator(expr_->L.ep);
  }

  iterator end() const {
    return iterator(expr_->R.ep);
  }
};

// A count expression.
class CountExpr : public Expr {
 private:
  explicit CountExpr(Expr e) : Expr(e) {
    assert(HasTypeOrNull(OPTYPE_COUNT));
  }

  template <typename Impl, typename Result, typename LResult>
  friend class ExprVisitor;

 public:
  typedef ArgIterator<LogicalExpr> iterator;

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

  template <typename Impl, typename Result, typename LResult>
  friend class ExprVisitor;

 public:
  LogicalExpr condition() const {
    return LogicalExpr(reinterpret_cast<expr_if*>(expr_)->e);
  }

  Expr true_expr() const {
    return Expr(reinterpret_cast<expr_if*>(expr_)->T);
  }

  Expr false_expr() const {
    return Expr(reinterpret_cast<expr_if*>(expr_)->F);
  }
};

// A piecewise-linear term.
class PiecewiseLinearTerm : public Expr {
 private:
  explicit PiecewiseLinearTerm(Expr e) : Expr(e) {
    assert(HasTypeOrNull(OPTYPE_PLTERM));
  }

  template <typename Impl, typename Result, typename LResult>
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

// A numeric constant.
class NumericConstant : public Expr {
 private:
  explicit NumericConstant(Expr e) : Expr(e) {
    assert(HasTypeOrNull(OPTYPE_NUMBER));
  }

  template <typename Impl, typename Result, typename LResult>
  friend class ExprVisitor;

  friend NumericConstant Cast<NumericConstant>(Expr);

 public:
  // Returns the value of this number.
  double value() const { return reinterpret_cast<expr_n*>(expr_)->v; }
};

template <>
inline NumericConstant Cast<NumericConstant>(Expr e) {
  return NumericConstant(e.opcode() == OPNUM ? e : Expr());
}

// A variable.
class Variable : public Expr {
 private:
  explicit Variable(Expr e) : Expr(e) {
    assert(HasTypeOrNull(OPTYPE_VARIABLE));
  }

  template <typename Impl, typename Result, typename LResult>
  friend class ExprVisitor;

  friend Variable Cast<Variable>(Expr);

 public:
  // Returns the index of this variable.
  int index() const { return expr_->a; }
};

template <>
inline Variable Cast<Variable>(Expr e) {
  return Variable(e.opcode() == OPVARVAL ? e : Expr());
}

// A numberof expression.
class NumberOfExpr : public Expr {
 private:
  explicit NumberOfExpr(Expr e) : Expr(e) {
    assert(opcode() == OPNUMBEROF);
  }

  template <typename Impl, typename Result, typename LResult>
  friend class ExprVisitor;

 public:
  Expr target() const { return Expr(*expr_->L.ep); }

  typedef ArgIterator<Expr> iterator;

  iterator begin() const {
    return iterator(expr_->L.ep + 1);
  }

  iterator end() const {
    return iterator(expr_->R.ep);
  }
};

class LogicalConstant : public LogicalExpr {
 private:
  explicit LogicalConstant(LogicalExpr e) : LogicalExpr(e) {
    assert(HasTypeOrNull(OPTYPE_NUMBER));
  }

  template <typename Impl, typename Result, typename LResult>
  friend class ExprVisitor;

 public:
  // Returns the value of this constant.
  bool value() const { return reinterpret_cast<expr_n*>(expr_)->v != 0; }
};

// A relational expression such as "less than".
class RelationalExpr : public LogicalExpr {
 private:
  explicit RelationalExpr(LogicalExpr e) : LogicalExpr(e) {
    assert(HasTypeOrNull(OPTYPE_BINARY));
  }

  template <typename Impl, typename Result, typename LResult>
  friend class ExprVisitor;

 public:
  // Returns the left-hand side (the first argument) of this expression.
  Expr lhs() const { return Expr(expr_->L.e); }

  // Returns the right-hand side (the second argument) of this expression.
  Expr rhs() const { return Expr(expr_->R.e); }
};

// A logical NOT expression.
class NotExpr : public LogicalExpr {
 private:
  explicit NotExpr(LogicalExpr e) : LogicalExpr(e) {
    assert(opcode() == OPNOT);
  }

  template <typename Impl, typename Result, typename LResult>
  friend class ExprVisitor;

 public:
  // Returns the argument of this expression.
  LogicalExpr arg() const { return LogicalExpr(expr_->L.e); }
};

// A binary logical expression such as logical OR.
class BinaryLogicalExpr : public LogicalExpr {
 private:
  explicit BinaryLogicalExpr(LogicalExpr e) : LogicalExpr(e) {
    assert(HasTypeOrNull(OPTYPE_BINARY));
  }

  template <typename Impl, typename Result, typename LResult>
  friend class ExprVisitor;

 public:
  // Returns the left-hand side (the first argument) of this expression.
  LogicalExpr lhs() const { return LogicalExpr(expr_->L.e); }

  // Returns the right-hand side (the second argument) of this expression.
  LogicalExpr rhs() const { return LogicalExpr(expr_->R.e); }
};

// An implication expression.
class ImplicationExpr : public LogicalExpr {
 private:
  explicit ImplicationExpr(LogicalExpr e) : LogicalExpr(e) {
    assert(HasTypeOrNull(OPTYPE_IF));
  }

  template <typename Impl, typename Result, typename LResult>
  friend class ExprVisitor;

 public:
  LogicalExpr condition() const {
    return LogicalExpr(reinterpret_cast<expr_if*>(expr_)->e);
  }

  LogicalExpr true_expr() const {
    return LogicalExpr(reinterpret_cast<expr_if*>(expr_)->T);
  }

  LogicalExpr false_expr() const {
    return LogicalExpr(reinterpret_cast<expr_if*>(expr_)->F);
  }
};

// An iterated logical expression such as exists or forall.
class IteratedLogicalExpr : public LogicalExpr {
 private:
  explicit IteratedLogicalExpr(LogicalExpr e) : LogicalExpr(e) {
    assert(HasTypeOrNull(OPTYPE_SUM) || HasTypeOrNull(OPTYPE_COUNT));
  }

  template <typename Impl, typename Result, typename LResult>
  friend class ExprVisitor;

 public:
  typedef ArgIterator<LogicalExpr> iterator;

  iterator begin() const {
    return iterator(expr_->L.ep);
  }

  iterator end() const {
    return iterator(expr_->R.ep);
  }
};

// An alldiff expression.
class AllDiffExpr : public LogicalExpr {
 private:
  explicit AllDiffExpr(LogicalExpr e) : LogicalExpr(e) {
    assert(HasTypeOrNull(OPTYPE_COUNT));
  }

  template <typename Impl, typename Result, typename LResult>
  friend class ExprVisitor;

 public:
  typedef ArgIterator<Expr> iterator;

  iterator begin() const {
    return iterator(expr_->L.ep);
  }

  iterator end() const {
    return iterator(expr_->R.ep);
  }
};

/// A general solver error.
class Error : public std::runtime_error {
public:
  explicit Error(const std::string &message) : std::runtime_error(message) {}
};

/// An exception that is thrown when an ASL expression not supported
/// by the solver is encountered.
class UnsupportedExprError : public Error {
public:
  explicit UnsupportedExprError(const char *expr) :
  Error(std::string("unsupported expression: ") + expr) {}
};

/// An exception that is thrown when an invalid arithmetic expression
/// is encountered.
class InvalidExprError : public Error {
public:
  explicit InvalidExprError(const char *expr) :
  Error(std::string("invalid arithmetic expression: ") + expr) {}
};

/// An exception that is thrown when an invalid logical or constraint
/// expression is encountered.
class InvalidLogicalExprError : public Error {
public:
  explicit InvalidLogicalExprError(const char *expr) :
  Error(std::string("invalid logical expression: ") + expr) {}
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
template <typename Impl, typename Result, typename LResult>
class ExprVisitor {
 public:
  Result Visit(Expr e);
  LResult Visit(LogicalExpr e);

  Result VisitUnhandledExpr(Expr e) {
    throw UnsupportedExprError(e.opname());
  }

  LResult VisitUnhandledExpr(LogicalExpr e) {
    throw UnsupportedExprError(e.opname());
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
  Result VisitCount(CountExpr e) { return VisitUnhandledExpr(e); }
  Result VisitNumberOf(NumberOfExpr e) { return VisitUnhandledExpr(e); }
  Result VisitConstExpPow(BinaryExpr e) { return VisitUnhandledExpr(e); }
  Result VisitPow2(UnaryExpr e) { return VisitUnhandledExpr(e); }
  Result VisitConstBasePow(BinaryExpr e) { return VisitUnhandledExpr(e); }
  Result VisitPLTerm(PiecewiseLinearTerm t) { return VisitUnhandledExpr(t); }
  Result VisitConstant(NumericConstant c) { return VisitUnhandledExpr(c); }
  Result VisitVariable(Variable v) { return VisitUnhandledExpr(v); }

  LResult VisitConstant(LogicalConstant c) { return VisitUnhandledExpr(c); }
  LResult VisitLess(RelationalExpr e) { return VisitUnhandledExpr(e); }
  LResult VisitLessEqual(RelationalExpr e) { return VisitUnhandledExpr(e); }
  LResult VisitEqual(RelationalExpr e) { return VisitUnhandledExpr(e); }
  LResult VisitGreaterEqual(RelationalExpr e) { return VisitUnhandledExpr(e); }
  LResult VisitGreater(RelationalExpr e) { return VisitUnhandledExpr(e); }
  LResult VisitNotEqual(RelationalExpr e) { return VisitUnhandledExpr(e); }
  LResult VisitAtMost(RelationalExpr e) { return VisitUnhandledExpr(e); }
  LResult VisitNotAtMost(RelationalExpr e) { return VisitUnhandledExpr(e); }
  LResult VisitAtLeast(RelationalExpr e) { return VisitUnhandledExpr(e); }
  LResult VisitNotAtLeast(RelationalExpr e) { return VisitUnhandledExpr(e); }
  LResult VisitExactly(RelationalExpr e) { return VisitUnhandledExpr(e); }
  LResult VisitNotExactly(RelationalExpr e) { return VisitUnhandledExpr(e); }
  LResult VisitOr(BinaryLogicalExpr e) { return VisitUnhandledExpr(e); }
  LResult VisitExists(IteratedLogicalExpr e) { return VisitUnhandledExpr(e); }
  LResult VisitAnd(BinaryLogicalExpr e) { return VisitUnhandledExpr(e); }
  LResult VisitForAll(IteratedLogicalExpr e) { return VisitUnhandledExpr(e); }
  LResult VisitNot(NotExpr e) { return VisitUnhandledExpr(e); }
  LResult VisitIff(BinaryLogicalExpr e) { return VisitUnhandledExpr(e); }
  LResult VisitImplication(ImplicationExpr e) { return VisitUnhandledExpr(e); }
  LResult VisitAllDiff(AllDiffExpr e) { return VisitUnhandledExpr(e); }
};

#define AMPL_DISPATCH(call) static_cast<Impl*>(this)->call

template <typename Impl, typename Result, typename LResult>
Result ExprVisitor<Impl, Result, LResult>::Visit(Expr e) {
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
    return AMPL_DISPATCH(VisitNumber(NumericConstant(e)));
  case OPVARVAL:
    return AMPL_DISPATCH(VisitVariable(Variable(e)));
  case OPCOUNT:
    return AMPL_DISPATCH(VisitCount(CountExpr(e)));
  case OPNUMBEROF:
    return AMPL_DISPATCH(VisitNumberOf(NumberOfExpr(e)));
  default:
    throw InvalidExprError(e.opname());
  }
}

template <typename Impl, typename Result, typename LResult>
LResult ExprVisitor<Impl, Result, LResult>::Visit(LogicalExpr e) {
  // TODO(viz): check that all opcodes are handled
  switch (e.opcode()) {
  case OPNUM:
    return AMPL_DISPATCH(VisitConstant(LogicalConstant(e)));
  case LT:
    return AMPL_DISPATCH(VisitLess(RelationalExpr(e)));
  case LE:
    return AMPL_DISPATCH(VisitLessEqual(RelationalExpr(e)));
  case EQ:
    return AMPL_DISPATCH(VisitEqual(RelationalExpr(e)));
  case GE:
    return AMPL_DISPATCH(VisitGreaterEqual(RelationalExpr(e)));
  case GT:
    return AMPL_DISPATCH(VisitGreater(RelationalExpr(e)));
  case NE:
    return AMPL_DISPATCH(VisitNotEqual(RelationalExpr(e)));
  case OPATMOST:
    return AMPL_DISPATCH(VisitAtMost(RelationalExpr(e)));
  case OPNOTATMOST:
    return AMPL_DISPATCH(VisitNotAtMost(RelationalExpr(e)));
  case OPATLEAST:
    return AMPL_DISPATCH(VisitAtLeast(RelationalExpr(e)));
  case OPNOTATLEAST:
    return AMPL_DISPATCH(VisitNotAtLeast(RelationalExpr(e)));
  case OPEXACTLY:
    return AMPL_DISPATCH(VisitExactly(RelationalExpr(e)));
  case OPNOTEXACTLY:
    return AMPL_DISPATCH(VisitNotExactly(RelationalExpr(e)));
  case OPOR:
    return AMPL_DISPATCH(VisitOr(BinaryLogicalExpr(e)));
  case ORLIST:
    return AMPL_DISPATCH(VisitExists(IteratedLogicalExpr(e)));
  case OPAND:
    return AMPL_DISPATCH(VisitAnd(BinaryLogicalExpr(e)));
  case ANDLIST:
    return AMPL_DISPATCH(VisitForAll(IteratedLogicalExpr(e)));
  case OPNOT:
    return AMPL_DISPATCH(VisitNot(NotExpr(e)));
  case OP_IFF:
    return AMPL_DISPATCH(VisitIff(BinaryLogicalExpr(e)));
  case OPIMPELSE:
    return AMPL_DISPATCH(VisitImplication(ImplicationExpr(e)));
  case OPALLDIFF:
    return AMPL_DISPATCH(VisitAllDiff(AllDiffExpr(e)));
  default:
    throw InvalidLogicalExprError(e.opname());
  }
}

#undef AMPL_DISPATCH
}

#endif  // SOLVERS_UTIL_EXPR_H_
