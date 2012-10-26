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
#include <string>

extern "C" {
#include "solvers/nlp.h"
#include "solvers/opcode.hd"
}

namespace ampl {

class Expr;

// Anything in namespace internal is AMPL's INTERNAL IMPLEMENTATION DETAIL
// and MUST NOT BE USED DIRECTLY in user code.
namespace internal {

// An expression proxy used for implementing operator-> in iterators.
template <typename ExprT>
class ExprProxy {
 private:
  ExprT expr_;

 public:
  explicit ExprProxy(expr *e) : expr_(e) {}

  const ExprT *operator->() const { return &expr_; }
};

// An expression array iterator.
template <typename ExprT>
class ExprArrayIterator :
  public std::iterator<std::forward_iterator_tag, ExprT> {
 private:
  expr *const *ptr_;

 public:
  explicit ExprArrayIterator(expr *const *p = 0) : ptr_(p) {}

  ExprT operator*() const { return ExprT(*ptr_); }

  internal::ExprProxy<ExprT> operator->() const {
    return internal::ExprProxy<ExprT>(*ptr_);
  }

  ExprArrayIterator &operator++() {
    ++ptr_;
    return *this;
  }

  ExprArrayIterator operator++(int) {
    ExprArrayIterator it(*this);
    ++ptr_;
    return it;
  }

  bool operator==(ExprArrayIterator other) const { return ptr_ == other.ptr_; }
  bool operator!=(ExprArrayIterator other) const { return ptr_ != other.ptr_; }
};
}

// An operation type.
// Numeric values for the operation types should be in sync with the ones in
// op_type.hd.
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
 private:
  static const char *const OP_NAMES[N_OPS];

  bool IsOpCodeInRange() const {
    return opcode() >= 0 && opcode() < N_OPS;
  }

  void True() const {}
  typedef void (Expr::*SafeBool)() const;

  friend class ExprBuilder;

 protected:
  expr *expr_;

  // Constructs an Expr object representing a reference to an AMPL
  // expression e. Only a minimal check is performed when assertions are
  // enabled to make sure that the opcode is within the valid range.
  explicit Expr(expr *e) : expr_(e) {
    assert(!expr_ || IsOpCodeInRange());
  }

  // Returns true iff this expression is null or has type t.
  bool HasTypeOrNull(OpType t) const {
    return !expr_ || optype() == t;
  }

 public:
  // Constructs an Expr object representing a null reference to an AMPL
  // expression. The only operation permitted for such expression is
  // copying, assignment and check whether it is null using operator SafeBool.
  Expr() : expr_() {}

  // Returns a value convertible to bool that can be used in conditions but not
  // in comparisons and evaluates to "true" if this expression is not null
  // and "false" otherwise.
  // Example:
  //   void foo(Expr e) {
  //     if (e) {
  //       // Do something if e is not null.
  //     }
  //   }
  operator SafeBool() const { return expr_ ? &Expr::True : 0; }

  // Returns the operation code (opcode) of this expression which should be
  // non-null. The opcodes are defined in opcode.hd.
  int opcode() const {
    return reinterpret_cast<std::size_t>(expr_->op);
  }

  // Returns the operation name of this expression which should be non-null.
  const char *opname() const {
    assert(IsOpCodeInRange());
    return OP_NAMES[opcode()];
  }

  // Returns the operation type of this expression which can be unary, binary,
  // etc. It is called "optype" rather than simply "type" to avoid confusion
  // with expression types such as logical or numeric. This expression should
  // be non-null.
  OpType optype() const {
    assert(IsOpCodeInRange());
    return static_cast<OpType>(::optype[opcode()]);
  }

  bool operator==(Expr other) const { return expr_ == other.expr_; }
  bool operator!=(Expr other) const { return expr_ != other.expr_; }

  // Recursively compares two expressions and returns true if they are equal.
  friend bool AreEqual(Expr e1, Expr e2);
};

// Casts an expression to type T. Returns a null expression if the cast
// is not possible.
template <typename T>
T Cast(Expr e);

// A numeric expression.
class NumericExpr : public Expr {
 private:
  // Returns true if this is a valid numeric expressions.
  bool IsValid() const;

 protected:
  NumericExpr(Expr e) : Expr(e) {}

 public:
  // Constructs a NumericExpr object.
  explicit NumericExpr(expr *e = 0) : Expr(e) {
    assert(IsValid());
  }
};

// A logical or constraint expression.
class LogicalExpr : public Expr {
 private:
  // Returns true if this is a valid logical expressions.
  bool IsValid() const;

 public:
  // Constructs a LogicalExpr object.
  explicit LogicalExpr(expr *e = 0) : Expr(e) {
    assert(IsValid());
  }
};

// A unary numeric expression.
// Examples: -x, sin(x), where x is a variable.
class UnaryExpr : public NumericExpr {
 private:
  explicit UnaryExpr(NumericExpr e) : NumericExpr(e) {
    assert(HasTypeOrNull(OPTYPE_UNARY));
  }

  template <typename Impl, typename Result, typename LResult>
  friend class ExprVisitor;
  friend class ExprBuilder;

 public:
  UnaryExpr() {}

  // Returns the argument of this expression.
  NumericExpr arg() const { return NumericExpr(expr_->L.e); }
};

// A binary numeric expression.
// Examples: x / y, atan2(x, y), where x and y are variables.
class BinaryExpr : public NumericExpr {
 private:
  explicit BinaryExpr(NumericExpr e) : NumericExpr(e) {
    assert(!expr_ || opcode() == OP1POW || opcode() == OPCPOW ||
           HasTypeOrNull(OPTYPE_BINARY));
  }

  template <typename Impl, typename Result, typename LResult>
  friend class ExprVisitor;
  friend class ExprBuilder;

 public:
  BinaryExpr() {}

  // Returns the left-hand side (the first argument) of this expression.
  NumericExpr lhs() const { return NumericExpr(expr_->L.e); }

  // Returns the right-hand side (the second argument) of this expression.
  NumericExpr rhs() const { return NumericExpr(expr_->R.e); }
};

// A numeric expression with a variable number of arguments.
// Example: min{i in I} x[i], where I is a set and x is a variable.
class VarArgExpr : public NumericExpr {
 private:
  explicit VarArgExpr(NumericExpr e) : NumericExpr(e) {
    assert(HasTypeOrNull(OPTYPE_VARARG));
  }

  template <typename Impl, typename Result, typename LResult>
  friend class ExprVisitor;
  friend class ExprBuilder;

  static const de END;

 public:
  VarArgExpr() {}

  // An argument iterator.
  class iterator :
    public std::iterator<std::forward_iterator_tag, NumericExpr> {
   private:
    const de *de_;

    friend class VarArgExpr;

    explicit iterator(const de *d) : de_(d) {}

   public:
    iterator() : de_(&END) {}

    NumericExpr operator*() const { return NumericExpr(de_->e); }

    internal::ExprProxy<NumericExpr> operator->() const {
      return internal::ExprProxy<NumericExpr>(de_->e);
    }

    iterator &operator++() {
      ++de_;
      return *this;
    }

    iterator operator++(int) {
      iterator it(*this);
      ++de_;
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

// A sum expression.
// Example: sum{i in I} x[i], where I is a set and x is a variable.
class SumExpr : public NumericExpr {
 private:
  explicit SumExpr(NumericExpr e) : NumericExpr(e) {
    assert(HasTypeOrNull(OPTYPE_SUM));
  }

  template <typename Impl, typename Result, typename LResult>
  friend class ExprVisitor;

 public:
  SumExpr() {}

  typedef internal::ExprArrayIterator<NumericExpr> iterator;

  iterator begin() const {
    return iterator(expr_->L.ep);
  }

  iterator end() const {
    return iterator(expr_->R.ep);
  }
};

// A count expression.
// Example: count{i in I} (x[i] >= 0), where I is a set and x is a variable.
class CountExpr : public NumericExpr {
 private:
  explicit CountExpr(NumericExpr e) : NumericExpr(e) {
    assert(HasTypeOrNull(OPTYPE_COUNT));
  }

  template <typename Impl, typename Result, typename LResult>
  friend class ExprVisitor;

 public:
  CountExpr() {}

  typedef internal::ExprArrayIterator<LogicalExpr> iterator;

  iterator begin() const {
    return iterator(expr_->L.ep);
  }

  iterator end() const {
    return iterator(expr_->R.ep);
  }
};

// An if-then-else expression.
// Example: if x != 0 then y else z, where x, y and z are variables.
class IfExpr : public NumericExpr {
 private:
  explicit IfExpr(NumericExpr e) : NumericExpr(e) {
    assert(HasTypeOrNull(OPTYPE_IF));
  }

  template <typename Impl, typename Result, typename LResult>
  friend class ExprVisitor;

 public:
  IfExpr() {}

  LogicalExpr condition() const {
    return LogicalExpr(reinterpret_cast<expr_if*>(expr_)->e);
  }

  NumericExpr true_expr() const {
    return NumericExpr(reinterpret_cast<expr_if*>(expr_)->T);
  }

  NumericExpr false_expr() const {
    return NumericExpr(reinterpret_cast<expr_if*>(expr_)->F);
  }
};

// A piecewise-linear term.
// Example: <<0; -1, 1>> x, where x is a variable.
class PiecewiseLinearTerm : public NumericExpr {
 private:
  explicit PiecewiseLinearTerm(NumericExpr e) : NumericExpr(e) {
    assert(HasTypeOrNull(OPTYPE_PLTERM));
  }

  template <typename Impl, typename Result, typename LResult>
  friend class ExprVisitor;

 public:
  PiecewiseLinearTerm() {}

  // Returns the number of slopes in this term.
  int num_slopes() const {
    assert(expr_->L.p->n >= 1);
    return expr_->L.p->n;
  }

  // Returns the number of breakpoints in this term.
  int num_breakpoints() const {
    return num_slopes() - 1;
  }

  // Returns the number of slopes in this term.
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
// Examples: 42, -1.23e-4
class NumericConstant : public NumericExpr {
 private:
  explicit NumericConstant(Expr e) : NumericExpr(e) {
    assert(HasTypeOrNull(OPTYPE_NUMBER));
  }

  template <typename Impl, typename Result, typename LResult>
  friend class ExprVisitor;

  friend NumericConstant Cast<NumericConstant>(Expr);

 public:
  NumericConstant() {}

  // Returns the value of this number.
  double value() const { return reinterpret_cast<expr_n*>(expr_)->v; }
};

template <>
inline NumericConstant Cast<NumericConstant>(Expr e) {
  return NumericConstant(e.opcode() == OPNUM ? e : Expr());
}

// A reference to a variable.
// Example: x
class Variable : public NumericExpr {
 private:
  explicit Variable(Expr e) : NumericExpr(e) {
    assert(HasTypeOrNull(OPTYPE_VARIABLE));
  }

  template <typename Impl, typename Result, typename LResult>
  friend class ExprVisitor;

  friend Variable Cast<Variable>(Expr);

 public:
  Variable() {}

  // Returns the index of the referenced variable.
  int index() const { return expr_->a; }
};

template <>
inline Variable Cast<Variable>(Expr e) {
  return Variable(e.opcode() == OPVARVAL ? e : Expr());
}

// A numberof expression.
// Example: numberof 42 in ({i in I} x[i]),
// where I is a set and x is a variable.
class NumberOfExpr : public NumericExpr {
 private:
  explicit NumberOfExpr(NumericExpr e) : NumericExpr(e) {
    assert(!expr_ || opcode() == OPNUMBEROF);
  }

  template <typename Impl, typename Result, typename LResult>
  friend class ExprVisitor;

 public:
  NumberOfExpr() {}

  NumericExpr target() const { return NumericExpr(*expr_->L.ep); }

  typedef internal::ExprArrayIterator<NumericExpr> iterator;

  iterator begin() const {
    return iterator(expr_->L.ep + 1);
  }

  iterator end() const {
    return iterator(expr_->R.ep);
  }
};

// A logical constant.
// Examples: 0, 1
class LogicalConstant : public LogicalExpr {
 private:
  explicit LogicalConstant(LogicalExpr e) : LogicalExpr(e) {
    assert(HasTypeOrNull(OPTYPE_NUMBER));
  }

  template <typename Impl, typename Result, typename LResult>
  friend class ExprVisitor;

 public:
  LogicalConstant() {}

  // Returns the value of this constant.
  bool value() const { return reinterpret_cast<expr_n*>(expr_)->v != 0; }
};

// A relational expression.
// Examples: x < y, x != y, where x and y are variables.
class RelationalExpr : public LogicalExpr {
 private:
  explicit RelationalExpr(LogicalExpr e) : LogicalExpr(e) {
    assert(HasTypeOrNull(OPTYPE_BINARY));
  }

  template <typename Impl, typename Result, typename LResult>
  friend class ExprVisitor;

 public:
  RelationalExpr() {}

  // Returns the left-hand side (the first argument) of this expression.
  NumericExpr lhs() const { return NumericExpr(expr_->L.e); }

  // Returns the right-hand side (the second argument) of this expression.
  NumericExpr rhs() const { return NumericExpr(expr_->R.e); }
};

// A logical NOT expression.
// Example: not a, where a is a logical expression.
class NotExpr : public LogicalExpr {
 private:
  explicit NotExpr(LogicalExpr e) : LogicalExpr(e) {
    assert(!expr_ || opcode() == OPNOT);
  }

  template <typename Impl, typename Result, typename LResult>
  friend class ExprVisitor;

 public:
  NotExpr() {}

  // Returns the argument of this expression.
  LogicalExpr arg() const { return LogicalExpr(expr_->L.e); }
};

// A binary logical expression.
// Examples: a || b, a && b, where a and b are logical expressions.
class BinaryLogicalExpr : public LogicalExpr {
 private:
  explicit BinaryLogicalExpr(LogicalExpr e) : LogicalExpr(e) {
    assert(HasTypeOrNull(OPTYPE_BINARY));
  }

  template <typename Impl, typename Result, typename LResult>
  friend class ExprVisitor;

 public:
  BinaryLogicalExpr() {}

  // Returns the left-hand side (the first argument) of this expression.
  LogicalExpr lhs() const { return LogicalExpr(expr_->L.e); }

  // Returns the right-hand side (the second argument) of this expression.
  LogicalExpr rhs() const { return LogicalExpr(expr_->R.e); }
};

// An implication expression.
// Example: a ==> b else c, where a, b and c are logical expressions.
class ImplicationExpr : public LogicalExpr {
 private:
  explicit ImplicationExpr(LogicalExpr e) : LogicalExpr(e) {
    assert(!expr_ || opcode() == OPIMPELSE);
  }

  template <typename Impl, typename Result, typename LResult>
  friend class ExprVisitor;

 public:
  ImplicationExpr() {}

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

// An iterated logical expression.
// Example: exists{i in I} x[i] >= 0, where I is a set and x is a variable.
class IteratedLogicalExpr : public LogicalExpr {
 private:
  explicit IteratedLogicalExpr(LogicalExpr e) : LogicalExpr(e) {
    assert(HasTypeOrNull(OPTYPE_SUM) || HasTypeOrNull(OPTYPE_COUNT));
  }

  template <typename Impl, typename Result, typename LResult>
  friend class ExprVisitor;

 public:
  IteratedLogicalExpr() {}

  typedef internal::ExprArrayIterator<LogicalExpr> iterator;

  iterator begin() const {
    return iterator(expr_->L.ep);
  }

  iterator end() const {
    return iterator(expr_->R.ep);
  }
};

// An alldiff expression.
// Example: alldiff{i in I} x[i], where I is a set and x is a variable.
class AllDiffExpr : public LogicalExpr {
 private:
  explicit AllDiffExpr(LogicalExpr e) : LogicalExpr(e) {
    assert(HasTypeOrNull(OPTYPE_COUNT));
  }

  template <typename Impl, typename Result, typename LResult>
  friend class ExprVisitor;

 public:
  AllDiffExpr() {}

  typedef internal::ExprArrayIterator<NumericExpr> iterator;

  iterator begin() const {
    return iterator(expr_->L.ep);
  }

  iterator end() const {
    return iterator(expr_->R.ep);
  }
};

// A general error.
class Error : public std::runtime_error {
public:
  explicit Error(const std::string &message) : std::runtime_error(message) {}
};

// An exception that is thrown when an ASL expression not supported
// by the solver is encountered.
class UnsupportedExprError : public Error {
public:
  explicit UnsupportedExprError(const char *expr) :
    Error(std::string("unsupported expression: ") + expr) {}
};

// An exception that is thrown when an invalid numeric expression
// is encountered.
class InvalidNumericExprError : public Error {
public:
  explicit InvalidNumericExprError(NumericExpr e) :
    Error(std::string("invalid numeric expression: ") + e.opname()) {}
};

// An exception that is thrown when an invalid logical or constraint
// expression is encountered.
class InvalidLogicalExprError : public Error {
public:
  explicit InvalidLogicalExprError(LogicalExpr e) :
    Error(std::string("invalid logical expression: ") + e.opname()) {}
};

#define AMPL_DISPATCH(call) static_cast<Impl*>(this)->call

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
//    double VisitConstant(NumericConstant n) { return n.value(); }
//  };
//
// ExprVisitor uses the curiously recurring template pattern:
// http://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
template <typename Impl, typename Result, typename LResult>
class ExprVisitor {
 public:
  Result Visit(NumericExpr e);
  LResult Visit(LogicalExpr e);

  Result VisitInvalidNumericExpr(NumericExpr e) {
    throw InvalidNumericExprError(e);
  }

  LResult VisitInvalidLogicalExpr(LogicalExpr e) {
    throw InvalidLogicalExprError(e);
  }

  Result VisitUnhandledNumericExpr(NumericExpr e) {
    throw UnsupportedExprError(e.opname());
  }

  LResult VisitUnhandledLogicalExpr(LogicalExpr e) {
    throw UnsupportedExprError(e.opname());
  }

  Result VisitPlus(BinaryExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitMinus(BinaryExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitMult(BinaryExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitDiv(BinaryExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitRem(BinaryExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitPow(BinaryExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitNumericLess(BinaryExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitMin(VarArgExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitMax(VarArgExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitFloor(UnaryExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitCeil(UnaryExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitAbs(UnaryExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitUnaryMinus(UnaryExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitIf(IfExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitTanh(UnaryExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitTan(UnaryExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitSqrt(UnaryExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitSinh(UnaryExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitSin(UnaryExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitLog10(UnaryExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitLog(UnaryExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitExp(UnaryExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitCosh(UnaryExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitCos(UnaryExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitAtanh(UnaryExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitAtan2(BinaryExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitAtan(UnaryExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitAsinh(UnaryExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitAsin(UnaryExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitAcosh(UnaryExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitAcos(UnaryExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitSum(SumExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitIntDiv(BinaryExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitPrecision(BinaryExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitRound(BinaryExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitTrunc(BinaryExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitCount(CountExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitNumberOf(NumberOfExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitPLTerm(PiecewiseLinearTerm t) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(t));
  }

  Result VisitConstExpPow(BinaryExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitPow2(UnaryExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitConstBasePow(BinaryExpr e) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitNumericConstant(NumericConstant c) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(c));
  }

  Result VisitVariable(Variable v) {
    return AMPL_DISPATCH(VisitUnhandledNumericExpr(v));
  }

  LResult VisitOr(BinaryLogicalExpr e) {
    return AMPL_DISPATCH(VisitUnhandledLogicalExpr(e));
  }

  LResult VisitAnd(BinaryLogicalExpr e) {
    return AMPL_DISPATCH(VisitUnhandledLogicalExpr(e));
  }

  LResult VisitLess(RelationalExpr e) {
    return AMPL_DISPATCH(VisitUnhandledLogicalExpr(e));
  }

  LResult VisitLessEqual(RelationalExpr e) {
    return AMPL_DISPATCH(VisitUnhandledLogicalExpr(e));
  }

  LResult VisitEqual(RelationalExpr e) {
    return AMPL_DISPATCH(VisitUnhandledLogicalExpr(e));
  }

  LResult VisitGreaterEqual(RelationalExpr e) {
    return AMPL_DISPATCH(VisitUnhandledLogicalExpr(e));
  }

  LResult VisitGreater(RelationalExpr e) {
    return AMPL_DISPATCH(VisitUnhandledLogicalExpr(e));
  }

  LResult VisitNotEqual(RelationalExpr e) {
    return AMPL_DISPATCH(VisitUnhandledLogicalExpr(e));
  }

  LResult VisitNot(NotExpr e) {
    return AMPL_DISPATCH(VisitUnhandledLogicalExpr(e));
  }

  LResult VisitAtLeast(RelationalExpr e) {
    return AMPL_DISPATCH(VisitUnhandledLogicalExpr(e));
  }

  LResult VisitAtMost(RelationalExpr e) {
    return AMPL_DISPATCH(VisitUnhandledLogicalExpr(e));
  }

  LResult VisitExactly(RelationalExpr e) {
    return AMPL_DISPATCH(VisitUnhandledLogicalExpr(e));
  }

  LResult VisitNotAtLeast(RelationalExpr e) {
    return AMPL_DISPATCH(VisitUnhandledLogicalExpr(e));
  }

  LResult VisitNotAtMost(RelationalExpr e) {
    return AMPL_DISPATCH(VisitUnhandledLogicalExpr(e));
  }

  LResult VisitNotExactly(RelationalExpr e) {
    return AMPL_DISPATCH(VisitUnhandledLogicalExpr(e));
  }

  LResult VisitForAll(IteratedLogicalExpr e) {
    return AMPL_DISPATCH(VisitUnhandledLogicalExpr(e));
  }

  LResult VisitExists(IteratedLogicalExpr e) {
    return AMPL_DISPATCH(VisitUnhandledLogicalExpr(e));
  }

  LResult VisitImplication(ImplicationExpr e) {
    return AMPL_DISPATCH(VisitUnhandledLogicalExpr(e));
  }

  LResult VisitIff(BinaryLogicalExpr e) {
    return AMPL_DISPATCH(VisitUnhandledLogicalExpr(e));
  }

  LResult VisitAllDiff(AllDiffExpr e) {
    return AMPL_DISPATCH(VisitUnhandledLogicalExpr(e));
  }

  LResult VisitLogicalConstant(LogicalConstant c) {
    return AMPL_DISPATCH(VisitUnhandledLogicalExpr(c));
  }
};

template <typename Impl, typename Result, typename LResult>
Result ExprVisitor<Impl, Result, LResult>::Visit(NumericExpr e) {
  // All expressions except OPNUMBEROFs, OPIFSYM, OPFUNCALL and OPHOL
  // are supported.
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
    return AMPL_DISPATCH(VisitNumericLess(BinaryExpr(e)));
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
    return AMPL_DISPATCH(VisitUnaryMinus(UnaryExpr(e)));
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
  case OPCOUNT:
    return AMPL_DISPATCH(VisitCount(CountExpr(e)));
  case OPNUMBEROF:
    return AMPL_DISPATCH(VisitNumberOf(NumberOfExpr(e)));
  case OPPLTERM:
    return AMPL_DISPATCH(VisitPLTerm(PiecewiseLinearTerm(e)));
  case OP1POW:
    return AMPL_DISPATCH(VisitConstExpPow(BinaryExpr(e)));
  case OP2POW:
    return AMPL_DISPATCH(VisitPow2(UnaryExpr(e)));
  case OPCPOW:
    return AMPL_DISPATCH(VisitConstBasePow(BinaryExpr(e)));
  case OPNUM:
    return AMPL_DISPATCH(VisitNumericConstant(NumericConstant(e)));
  case OPVARVAL:
    return AMPL_DISPATCH(VisitVariable(Variable(e)));
  default:
    return AMPL_DISPATCH(VisitInvalidNumericExpr(e));
  }
}

template <typename Impl, typename Result, typename LResult>
LResult ExprVisitor<Impl, Result, LResult>::Visit(LogicalExpr e) {
  switch (e.opcode()) {
  case OPOR:
    return AMPL_DISPATCH(VisitOr(BinaryLogicalExpr(e)));
  case OPAND:
    return AMPL_DISPATCH(VisitAnd(BinaryLogicalExpr(e)));
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
  case OPNOT:
    return AMPL_DISPATCH(VisitNot(NotExpr(e)));
  case OPATLEAST:
    return AMPL_DISPATCH(VisitAtLeast(RelationalExpr(e)));
  case OPATMOST:
    return AMPL_DISPATCH(VisitAtMost(RelationalExpr(e)));
  case OPEXACTLY:
    return AMPL_DISPATCH(VisitExactly(RelationalExpr(e)));
  case OPNOTATLEAST:
    return AMPL_DISPATCH(VisitNotAtLeast(RelationalExpr(e)));
  case OPNOTATMOST:
    return AMPL_DISPATCH(VisitNotAtMost(RelationalExpr(e)));
  case OPNOTEXACTLY:
    return AMPL_DISPATCH(VisitNotExactly(RelationalExpr(e)));
  case ANDLIST:
    return AMPL_DISPATCH(VisitForAll(IteratedLogicalExpr(e)));
  case ORLIST:
    return AMPL_DISPATCH(VisitExists(IteratedLogicalExpr(e)));
  case OPIMPELSE:
    return AMPL_DISPATCH(VisitImplication(ImplicationExpr(e)));
  case OP_IFF:
    return AMPL_DISPATCH(VisitIff(BinaryLogicalExpr(e)));
  case OPALLDIFF:
    return AMPL_DISPATCH(VisitAllDiff(AllDiffExpr(e)));
  case OPNUM:
    return AMPL_DISPATCH(VisitLogicalConstant(LogicalConstant(e)));
  default:
    return AMPL_DISPATCH(VisitInvalidLogicalExpr(e));
  }
}

#undef AMPL_DISPATCH

class ExprChecker : public ExprVisitor<ExprChecker, bool, bool> {
 public:
  bool VisitInvalidNumericExpr(NumericExpr) { return false; }
  bool VisitInvalidLogicalExpr(LogicalExpr) { return false; }

  bool VisitUnhandledNumericExpr(NumericExpr) { return true; }
  bool VisitUnhandledLogicalExpr(LogicalExpr) { return true; }
};

inline bool NumericExpr::IsValid() const {
  return !expr_ || ExprChecker().Visit(*this);
}

inline bool LogicalExpr::IsValid() const {
  return !expr_ || ExprChecker().Visit(*this);
}
}

#endif  // SOLVERS_UTIL_EXPR_H_
