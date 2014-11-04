/*
 Expression classes

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

#ifndef MP_EXPR_H_
#define MP_EXPR_H_

#include <memory>
#include <vector>

#include "mp/format.h"
#include "mp/problem-base.h"

namespace mp {

namespace internal {

struct ExprImpl {
  expr::Kind kind;
};

struct NumericConstantImpl : ExprImpl {
  double value;
};

struct VariableImpl : ExprImpl {
  int index;
};

struct UnaryExprImpl : ExprImpl {
  const ExprImpl *arg;
};

struct BinaryExprImpl : ExprImpl {
  const ExprImpl *lhs;
  const ExprImpl *rhs;
};

struct IfExprImpl : ExprImpl {
  const ExprImpl *condition;
  const ExprImpl *true_expr;
  const ExprImpl *false_expr;
};

struct PLTermImpl : ExprImpl {
  int num_breakpoints;
  int var_index;
  double data[1];
};

struct LogicalConstantImpl : ExprImpl {
  bool value;
};

// Returns true if the non-null expression e is of type ExprClass.
template <typename ExprClass>
bool Is(expr::Kind k);

}  // namespace internal

// Specialize Is<Expr> for the class ExprClass corresponding to a single
// expression kind.
#define MP_SPECIALIZE_IS(ExprClass, expr_kind) \
namespace internal { \
template <> \
inline bool Is<ExprClass>(expr::Kind k) { return k == expr::expr_kind; } \
}

// Specialize Is<Expr> for the class ExprClass corresponding to a range
// of expression kinds [start, end].
#define MP_SPECIALIZE_IS_RANGE(ExprClass, expr_kind) \
namespace internal { \
template <> \
inline bool Is<ExprClass>(expr::Kind k) { \
  return k >= expr::FIRST_##expr_kind && k <= expr::LAST_##expr_kind; \
} \
}

// An expression.
// An Expr object represents a reference to an expression so
// it is cheap to construct and pass by value. A type safe way to
// process expressions of different types is by using ExprVisitor.
class Expr {
 protected:
  const internal::ExprImpl *impl_;

  // A member function representing the true value of SafeBool.
  void True() const {}

  // Safe bool type.
  typedef void (Expr::*SafeBool)() const;

  // Creates an expression from an implementation.
  template <typename TargetExpr>
  static TargetExpr Create(const internal::ExprImpl *impl) {
    assert((!impl || internal::Is<TargetExpr>(impl->kind)) &&
           "invalid expression kind");
    TargetExpr expr;
    expr.impl_ = impl;
    return expr;
  }

  friend class ExprFactory;

 public:
  // Constructs an Expr object representing a null reference to an AMPL
  // expression. The only operation permitted for such expression is
  // copying, assignment and check whether it is null using operator SafeBool.
  Expr() : impl_(0) {}

  // Returns the expression kind.
  expr::Kind kind() const { return impl_->kind; }

  // Returns a value convertible to bool that can be used in conditions but not
  // in comparisons and evaluates to "true" if this expression is not null
  // and "false" otherwise.
  // Example:
  //   void foo(Expr e) {
  //     if (e) {
  //       // Do something if e is not null.
  //     }
  //   }
  operator SafeBool() const { return impl_ != 0 ? &Expr::True : 0; }
};

#define MP_EXPR(Impl) \
 private: \
  const internal::Impl *impl() const { \
    return static_cast<const internal::Impl*>(this->impl_); \
  }

MP_SPECIALIZE_IS_RANGE(Expr, EXPR)

// A numeric expression.
class NumericExpr : public Expr {};
MP_SPECIALIZE_IS_RANGE(NumericExpr, NUMERIC)

// A logical expression.
class LogicalExpr : public Expr {};
MP_SPECIALIZE_IS_RANGE(LogicalExpr, LOGICAL)

// A numeric constant.
// Examples: 42, -1.23e-4
class NumericConstant : public NumericExpr {
  MP_EXPR(NumericConstantImpl)
 public:
  // Returns the value of this constant.
  double value() const { return impl()->value; }
};

MP_SPECIALIZE_IS(NumericConstant, CONSTANT)

// A reference to a variable.
// Example: x
class Variable : public NumericExpr {
  MP_EXPR(VariableImpl)
 public:
  // Returns the index of the referenced variable.
  int index() const { return impl()->index; }
};

MP_SPECIALIZE_IS(Variable, VARIABLE)

// A unary expression.
// Base: base expression class.
template <typename Base>
class BasicUnaryExpr : public Base {
  MP_EXPR(UnaryExprImpl)
 public:
  // Returns the argument of this expression.
  Base arg() const { return Expr::Create<Base>(impl()->arg); }
};

// A unary numeric expression.
// Examples: -x, sin(x), where x is a variable.
typedef BasicUnaryExpr<NumericExpr> UnaryExpr;
MP_SPECIALIZE_IS_RANGE(UnaryExpr, UNARY)

// A binary expression.
// Base: base expression class.
// Arg: argument expression class.
template <typename Base, typename Arg = Base>
class BasicBinaryExpr : public Base {
  MP_EXPR(BinaryExprImpl)
 public:
  // Returns the left-hand side (the first argument) of this expression.
  Arg lhs() const { return Expr::Create<Arg>(impl()->lhs); }

  // Returns the right-hand side (the second argument) of this expression.
  Arg rhs() const { return Expr::Create<Arg>(impl()->rhs); }
};

// A binary numeric expression.
// Examples: x / y, atan2(x, y), where x and y are variables.
typedef BasicBinaryExpr<NumericExpr> BinaryExpr;
MP_SPECIALIZE_IS_RANGE(BinaryExpr, BINARY)

template <typename Base>
class BasicIfExpr : public Base {
  MP_EXPR(IfExprImpl)
 public:
  LogicalExpr condition() const {
    return Expr::Create<LogicalExpr>(impl()->condition);
  }

  Base true_expr() const { return Expr::Create<Base>(impl()->true_expr); }
  Base false_expr() const { return Expr::Create<Base>(impl()->false_expr); }
};

// An if-then-else expression.
// Example: if x != 0 then y else z, where x, y and z are variables.
typedef BasicIfExpr<NumericExpr> IfExpr;
MP_SPECIALIZE_IS(IfExpr, IF)

// A piecewise-linear term.
// Example: <<0; -1, 1>> x, where x is a variable.
class PLTerm : public NumericExpr {
  MP_EXPR(PLTermImpl)
 public:
  // Returns the number of breakpoints in this term.
  int num_breakpoints() const { return impl()->num_breakpoints; }

  // Returns the number of slopes in this term.
  int num_slopes() const { return num_breakpoints() + 1; }

  // Returns a breakpoint with the specified index.
  double breakpoint(int index) const {
    assert(index >= 0 && index < num_breakpoints());
    return impl()->data[2 * index + 1];
  }

  // Returns a slope with the specified index.
  double slope(int index) const {
    assert(index >= 0 && index < num_slopes());
    return impl()->data[2 * index];
  }

  int var_index() const { return impl()->var_index; }
};

MP_SPECIALIZE_IS(PLTerm, PLTERM)

// TODO: numeric expressions

// A logical constant.
// Examples: 0, 1
class LogicalConstant : public LogicalExpr {
  MP_EXPR(LogicalConstantImpl)
 public:
  // Returns the value of this constant.
  bool value() const { return impl()->value; }
};

MP_SPECIALIZE_IS(LogicalConstant, CONSTANT)

class ExprFactory {
 private:
  std::vector<internal::ExprImpl*> exprs_;

  FMT_DISALLOW_COPY_AND_ASSIGN(ExprFactory);

  // Allocates memory for an object of type Impl which is a subclass
  // of ExprImpl.
  // extra_bytes: extra bytes to allocate at the end.
  template <typename Impl>
  Impl *Allocate(expr::Kind kind, int extra_bytes = 0) {
    // Call push_back first to make sure that the impl pointer doesn't leak
    // if push_back throws an exception.
    exprs_.push_back(0);
    Impl *impl = reinterpret_cast<Impl*>(
          ::operator new(sizeof(Impl) + extra_bytes));
    impl->kind = kind;
    exprs_.back() = impl;
    return impl;
  }

 public:
  ExprFactory() {}
  ~ExprFactory() {
    // TODO: delete expressions
  }

  // Makes a numeric constant.
  NumericConstant MakeNumericConstant(double value) {
    internal::NumericConstantImpl *impl =
        Allocate<internal::NumericConstantImpl>(expr::CONSTANT);
    impl->value = value;
    return Expr::Create<NumericConstant>(impl);
  }

  // Makes a variable reference.
  Variable MakeVariable(int index) {
    internal::VariableImpl *impl =
        Allocate<internal::VariableImpl>(expr::VARIABLE);
    impl->index = index;
    return Expr::Create<Variable>(impl);
  }

  // Makes a unary expression.
  UnaryExpr MakeUnary(expr::Kind kind, NumericExpr arg) {
    internal::UnaryExprImpl *impl = Allocate<internal::UnaryExprImpl>(kind);
    impl->arg = arg.impl_;
    return Expr::Create<UnaryExpr>(impl);
  }

  // Makes a binary expression.
  BinaryExpr MakeBinary(expr::Kind kind, NumericExpr lhs, NumericExpr rhs) {
    internal::BinaryExprImpl *impl = Allocate<internal::BinaryExprImpl>(kind);
    impl->lhs = lhs.impl_;
    impl->rhs = rhs.impl_;
    return Expr::Create<BinaryExpr>(impl);
  }

  // Makes an if expression.
  IfExpr MakeIf(LogicalExpr condition,
                NumericExpr true_expr, NumericExpr false_expr) {
    internal::IfExprImpl *impl = Allocate<internal::IfExprImpl>(expr::IF);
    impl->condition = condition.impl_;
    impl->true_expr = true_expr.impl_;
    impl->false_expr = false_expr.impl_;
    return Expr::Create<IfExpr>(impl);
  }

  // A piecewise-linear term builder.
  class PLTermBuilder {
   private:
    internal::PLTermImpl *impl_;
    int slope_index_;
    int breakpoint_index_;

    friend class ExprFactory;

   public:
    explicit PLTermBuilder(internal::PLTermImpl *impl)
      : impl_(impl), slope_index_(0), breakpoint_index_(0) {}

    void AddSlope(double slope) {
      assert(slope_index_ < impl_->num_breakpoints + 1 && "too many slopes");
      impl_->data[2 * slope_index_] = slope;
      ++slope_index_;
    }

    void AddBreakpoint(double breakpoint) {
      assert(breakpoint_index_ < impl_->num_breakpoints &&
             "too many breakpoints");
      impl_->data[2 * breakpoint_index_ + 1] = breakpoint;
      ++breakpoint_index_;
    }
  };

  // Begins building a piecewise-linear term.
  PLTermBuilder BeginPLTerm(int num_breakpoints) {
    assert(num_breakpoints > 0 && "invalid number of breakpoints");
    internal::PLTermImpl *impl = Allocate<internal::PLTermImpl>(
          expr::PLTERM, sizeof(double) * num_breakpoints * 2);
    impl->num_breakpoints = num_breakpoints;
    return PLTermBuilder(impl);
  }

  // Ends building a piecewise-linear term.
  PLTerm EndPLTerm(PLTermBuilder builder, Variable var) {
    internal::PLTermImpl *impl = builder.impl_;
    impl->var_index = var.index();
    // Make sure that all slopes and breakpoints provided.
    assert(builder.slope_index_ == impl->num_breakpoints + 1 &&
           "too few slopes");
    assert(builder.breakpoint_index_ == impl->num_breakpoints &&
           "too few breakpoints");
    return Expr::Create<PLTerm>(impl);
  }

  // TODO: more numeric expressions

  // Makes a logical constant.
  LogicalConstant MakeLogicalConstant(bool value) {
    internal::LogicalConstantImpl *impl =
        Allocate<internal::LogicalConstantImpl>(expr::CONSTANT);
    impl->value = value;
    return Expr::Create<LogicalConstant>(impl);
  }
};
}  // namespace mp

#endif  // MP_EXPR_H_
