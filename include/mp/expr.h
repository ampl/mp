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

struct UnaryExprImpl : ExprImpl {
  const ExprImpl *arg;
};

struct NumericConstantImpl : ExprImpl {
  double value;
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
    assert(!impl || internal::Is<TargetExpr>(impl->kind));
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

MP_SPECIALIZE_IS_RANGE(Expr, EXPR)

// A numeric expression.
class NumericExpr : public Expr {};

MP_SPECIALIZE_IS_RANGE(NumericExpr, NUMERIC)

// An unary numeric expression.
class UnaryExpr : public NumericExpr {
 private:
  const internal::UnaryExprImpl *impl() const {
    return static_cast<const internal::UnaryExprImpl*>(impl_);
  }

 public:
  // Returns the argument of this expression.
  NumericExpr arg() const { return Create<NumericExpr>(impl()->arg); }
};

MP_SPECIALIZE_IS_RANGE(UnaryExpr, UNARY)

// A numeric constant.
// Examples: 42, -1.23e-4
class NumericConstant : public NumericExpr {
 private:
  const internal::NumericConstantImpl *impl() const {
    return static_cast<const internal::NumericConstantImpl*>(impl_);
  }

 public:
  // Returns the value of this constant.
  double value() const { return impl()->value; }
};

MP_SPECIALIZE_IS(NumericConstant, CONSTANT)

class ExprFactory {
 private:
  std::vector<internal::ExprImpl*> exprs_;

  FMT_DISALLOW_COPY_AND_ASSIGN(ExprFactory);

  // Allocates memory for an object of type Impl which is a subclass
  // of ExprImpl.
  template <typename Impl>
  Impl *Allocate() {
    // Call push_back first to make sure that the impl pointer doesn't leak
    // if push_back throws an exception.
    exprs_.push_back(0);
    Impl *impl = new Impl();
    exprs_.back() = impl;
    return impl;
  }

 public:
  ExprFactory() {}
  ~ExprFactory() {
    // TODO: delete expressions
  }

  // Makes an unary expression.
  UnaryExpr MakeUnary(expr::Kind kind, NumericExpr arg) {
    internal::UnaryExprImpl *impl = Allocate<internal::UnaryExprImpl>();
    impl->kind = kind;
    impl->arg = arg.impl_;
    return Expr::Create<UnaryExpr>(impl);
  }

  // Makes a numeric constant.
  NumericConstant MakeNumericConstant(double value) {
    internal::NumericConstantImpl *impl =
        Allocate<internal::NumericConstantImpl>();
    impl->value = value;
    return Expr::Create<NumericConstant>(impl);
  }
};
}  // namespace mp

#endif  // MP_EXPR_H_
