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

#include <cassert>
#include <memory>
#include <vector>

#include "mp/format.h"
#include "mp/problem-base.h"

#ifndef MP_ASSERT
# define MP_ASSERT(condition, message) assert((condition) && message)
#endif

namespace mp {

class ExprFactory;

namespace internal {
// Returns true if the non-null expression e is of type ExprType.
template <typename ExprType>
bool Is(expr::Kind k);
}

// Specialize Is<Expr> for the class ExprType corresponding to a single
// expression kind.
#define MP_SPECIALIZE_IS(ExprType, expr_kind) \
namespace internal { \
template <> \
inline bool Is<ExprType>(expr::Kind k) { return k == expr::expr_kind; } \
}

// Specialize Is<Expr> for the class ExprType corresponding to a range
// of expression kinds [start, end].
#define MP_SPECIALIZE_IS_RANGE(ExprType, expr_kind) \
namespace internal { \
template <> \
inline bool Is<ExprType>(expr::Kind k) { \
  return k >= expr::FIRST_##expr_kind && k <= expr::LAST_##expr_kind; \
} \
}

// An expression.
// An Expr object represents a reference to an expression so
// it is cheap to construct and pass by value. A type safe way to
// process expressions of different types is by using ExprVisitor.
class Expr {
 protected:
  class Impl {
   private:
    // Only ExprFactory should be able to set Impl::kind_.
    expr::Kind kind_;

    friend class mp::ExprFactory;

   public:
    expr::Kind kind() const { return kind_; }
  };

 private:
  const Impl *impl_;

  // A member function representing the true value of SafeBool.
  void True() const {}

  // Safe bool type.
  typedef void (Expr::*SafeBool)() const;

  // An expression proxy used for implementing operator-> in iterators.
  template <typename ExprType>
  class Proxy {
   private:
    ExprType expr_;

   public:
    explicit Proxy(const Expr::Impl *e) : expr_(Create<ExprType>(e)) {}

    const ExprType *operator->() const { return &expr_; }
  };

 protected:
  // Returns a pointer to the implementation.
  const Impl *impl() const { return impl_; }

  // Creates an expression from an implementation.
  template <typename TargetExpr>
  static TargetExpr Create(const Expr::Impl *impl) {
    MP_ASSERT((!impl || internal::Is<TargetExpr>(impl->kind())),
              "invalid expression kind");
    TargetExpr expr;
    expr.impl_ = impl;
    return expr;
  }

  friend class ExprFactory;

  // An expression iterator.
  template <typename ExprType>
  class BasicIterator :
      public std::iterator<std::forward_iterator_tag, ExprType> {
   private:
    const Expr::Impl *const *ptr_;

   public:
    explicit BasicIterator(const Expr::Impl *const *p = 0) : ptr_(p) {}

    ExprType operator*() const { return Create<ExprType>(*ptr_); }

    Proxy<ExprType> operator->() const { return Proxy<ExprType>(*ptr_); }

    BasicIterator &operator++() {
      ++ptr_;
      return *this;
    }

    BasicIterator operator++(int ) {
      BasicIterator it(*this);
      ++ptr_;
      return it;
    }

    bool operator==(BasicIterator other) const { return ptr_ == other.ptr_; }
    bool operator!=(BasicIterator other) const { return ptr_ != other.ptr_; }
  };

 public:
  // Constructs an Expr object representing a null reference to an
  // expression. The only operation permitted for such object is copying,
  // assignment and check whether it is null using operator SafeBool.
  Expr() : impl_(0) {}

  // Returns the expression kind.
  expr::Kind kind() const { return impl_->kind(); }

  // Returns a value convertible to bool that can be used in conditions but not
  // in comparisons and evaluates to "true" if this expression is not null
  // and "false" otherwise.
  // Example:
  //   if (e) {
  //     // Do something if e is not null.
  //   }
  operator SafeBool() const { return impl_ != 0 ? &Expr::True : 0; }
};

#define MP_EXPR \
  const Impl *impl() const { return static_cast<const Impl*>(Expr::impl()); } \
  friend class ExprFactory

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
 private:
  struct Impl : Expr::Impl {
    double value;
  };
  MP_EXPR;

 public:
  // Returns the value of this constant.
  double value() const { return impl()->value; }
};

MP_SPECIALIZE_IS(NumericConstant, CONSTANT)

// A reference to a variable.
// Example: x
class Variable : public NumericExpr {
 private:
  struct Impl : Expr::Impl {
    int index;
  };
  MP_EXPR;

 public:
  // Returns the index of the referenced variable.
  int index() const { return impl()->index; }
};

MP_SPECIALIZE_IS(Variable, VARIABLE)

// A unary expression.
// Base: base expression class.
template <typename Base>
class BasicUnaryExpr : public Base {
 private:
  struct Impl : Expr::Impl {
    const Expr::Impl *arg;
  };
  MP_EXPR;

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
 private:
  struct Impl : Expr::Impl {
    const Expr::Impl *lhs;
    const Expr::Impl *rhs;
  };
  MP_EXPR;

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

template <typename Base, expr::Kind K>
class BasicIfExpr : public Base {
 private:
  struct Impl : Expr::Impl {
    const Expr::Impl *condition;
    const Expr::Impl *true_expr;
    const Expr::Impl *false_expr;
  };
  MP_EXPR;

 public:
  static const expr::Kind KIND = K;

  LogicalExpr condition() const {
    return Expr::Create<LogicalExpr>(impl()->condition);
  }

  Base true_expr() const { return Expr::Create<Base>(impl()->true_expr); }
  Base false_expr() const { return Expr::Create<Base>(impl()->false_expr); }
};

// An if-then-else expression.
// Example: if x != 0 then y else z, where x, y and z are variables.
typedef BasicIfExpr<NumericExpr, expr::IF> IfExpr;
MP_SPECIALIZE_IS(IfExpr, IF)

// A piecewise-linear term.
// Example: <<0; -1, 1>> x, where x is a variable.
class PLTerm : public NumericExpr {
 private:
  struct Impl : Expr::Impl {
    int num_breakpoints;
    int var_index;
    double data[1];
  };
  MP_EXPR;

 public:
  // Returns the number of breakpoints in this term.
  int num_breakpoints() const { return impl()->num_breakpoints; }

  // Returns the number of slopes in this term.
  int num_slopes() const { return num_breakpoints() + 1; }

  // Returns a breakpoint with the specified index.
  double breakpoint(int index) const {
    MP_ASSERT(index >= 0 && index < num_breakpoints(), "index out of bounds");
    return impl()->data[2 * index + 1];
  }

  // Returns a slope with the specified index.
  double slope(int index) const {
    MP_ASSERT(index >= 0 && index < num_slopes(), "index out of bounds");
    return impl()->data[2 * index];
  }

  int var_index() const { return impl()->var_index; }
};

MP_SPECIALIZE_IS(PLTerm, PLTERM)

// A reference to a function.
class Function {
 private:
  struct Impl {};

  const Impl *impl_;

  // A member function representing the true value of SafeBool.
  void True() const {}

  // Safe bool type.
  typedef void (Function::*SafeBool)() const;

  explicit Function(const Impl *impl) : impl_(impl) {}

  friend class CallExpr;
  friend class ExprFactory;

 public:
  // Constructs a Function object representing a null reference to a
  // function. The only operation permitted for such object is copying,
  // assignment and check whether it is null using operator SafeBool.
  Function() : impl_(0) {}

  // Returns a value convertible to bool that can be used in conditions but not
  // in comparisons and evaluates to "true" if this function is not null
  // and "false" otherwise.
  // Example:
  //   if (f) {
  //     // Do something if f is not null.
  //   }
  operator SafeBool() const { return impl_ != 0 ? &Function::True : 0; }
};

// A function call expression.
// Example: f(x), where f is a function and x is a variable.
class CallExpr : public NumericExpr {
 private:
  struct Impl : Expr::Impl {
    const Function::Impl *func;
    int num_args;
    const Expr::Impl *args[1];
  };
  MP_EXPR;

 public:
  Function function() const { return Function(impl()->func); }

  typedef Expr Arg;

  // Returns the number of arguments.
  int num_args() const { return impl()->num_args; }

  // Returns an argument with the specified index.
  Expr arg(int index) {
    MP_ASSERT(index >= 0 && index < num_args(), "index out of bounds");
    return Create<Expr>(impl()->args[index]);
  }

  // An argument iterator.
  typedef BasicIterator<Expr> iterator;

  iterator begin() const { return iterator(impl()->args); }
  iterator end() const { return iterator(impl()->args + num_args()); }
};

MP_SPECIALIZE_IS(CallExpr, CALL)

template <expr::Kind KIND, typename Base = NumericExpr, typename ArgType = Base>
class BasicIteratedExpr : public Base {
 private:
  struct Impl : Expr::Impl {
    int num_args;
    const Expr::Impl *args[1];
  };
  MP_EXPR;

 public:
  typedef ArgType Arg;

  // Returns the number of arguments.
  int num_args() const { return impl()->num_args; }

  // Returns an argument with the specified index.
  Arg arg(int index) {
    MP_ASSERT(index >= 0 && index < num_args(), "index out of bounds");
    return Expr::Create<Arg>(impl()->args[index]);
  }

  // An argument iterator.
  typedef Expr::BasicIterator<Arg> iterator;

  iterator begin() const { return iterator(impl()->args); }
  iterator end() const { return iterator(impl()->args + num_args()); }
};

// A numeric expression with a variable number of arguments.
// The min and max functions always have at least one argument.
// Example: min{i in I} x[i], where I is a set and x is a variable.
typedef BasicIteratedExpr<expr::FIRST_VARARG> VarArgExpr;
MP_SPECIALIZE_IS_RANGE(VarArgExpr, VARARG)

// A sum expression.
// Example: sum{i in I} x[i], where I is a set and x is a variable.
typedef BasicIteratedExpr<expr::SUM> SumExpr;
MP_SPECIALIZE_IS(SumExpr, SUM)

// A count expression.
// Example: count{i in I} (x[i] >= 0), where I is a set and x is a variable.
typedef BasicIteratedExpr<expr::COUNT, NumericExpr, LogicalExpr> CountExpr;
MP_SPECIALIZE_IS(CountExpr, COUNT)

// A numberof expression.
// Example: numberof 42 in ({i in I} x[i]),
// where I is a set and x is a variable.
typedef BasicIteratedExpr<expr::NUMBEROF> NumberOfExpr;
MP_SPECIALIZE_IS(NumberOfExpr, NUMBEROF)

// A logical constant.
// Examples: 0, 1
class LogicalConstant : public LogicalExpr {
 private:
  struct Impl : Expr::Impl {
    bool value;
  };
  MP_EXPR;

 public:
  // Returns the value of this constant.
  bool value() const { return impl()->value; }
};

MP_SPECIALIZE_IS(LogicalConstant, CONSTANT)

// A logical NOT expression.
// Example: not a, where a is a logical expression.
typedef BasicUnaryExpr<LogicalExpr> NotExpr;
MP_SPECIALIZE_IS(NotExpr, NOT)

// A binary logical expression.
// Examples: a || b, a && b, where a and b are logical expressions.
typedef BasicBinaryExpr<LogicalExpr> BinaryLogicalExpr;
MP_SPECIALIZE_IS_RANGE(BinaryLogicalExpr, BINARY_LOGICAL)

// A relational expression.
// Examples: x < y, x != y, where x and y are variables.
typedef BasicBinaryExpr<LogicalExpr, NumericExpr> RelationalExpr;
MP_SPECIALIZE_IS_RANGE(RelationalExpr, RELATIONAL)

// A logical count expression.
// Examples: atleast 1 (x < y, x != y), where x and y are variables.
class LogicalCountExpr : public LogicalExpr {
 private:
  struct Impl : Expr::Impl {
    const Expr::Impl *lhs;
    const Expr::Impl *rhs;
  };
  MP_EXPR;

 public:
  // Returns the left-hand side (the first argument) of this expression.
  NumericExpr lhs() const { return Create<NumericExpr>(impl()->lhs); }

  // Returns the right-hand side (the second argument) of this expression.
  CountExpr rhs() const { return Create<CountExpr>(impl()->rhs); }
};

MP_SPECIALIZE_IS_RANGE(LogicalCountExpr, LOGICAL_COUNT)

// An implication expression.
// Example: a ==> b else c, where a, b and c are logical expressions.
typedef BasicIfExpr<LogicalExpr, expr::IMPLICATION> ImplicationExpr;
MP_SPECIALIZE_IS(ImplicationExpr, IMPLICATION)

// An iterated logical expression.
// Example: exists{i in I} x[i] >= 0, where I is a set and x is a variable.
typedef BasicIteratedExpr<
  expr::FIRST_ITERATED_LOGICAL, LogicalExpr> IteratedLogicalExpr;
MP_SPECIALIZE_IS_RANGE(IteratedLogicalExpr, ITERATED_LOGICAL)

// TODO: expressions: alldiff and string literal

class ExprFactory {
 private:
  std::vector<Function::Impl*> funcs_;
  std::vector<Expr::Impl*> exprs_;

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
    impl->kind_ = kind;
    exprs_.back() = impl;
    return impl;
  }

  template <typename ExprType, typename Arg>
  ExprType MakeUnary(expr::Kind kind, Arg arg) {
    MP_ASSERT(arg != 0, "invalid argument");
    typename ExprType::Impl *impl = Allocate<typename  ExprType::Impl>(kind);
    impl->arg = arg.impl_;
    return Expr::Create<ExprType>(impl);
  }

  template <typename ExprType, typename LHS, typename RHS>
  ExprType MakeBinary(expr::Kind kind, LHS lhs, RHS rhs) {
    MP_ASSERT(internal::Is<ExprType>(kind), "invalid expression kind");
    MP_ASSERT(lhs != 0 && rhs != 0, "invalid argument");
    typename ExprType::Impl *impl = Allocate<typename ExprType::Impl>(kind);
    impl->lhs = lhs.impl_;
    impl->rhs = rhs.impl_;
    return Expr::Create<ExprType>(impl);
  }

  template <typename ExprType, typename Arg>
  ExprType MakeIf(LogicalExpr condition, Arg true_expr, Arg false_expr) {
    // false_expr can be null.
    MP_ASSERT(condition != 0 && true_expr != 0, "invalid argument");
    typename ExprType::Impl *impl =
        Allocate<typename ExprType::Impl>(ExprType::KIND);
    impl->condition = condition.impl_;
    impl->true_expr = true_expr.impl_;
    impl->false_expr = false_expr.impl_;
    return Expr::Create<ExprType>(impl);
  }

  // A variable argument expression builder.
  template <typename ExprType>
  class IteratedExprBuilder {
   private:
    typename ExprType::Impl *impl_;
    int arg_index_;

    friend class ExprFactory;

    explicit IteratedExprBuilder(typename ExprType::Impl *impl)
      : impl_(impl), arg_index_(0) {}

   public:
    void AddArg(typename ExprType::Arg arg) {
      MP_ASSERT(arg_index_ < impl_->num_args, "too many arguments");
      MP_ASSERT(arg != 0, "invalid argument");
      impl_->args[arg_index_++] = arg.impl_;
    }
  };

  template <typename ExprType>
  IteratedExprBuilder<ExprType> BeginIteratedExpr(
        expr::Kind kind, int num_args) {
    MP_ASSERT(num_args >= 0, "invalid number of arguments");
    typename ExprType::Impl *impl = Allocate<typename ExprType::Impl>(
          kind, sizeof(Expr::Impl*) * (num_args - 1));
    impl->num_args = num_args;
    return IteratedExprBuilder<ExprType>(impl);
  }

  template <typename ExprType>
  ExprType EndIteratedExpr(IteratedExprBuilder<ExprType> builder) {
    typename ExprType::Impl *impl = builder.impl_;
    // Check that all arguments provided.
    MP_ASSERT(builder.arg_index_ == impl->num_args, "too few arguments");
    return Expr::Create<ExprType>(impl);
  }

 public:
  ExprFactory() {}
  ~ExprFactory() {
    // TODO: delete expressions
  }

  Function AddFunction(const char *) {
    // Call push_back first to make sure that the impl pointer doesn't leak
    // if push_back throws an exception.
    funcs_.push_back(0);
    Function::Impl *impl = new Function::Impl();
    funcs_.back() = impl;
    return Function(impl);
  }

  // Makes a numeric constant.
  NumericConstant MakeNumericConstant(double value) {
    NumericConstant::Impl *impl =
        Allocate<NumericConstant::Impl>(expr::CONSTANT);
    impl->value = value;
    return Expr::Create<NumericConstant>(impl);
  }

  // Makes a variable reference.
  Variable MakeVariable(int index) {
    Variable::Impl *impl = Allocate<Variable::Impl>(expr::VARIABLE);
    impl->index = index;
    return Expr::Create<Variable>(impl);
  }

  // Makes a unary expression.
  UnaryExpr MakeUnary(expr::Kind kind, NumericExpr arg) {
    MP_ASSERT(internal::Is<UnaryExpr>(kind), "invalid expression kind");
    return MakeUnary<UnaryExpr>(kind, arg);
  }

  // Makes a binary expression.
  BinaryExpr MakeBinary(expr::Kind kind, NumericExpr lhs, NumericExpr rhs) {
    return MakeBinary<BinaryExpr>(kind, lhs, rhs);
  }

  // Makes an if expression.
  IfExpr MakeIf(LogicalExpr condition,
                NumericExpr true_expr, NumericExpr false_expr) {
    return MakeIf<IfExpr>(condition, true_expr, false_expr);
  }

  // A piecewise-linear term builder.
  class PLTermBuilder {
   private:
    PLTerm::Impl *impl_;
    int slope_index_;
    int breakpoint_index_;

    friend class ExprFactory;

    explicit PLTermBuilder(PLTerm::Impl *impl)
      : impl_(impl), slope_index_(0), breakpoint_index_(0) {}

   public:
    void AddSlope(double slope) {
      MP_ASSERT(slope_index_ < impl_->num_breakpoints + 1, "too many slopes");
      impl_->data[2 * slope_index_] = slope;
      ++slope_index_;
    }

    void AddBreakpoint(double breakpoint) {
      MP_ASSERT(breakpoint_index_ < impl_->num_breakpoints,
                "too many breakpoints");
      impl_->data[2 * breakpoint_index_ + 1] = breakpoint;
      ++breakpoint_index_;
    }
  };

  // Begins building a piecewise-linear term.
  PLTermBuilder BeginPLTerm(int num_breakpoints) {
    MP_ASSERT(num_breakpoints > 0, "invalid number of breakpoints");
    PLTerm::Impl *impl = Allocate<PLTerm::Impl>(
          expr::PLTERM, sizeof(double) * num_breakpoints * 2);
    impl->num_breakpoints = num_breakpoints;
    return PLTermBuilder(impl);
  }

  // Ends building a piecewise-linear term.
  PLTerm EndPLTerm(PLTermBuilder builder, Variable var) {
    PLTerm::Impl *impl = builder.impl_;
    // Check that all slopes and breakpoints provided.
    MP_ASSERT(builder.slope_index_ == impl->num_breakpoints + 1,
              "too few slopes");
    MP_ASSERT(builder.breakpoint_index_ == impl->num_breakpoints,
              "too few breakpoints");
    MP_ASSERT(var != 0, "invalid argument");
    impl->var_index = var.index();
    return Expr::Create<PLTerm>(impl);
  }

  typedef IteratedExprBuilder<CallExpr> CallExprBuilder;

  // Begins building a call expression.
  CallExprBuilder BeginCall(Function func, int num_args) {
    MP_ASSERT(func != 0, "invalid function");
    CallExprBuilder builder = BeginIteratedExpr<CallExpr>(expr::CALL, num_args);
    builder.impl_->func = func.impl_;
    return builder;
  }

  // Ends building a call expression.
  CallExpr EndCall(CallExprBuilder builder) {
    return EndIteratedExpr<CallExpr>(builder);
  }

  typedef IteratedExprBuilder<VarArgExpr> VarArgExprBuilder;

  // Begins building a variable argument expression.
  VarArgExprBuilder BeginVarArg(expr::Kind kind, int num_args) {
    MP_ASSERT(internal::Is<VarArgExpr>(kind), "invalid expression kind");
    return BeginIteratedExpr<VarArgExpr>(kind, num_args);
  }

  // Ends building a variable argument expression.
  VarArgExpr EndVarArg(VarArgExprBuilder builder) {
    return EndIteratedExpr<VarArgExpr>(builder);
  }

  typedef IteratedExprBuilder<SumExpr> SumExprBuilder;

  // Begins building a sum expression.
  SumExprBuilder BeginSum(int num_args) {
    return BeginIteratedExpr<SumExpr>(expr::SUM, num_args);
  }

  // Ends building a sum expression.
  SumExpr EndSum(SumExprBuilder builder) {
    return EndIteratedExpr<SumExpr>(builder);
  }

  typedef IteratedExprBuilder<CountExpr> CountExprBuilder;

  // Begins building a count expression.
  CountExprBuilder BeginCount(int num_args) {
    return BeginIteratedExpr<CountExpr>(expr::COUNT, num_args);
  }

  // Ends building a count expression.
  CountExpr EndCount(CountExprBuilder builder) {
    return EndIteratedExpr<CountExpr>(builder);
  }

  typedef IteratedExprBuilder<NumberOfExpr> NumberOfExprBuilder;

  // Begins building a numberof expression.
  NumberOfExprBuilder BeginNumberOf(int num_args, NumericExpr arg0) {
    MP_ASSERT(num_args >= 1, "invalid number of arguments");
    NumberOfExprBuilder builder =
        BeginIteratedExpr<NumberOfExpr>(expr::NUMBEROF, num_args);
    builder.AddArg(arg0);
    return builder;
  }

  // Ends building a numberof expression.
  NumberOfExpr EndNumberOf(NumberOfExprBuilder builder) {
    return EndIteratedExpr<NumberOfExpr>(builder);
  }

  // Makes a logical constant.
  LogicalConstant MakeLogicalConstant(bool value) {
    LogicalConstant::Impl *impl =
        Allocate<LogicalConstant::Impl>(expr::CONSTANT);
    impl->value = value;
    return Expr::Create<LogicalConstant>(impl);
  }

  // Makes a logical NOT expression.
  NotExpr MakeNot(LogicalExpr arg) {
    return MakeUnary<NotExpr>(expr::NOT, arg);
  }

  // Makes a binary logical expression.
  BinaryLogicalExpr MakeBinaryLogical(
        expr::Kind kind, LogicalExpr lhs, LogicalExpr rhs) {
    return MakeBinary<BinaryLogicalExpr>(kind, lhs, rhs);
  }

  // Makes a relational expression.
  RelationalExpr MakeRelational(
        expr::Kind kind, NumericExpr lhs, NumericExpr rhs) {
    return MakeBinary<RelationalExpr>(kind, lhs, rhs);
  }

  // Makes a logical count expression.
  LogicalCountExpr MakeLogicalCount(
        expr::Kind kind, NumericExpr lhs, CountExpr rhs) {
    return MakeBinary<LogicalCountExpr>(kind, lhs, rhs);
  }

  // Makes an implication expression.
  ImplicationExpr MakeImplication(LogicalExpr condition, LogicalExpr true_expr,
                                  LogicalExpr false_expr) {
    return MakeIf<ImplicationExpr>(condition, true_expr, false_expr);
  }

  typedef IteratedExprBuilder<IteratedLogicalExpr> IteratedLogicalExprBuilder;

  // Begins building a sum expression.
  IteratedLogicalExprBuilder BeginIteratedLogical(
        expr::Kind kind, int num_args) {
    return BeginIteratedExpr<IteratedLogicalExpr>(kind, num_args);
  }

  // Ends building a sum expression.
  IteratedLogicalExpr EndIteratedLogical(IteratedLogicalExprBuilder builder) {
    return EndIteratedExpr<IteratedLogicalExpr>(builder);
  }
};
}  // namespace mp

#endif  // MP_EXPR_H_
