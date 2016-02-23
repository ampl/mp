/*
 ASL expression classes

 Copyright (C) 2012 AMPL Optimization Inc

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

#ifndef MP_ASLEXPR_H_
#define MP_ASLEXPR_H_

#include <cassert>
#include <cstddef>
#include <algorithm>
#include <iterator>
#include <map>
#include <string>
#include <utility>
#include <vector>

extern "C" {
#include "nlp.h"
}

// Undefine ASL macros because they often clash with names used in solver
// libraries.
#undef solve_code
#undef var_name
#undef con_name
#undef filename
#undef ampl_vbtol
#undef strtod
#ifndef ASL_PRESERVE_DEFINES
# undef Char
# undef A_colstarts
# undef A_rownos
# undef A_vals
# undef Cgrad
# undef Fortran
# undef LUrhs
# undef LUv
# undef Lastx
# undef Ograd
# undef Urhsx
# undef Uvx
# undef X0
# undef adjoints
# undef adjoints_nv1
# undef amax
# undef ampl_options
# undef amplflag
# undef archan
# undef awchan
# undef binary_nl
# undef c_cexp1st
# undef c_vars
# undef co_index
# undef comb
# undef combc
# undef comc
# undef comc1
# undef como
# undef como1
# undef cv_index
# undef cvar
# undef err_jmp
# undef err_jmp1
# undef fhash
# undef funcs
# undef funcsfirst
# undef funcslast
# undef havepi0
# undef havex0
# undef lnc
# undef maxcolnamelen
# undef maxrownamelen
# undef n_cc
# undef n_con
# undef n_conjac
# undef n_obj
# undef n_var
# undef nbv
# undef nclcon
# undef ncom0
# undef ncom1
# undef nderps
# undef need_nl
# undef n_eqn
# undef nfunc
# undef niv
# undef nlc
# undef nlcc
# undef nlnc
# undef nlo
# undef n_lcon
# undef nlogv
# undef nlvb
# undef nlvbi
# undef nlvc
# undef nlvci
# undef nlvo
# undef nlvoi
# undef nranges
# undef nwv
# undef nzc
# undef nzjac
# undef nzo
# undef o_cexp1st
# undef o_vars
# undef obj_no
# undef objtype
# undef pi0
# undef plterms
# undef real
# undef return_nofile
# undef size_expr_n
# undef skip_int_derivs
# undef sputinfo
# undef stub_end
# undef want_deriv
# undef want_xpi0
# undef x0kind
# undef x0len
# undef xscanf
# undef zaC
# undef zac
# undef zao
# undef zerograds
# undef write_sol
#endif  // ASL_PRESERVE_DEFINES

#include "mp/common.h"
#include "mp/error.h"

namespace mp {

class ASLProblem;

namespace asl {

template <expr::Kind FIRST, expr::Kind LAST>
class BasicExpr;

typedef BasicExpr<expr::FIRST_EXPR, expr::LAST_EXPR> Expr;

// A numeric expression.
typedef BasicExpr<expr::FIRST_NUMERIC, expr::LAST_NUMERIC> NumericExpr;

// A logical or constraint expression.
typedef BasicExpr<expr::FIRST_LOGICAL, expr::LAST_LOGICAL> LogicalExpr;

class LogicalConstant;
}  // namespace asl

namespace internal {
template <>
inline bool Is<asl::LogicalExpr>(expr::Kind k) {
  return (k >= expr::FIRST_LOGICAL && k <= expr::LAST_LOGICAL) ||
          k == expr::NUMBER;
}

template <>
inline bool Is<asl::LogicalConstant>(expr::Kind k) {
  return k == expr::BOOL || k == expr::NUMBER;
}
}
}  // namespace mp

namespace mp {

namespace internal {
// Casts expression to type ExprType.
// If assertions are enabled, it generates an assertion failure when
// e is not of runtime type ExprType. Otherwise no runtime check is
// performed.
template <typename ExprType>
ExprType UncheckedCast(asl::Expr e);
}

namespace asl {

namespace internal {

class ASLBuilder;

template <typename ExprType>
class ExprIterator;

class ExprBase {
 private:
  typedef ::expr Impl;

  void True() const {}

  template <typename ExprType>
  friend class internal::ExprIterator;

  // Constructs an Expr object from the implementation impl. Only a minimal
  // check is performed when assertions are enabled to make sure that the
  // opcode is within the valid range.
  explicit ExprBase(Impl *impl) : impl_(impl) {
    assert(!impl_ || (kind() >= expr::FIRST_EXPR && kind() <= expr::LAST_EXPR));
  }

 protected:
  Impl *impl_;

  // Returns a pointer to the implementation.
  Impl *impl() const { return impl_; }

  template <typename ExprType>
  ExprBase(ExprType e) : impl_(e.impl_) {}

  // An expression proxy used for implementing operator-> in iterators.
  template <typename ExprType>
  class Proxy {
   private:
    ExprType impl_;

   public:
    explicit Proxy(Impl *e) : impl_(Create<ExprType>(e)) {}

    const ExprType *operator->() const { return &impl_; }
  };

  // Creates an expression object from a raw expr pointer.
  // For safety reason expression classes don't provide constructors
  // taking raw pointers and this method should be used instead.
  template <typename ExprType>
  static ExprType Create(Impl *e) {
    assert(!e || mp::internal::Is<ExprType>(ExprBase(e).kind()));
    ExprType expr;
    expr.impl_ = e;
    return expr;
  }

  typedef void (ExprBase::*SafeBool)() const;

 public:
  // Constructs an Expr object representing a null reference to an AMPL
  // expression. The only operation permitted for such expression is
  // copying, assignment and check whether it is null using operator SafeBool.
  ExprBase() : impl_() {}

  // Returns a value convertible to bool that can be used in conditions but not
  // in comparisons and evaluates to "true" if this expression is not null
  // and "false" otherwise.
  // Example:
  //   void foo(Expr e) {
  //     if (e) {
  //       // Do something if e is not null.
  //     }
  //   }
  operator SafeBool() const { return impl_ ? &ExprBase::True : 0; }

  // Returns the expression kind.
  expr::Kind kind() const {
    std::size_t opcode = reinterpret_cast<std::size_t>(impl_->op);
    return mp::internal::GetOpCodeInfo(static_cast<int>(opcode)).kind;
  }
};
}  // namespace internal

template <typename Impl, typename Result, typename LResult>
class ExprConverter;

// An expression.
// An Expr object represents a reference to an expression so
// it is cheap to construct and pass by value. A type safe way to
// process expressions of different types is by using ExprVisitor.
template <expr::Kind FIRST, expr::Kind LAST = FIRST>
class BasicExpr : private internal::ExprBase {
  friend class internal::ExprBase;
  friend class mp::ASLProblem;
  friend class internal::ASLBuilder;

  template <typename Impl, typename Result, typename LResult>
  friend class ExprConverter;

  template <typename ExprType>
  friend ExprType mp::internal::UncheckedCast(asl::Expr e);

 protected:
  using ExprBase::impl;
  using ExprBase::Create;
  using ExprBase::Proxy;

 public:
  enum { FIRST_KIND = FIRST, LAST_KIND = LAST };

  BasicExpr() {}

  template <typename Expr>
  BasicExpr(
      Expr other,
      typename fmt::internal::EnableIf<
        static_cast<expr::Kind>(Expr::FIRST_KIND) >= FIRST &&
        static_cast<expr::Kind>(Expr::LAST_KIND)  <= LAST, int>::type = 0)
    : ExprBase(other) {}

  using ExprBase::kind;
  using ExprBase::operator SafeBool;

  bool operator==(BasicExpr rhs) const { return impl() == rhs.impl(); }
  bool operator!=(BasicExpr rhs) const { return impl() != rhs.impl(); }
};

// Casts an expression to type ExprType. Returns a null expression if the cast
// is not possible.
template <typename ExprType>
ExprType Cast(Expr e) {
  return mp::internal::Is<ExprType>(e.kind()) ?
        mp::internal::UncheckedCast<ExprType>(e) : ExprType();
}

// A numeric constant.
// Examples: 42, -1.23e-4
class NumericConstant : public BasicExpr<expr::NUMBER> {
 public:
  // Returns the value of this number.
  double value() const { return reinterpret_cast<const expr_n*>(impl())->v; }
};

// A reference to a variable or a common expression.
// Example: x
class Reference :
  public BasicExpr<expr::FIRST_REFERENCE, expr::LAST_REFERENCE> {
 public:
  // Returns the index of the referenced object.
  int index() const { return impl()->a; }
};

template <typename Arg, expr::Kind FIRST, expr::Kind LAST = FIRST>
class BasicUnaryExpr : public BasicExpr<FIRST, LAST> {
 public:
  // Returns the argument of this expression.
  Arg arg() const {
    return BasicUnaryExpr::template Create<Arg>(this->impl()->L.e);
  }
};

// A unary numeric expression.
// Examples: -x, sin(x), where x is a variable.
typedef BasicUnaryExpr<
  NumericExpr, expr::FIRST_UNARY, expr::LAST_UNARY> UnaryExpr;

// A binary expression.
// Base: base expression class.
// Arg: argument expression class.
template <typename Arg, expr::Kind FIRST, expr::Kind LAST = FIRST>
class BasicBinaryExpr : public BasicExpr<FIRST, LAST> {
 public:
  // Returns the left-hand side (the first argument) of this expression.
  Arg lhs() const {
    return BasicBinaryExpr::template Create<Arg>(this->impl()->L.e);
  }

  // Returns the right-hand side (the second argument) of this expression.
  Arg rhs() const {
    return BasicBinaryExpr::template Create<Arg>(this->impl()->R.e);
  }
};

// A binary numeric expression.
// Examples: x / y, atan2(x, y), where x and y are variables.
typedef BasicBinaryExpr<
  NumericExpr, expr::FIRST_BINARY, expr::LAST_BINARY> BinaryExpr;

template <typename Arg, expr::Kind KIND>
class BasicIfExpr : public BasicExpr<KIND> {
 private:
  typedef BasicExpr<KIND> Base;

  const expr_if *impl() const {
    return reinterpret_cast<const expr_if*>(Base::impl());
  }

 public:
  LogicalExpr condition() const {
    return Base::template Create<LogicalExpr>(impl()->e);
  }

  Arg then_expr() const { return Base::template Create<Arg>(impl()->T); }
  Arg else_expr() const { return Base::template Create<Arg>(impl()->F); }
};

// An if-then-else expression.
// Example: if x != 0 then y else z, where x, y and z are variables.
typedef BasicIfExpr<NumericExpr, expr::IF> IfExpr;

// A piecewise-linear expression.
// Example: <<0; -1, 1>> x, where x is a variable.
class PiecewiseLinearExpr : public BasicExpr<expr::PLTERM> {
 public:
  // Returns the number of breakpoints in this term.
  int num_breakpoints() const {
    return num_slopes() - 1;
  }

  // Returns the number of slopes in this term.
  int num_slopes() const {
    assert(impl()->L.p->n >= 1);
    return impl()->L.p->n;
  }

  // Returns a breakpoint with the specified index.
  double breakpoint(int index) const {
    assert(index >= 0 && index < num_breakpoints());
    return impl()->L.p->bs[2 * index + 1];
  }

  // Returns a slope with the specified index.
  double slope(int index) const {
    assert(index >= 0 && index < num_slopes());
    return impl()->L.p->bs[2 * index];
  }

  NumericExpr arg() const {
    return Create<NumericExpr>(impl()->R.e);
  }
};

class CallExpr;

class Function {
 private:
  func_info *fi_;

  friend class CallExpr;
  friend class internal::ASLBuilder;

  explicit Function(func_info *fi) : fi_(fi) {}

  void True() const {}
  typedef void (Function::*SafeBool)() const;

 public:
  Function() : fi_(0) {}

  // Returns the function name.
  const char *name() const { return fi_->name; }

  // Returns the number of arguments.
  int num_args() const { return fi_->nargs; }

  // Returns a value convertible to bool that can be used in conditions but not
  // in comparisons and evaluates to "true" if this function is not null
  // and "false" otherwise.
  // Example:
  //   void foo(Function f) {
  //     if (f) {
  //       // Do something if e is not null.
  //     }
  //   }
  operator SafeBool() const { return fi_ ? &Function::True : 0; }

  bool operator==(const Function &other) const { return fi_ == other.fi_; }
  bool operator!=(const Function &other) const { return fi_ != other.fi_; }
};

namespace internal {
// An expression iterator.
template <typename ExprType>
class ExprIterator :
  public std::iterator<std::forward_iterator_tag, ExprType> {
 private:
  ::expr *const *ptr_;

 public:
  explicit ExprIterator(::expr *const *p = 0) : ptr_(p) {}

  ExprType operator*() const { return ExprBase::Create<ExprType>(*ptr_); }

  ExprBase::Proxy<ExprType> operator->() const {
    return ExprBase::Proxy<ExprType>(*ptr_);
  }

  ExprIterator &operator++() {
    ++ptr_;
    return *this;
  }

  ExprIterator operator++(int ) {
    ExprIterator it(*this);
    ++ptr_;
    return it;
  }

  bool operator==(ExprIterator other) const { return ptr_ == other.ptr_; }
  bool operator!=(ExprIterator other) const { return ptr_ != other.ptr_; }
};
}  // namespace internal

// A function call expression.
// Example: f(x), where f is a function and x is a variable.
class CallExpr : public BasicExpr<expr::CALL> {
 private:
  const expr_f *impl() const {
    return reinterpret_cast<const expr_f*>(BasicExpr<expr::CALL>::impl());
  }

 public:
  Function function() const { return Function(impl()->fi); }

  int num_args() const { return impl()->al->n; }

  Expr arg(int index) {
    assert(index >= 0 && index < num_args());
    return Create<Expr>(impl()->args[index]);
  }

  // An argument iterator.
  typedef internal::ExprIterator<Expr> iterator;

  iterator begin() const { return iterator(impl()->args); }
  iterator end() const { return iterator(impl()->args + num_args()); }
};

// A numeric expression with a variable number of arguments.
// The min and max functions always have at least one argument.
// Example: min{i in I} x[i], where I is a set and x is a variable.
class VarArgExpr : public BasicExpr<expr::FIRST_VARARG, expr::LAST_VARARG> {
 private:
  static const de END;

 public:
  typedef NumericExpr Arg;

  // An argument iterator.
  class iterator :
    public std::iterator<std::forward_iterator_tag, NumericExpr> {
   private:
    const de *de_;

    friend class VarArgExpr;

    explicit iterator(const de *d) : de_(d) {}

   public:
    iterator() : de_(&END) {}

    NumericExpr operator*() const { return Create<NumericExpr>(de_->e); }

    Proxy<NumericExpr> operator->() const {
      return Proxy<NumericExpr>(de_->e);
    }

    iterator &operator++() {
      ++de_;
      return *this;
    }

    iterator operator++(int ) {
      iterator it(*this);
      ++de_;
      return it;
    }

    bool operator==(iterator other) const { return de_->e == other.de_->e; }
    bool operator!=(iterator other) const { return de_->e != other.de_->e; }
  };

  iterator begin() const {
    return iterator(reinterpret_cast<const expr_va*>(impl())->L.d);
  }

  iterator end() const {
    return iterator();
  }
};

// An iterated expression.
// ID is used to distinguish between expression types that have the same
// base and argument type, namely SumExpr and NumberOfExpr.
template <typename ArgType, expr::Kind FIRST, expr::Kind LAST = FIRST>
class BasicIteratedExpr : public BasicExpr<FIRST, LAST> {
 public:
  typedef ArgType Arg;

  // Returns the number of arguments.
  int num_args() const {
    return static_cast<int>(this->impl()->R.ep - this->impl()->L.ep);
  }

  Arg arg(int index) const {
    assert(index >= 0 && index < num_args());
    return BasicIteratedExpr::template Create<Arg>(this->impl()->L.ep[index]);
  }

  typedef internal::ExprIterator<Arg> iterator;

  iterator begin() const { return iterator(this->impl()->L.ep); }
  iterator end() const { return iterator(this->impl()->R.ep); }
};

// A sum expression.
// Example: sum{i in I} x[i], where I is a set and x is a variable.
typedef BasicIteratedExpr<NumericExpr, expr::SUM> SumExpr;

// A count expression.
// Example: count{i in I} (x[i] >= 0), where I is a set and x is a variable.
typedef BasicIteratedExpr<LogicalExpr, expr::COUNT> CountExpr;

// A numberof expression.
// Example: numberof 42 in ({i in I} x[i]),
// where I is a set and x is a variable.
typedef BasicIteratedExpr<NumericExpr, expr::NUMBEROF> NumberOfExpr;

// A symbolic numberof expression.
// Example: numberof (if x then 'a' else 'b') in ('a', 'b', 'c'),
// where x is a variable.
typedef BasicIteratedExpr<Expr, expr::NUMBEROF_SYM> SymbolicNumberOfExpr;

// A logical constant.
// Examples: 0, 1
class LogicalConstant : public BasicExpr<expr::BOOL> {
 public:
  // Returns the value of this constant.
  bool value() const { return reinterpret_cast<const expr_n*>(impl())->v != 0; }
};

// A logical NOT expression.
// Example: not a, where a is a logical expression.
typedef BasicUnaryExpr<LogicalExpr, expr::NOT> NotExpr;

// A binary logical expression.
// Examples: a || b, a && b, where a and b are logical expressions.
typedef BasicBinaryExpr<LogicalExpr,
  expr::FIRST_BINARY_LOGICAL, expr::LAST_BINARY_LOGICAL> BinaryLogicalExpr;

// A relational expression.
// Examples: x < y, x != y, where x and y are variables.
typedef BasicBinaryExpr<
  NumericExpr, expr::FIRST_RELATIONAL, expr::LAST_RELATIONAL> RelationalExpr;

// A logical count expression.
// Examples: atleast 1 (x < y, x != y), where x and y are variables.
class LogicalCountExpr :
  public BasicExpr<expr::FIRST_LOGICAL_COUNT, expr::LAST_LOGICAL_COUNT> {
 public:
  // Returns the left-hand side (the first argument) of this expression.
  NumericExpr lhs() const { return Create<NumericExpr>(impl()->L.e); }

  // Returns the right-hand side (the second argument) of this expression.
  CountExpr rhs() const { return Create<CountExpr>(impl()->R.e); }
};

// An implication expression.
// Example: a ==> b else c, where a, b and c are logical expressions.
typedef BasicIfExpr<LogicalExpr, expr::IMPLICATION> ImplicationExpr;

// An iterated logical expression.
// Example: exists{i in I} x[i] >= 0, where I is a set and x is a variable.
typedef BasicIteratedExpr<LogicalExpr,
  expr::FIRST_ITERATED_LOGICAL, expr::LAST_ITERATED_LOGICAL>
  IteratedLogicalExpr;

// A pairwise expression.
// Example: alldiff{i in I} x[i], where I is a set and x is a variable.
typedef BasicIteratedExpr<
  NumericExpr, expr::FIRST_PAIRWISE, expr::LAST_PAIRWISE> PairwiseExpr;

class StringLiteral : public BasicExpr<expr::STRING> {
 public:
  const char *value() const {
    return reinterpret_cast<const expr_h*>(impl())->sym;
  }
};

// A symbolic if-then-else expression.
// Example: if x != 0 then 'a' else 0, where x is a variable.
typedef BasicIfExpr<Expr, expr::IFSYM> SymbolicIfExpr;

template <typename Grad>
class LinearExpr;

// A linear term.
template <typename GradT>
class LinearTerm {
 private:
  typedef GradT Grad;
  Grad *grad_;

  friend class LinearExpr< LinearTerm<Grad> >;

  explicit LinearTerm(Grad *g = 0) : grad_(g) {}

 public:
  // Returns the coefficient.
  double coef() const { return grad_->coef; }

  // Returns the variable index.
  int var_index() const { return grad_->varno; }
};

typedef LinearTerm<ograd> LinearObjTerm;
typedef LinearTerm<cgrad> LinearConTerm;

// A linear expression.
template <typename Term>
class LinearExpr {
 private:
  Term first_term_;

  friend class mp::ASLProblem;

  explicit LinearExpr(typename Term::Grad *first_term)
  : first_term_(Term(first_term)) {}

 public:
  LinearExpr() {}

  class iterator : public std::iterator<std::forward_iterator_tag, Term> {
   private:
    Term term_;

   public:
    explicit iterator(Term t) : term_(t) {}

    Term operator*() const { return term_; }
    const Term *operator->() const { return &term_; }

    iterator &operator++() {
      term_ = Term(term_.grad_->next);
      return *this;
    }

    iterator operator++(int ) {
      iterator it(*this);
      term_ = Term(term_.grad_->next);
      return it;
    }

    bool operator==(iterator other) const {
      return term_.grad_ == other.term_.grad_;
    }
    bool operator!=(iterator other) const {
      return term_.grad_ != other.term_.grad_;
    }
  };

  iterator begin() const { return iterator(first_term_); }
  iterator end() const { return iterator(Term()); }
};

typedef LinearExpr<LinearObjTerm> LinearObjExpr;
typedef LinearExpr<LinearConTerm> LinearConExpr;

class LinearCommonExpr {
 private:
  const cexp *expr_;

 public:
  explicit LinearCommonExpr(const cexp *e = 0) : expr_(e) {}
  int num_terms() const { return expr_->nlin; }
  // TODO
};

template <typename LinearExpr>
void WriteExpr(fmt::Writer &w, LinearExpr linear, NumericExpr nonlinear);

namespace internal {

// ASL expression types.
struct ExprTypes {
  typedef mp::asl::Expr Expr;
  typedef mp::asl::NumericExpr NumericExpr;
  typedef mp::asl::LogicalExpr LogicalExpr;
  typedef mp::asl::NumericConstant NumericConstant;
  typedef mp::asl::Reference Variable;
  typedef mp::asl::Reference CommonExpr;
  typedef mp::asl::UnaryExpr UnaryExpr;
  typedef mp::asl::BinaryExpr BinaryExpr;
  typedef mp::asl::IfExpr IfExpr;
  typedef mp::asl::PiecewiseLinearExpr PLTerm;
  typedef mp::asl::CallExpr CallExpr;
  typedef mp::asl::VarArgExpr VarArgExpr;
  typedef mp::asl::SumExpr SumExpr;
  typedef mp::asl::NumberOfExpr NumberOfExpr;
  typedef mp::asl::SymbolicNumberOfExpr SymbolicNumberOfExpr;
  typedef mp::asl::CountExpr CountExpr;
  typedef mp::asl::LogicalConstant LogicalConstant;
  typedef mp::asl::NotExpr NotExpr;
  typedef mp::asl::BinaryLogicalExpr BinaryLogicalExpr;
  typedef mp::asl::RelationalExpr RelationalExpr;
  typedef mp::asl::LogicalCountExpr LogicalCountExpr;
  typedef mp::asl::ImplicationExpr ImplicationExpr;
  typedef mp::asl::IteratedLogicalExpr IteratedLogicalExpr;
  typedef mp::asl::PairwiseExpr PairwiseExpr;
  typedef mp::asl::StringLiteral StringLiteral;
  typedef mp::asl::SymbolicIfExpr SymbolicIfExpr;

  // Checked cast. See mp::Cast.
  template <typename ExprType>
  static ExprType Cast(Expr e) {
    return mp::asl::Cast<ExprType>(e);
  }

  // Unchecked cast. See mp::internal::Cast.
  template <typename ExprType>
  static ExprType UncheckedCast(Expr e) {
    return mp::internal::UncheckedCast<ExprType>(e);
  }
};
}  // namespace internal
}  // namespace asl
}  // namespace mp

template <typename ExprType>
ExprType mp::internal::UncheckedCast(asl::Expr e) {
  assert(Is<ExprType>(e.kind()));
  ExprType expr;
  expr.impl_ = e.impl_;
  return expr;
}

#endif  // MP_ASLEXPR_H_
