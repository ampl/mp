/*
 AMPL function testing infrastructure.

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

#ifndef TESTS_FUNCTION_H_
#define TESTS_FUNCTION_H_

#include <algorithm>
#include <iosfwd>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>

struct ASL;
struct func_info;

namespace fun {

// A tuple of doubles.
class Tuple {
 private:
  std::vector<double> items_;

  Tuple &operator<<(double arg) {
    items_.push_back(arg);
    return *this;
  }

 public:
  static Tuple GetTupleWithSize(unsigned size) {
    Tuple t(0);
    t.items_.resize(size);
    return t;
  }

  explicit Tuple(double a0) { *this << a0; }
  Tuple(double a0, double a1) { *this << a0 << a1; }
  Tuple(double a0, double a1, double a2) { *this << a0 << a1 << a2; }
  Tuple(double a0, double a1, double a2, double a3) {
    *this << a0 << a1 << a2 << a3;
  }
  Tuple(double a0, double a1, double a2, double a3, double a4) {
    *this << a0 << a1 << a2 << a3 << a4;
  }
  Tuple(double a0, double a1, double a2, double a3, double a4, double a5) {
    *this << a0 << a1 << a2 << a3 << a4 << a5;
  }
  Tuple(double a0, double a1, double a2, double a3,
      double a4, double a5, double a6, double a7, double a8) {
    *this << a0 << a1 << a2 << a3 << a4 << a5 << a6 << a7 << a8;
  }

  unsigned size() const { return items_.size(); }
  double &operator[](unsigned index) { return items_.at(index); }
  double operator[](unsigned index) const { return items_.at(index); }
};

std::ostream &operator<<(std::ostream &os, const Tuple &t);

// A dynamic bit set.
class BitSet {
 private:
  std::vector<bool> store_;

 public:
  typedef std::vector<bool>::reference reference;
  typedef std::vector<bool>::const_reference const_reference;

  BitSet() {}

  BitSet(unsigned size, bool value) : store_(size, value) {}

  explicit BitSet(const char *s);

  unsigned size() const { return store_.size(); }
  reference operator[](unsigned index) { return store_.at(index); }
  const_reference operator[](unsigned index) const { return store_.at(index); }
};

enum Type { VOID, INT, UINT, DOUBLE };

template <typename T>
struct GetType;

template <>
struct GetType<void> {
  static const Type VALUE = VOID;
};

template <>
struct GetType<int> {
  static const Type VALUE = INT;
};

template <>
struct GetType<unsigned> {
  static const Type VALUE = UINT;
};

template <>
struct GetType<double> {
  static const Type VALUE = DOUBLE;
};

template <typename Arg1, typename Arg2 = void,
    typename Arg3 = void, typename Arg4 = void>
class FunctionWithTypes {
 private:
  static const Type ARG_TYPES[];

 protected:
  void CheckArgs(const Tuple &args) const {
    if (args.size() != GetNumArgs())
      throw std::invalid_argument("invalid number of arguments");
  }

 public:
  unsigned GetNumArgs() const;

  Type GetArgType(unsigned arg_index) const {
    if (arg_index >= GetNumArgs())
      throw std::out_of_range("argument index is out of range");
    return ARG_TYPES[arg_index];
  }
};

template <typename Arg1, typename Arg2, typename Arg3, typename Arg4>
unsigned FunctionWithTypes<Arg1, Arg2, Arg3, Arg4>::GetNumArgs() const {
  if (GetType<Arg2>::VALUE == VOID)
    return 1;
  if (GetType<Arg3>::VALUE == VOID)
    return 2;
  if (GetType<Arg4>::VALUE == VOID)
    return 3;
  return 4;
}

template <typename Arg1, typename Arg2, typename Arg3, typename Arg4>
const Type FunctionWithTypes<Arg1, Arg2, Arg3, Arg4>::ARG_TYPES[] = {
    GetType<Arg1>::VALUE, GetType<Arg2>::VALUE,
    GetType<Arg3>::VALUE, GetType<Arg4>::VALUE
};

template <typename Arg1, typename Result>
class FunctionPointer1 : public FunctionWithTypes<Arg1> {
 private:
  Result (*f_)(Arg1);

 public:
  explicit FunctionPointer1(Result (*f)(Arg1)) : f_(f) {}

  Result operator()(const Tuple &args) const {
    this->CheckArgs(args);
    return f_(args[0]);
  }
};

template <typename Arg1, typename Result>
FunctionPointer1<Arg1, Result> FunctionPointer(Result (*f)(Arg1)) {
  return FunctionPointer1<Arg1, Result>(f);
}

template <typename Arg1, typename Arg2, typename Result>
class FunctionPointer2 : public FunctionWithTypes<Arg1, Arg2> {
 private:
  Result (*f_)(Arg1, Arg2);

 public:
  explicit FunctionPointer2(Result (*f)(Arg1, Arg2)) : f_(f) {}

  Result operator()(const Tuple &args) const {
    this->CheckArgs(args);
    return f_(args[0], args[1]);
  }
};

template <typename Arg1, typename Arg2, typename Result>
FunctionPointer2<Arg1, Arg2, Result> FunctionPointer(Result (*f)(Arg1, Arg2)) {
  return FunctionPointer2<Arg1, Arg2, Result>(f);
}

template <typename Arg1, typename Arg2, typename Arg3, typename Result>
class FunctionPointer3 : public FunctionWithTypes<Arg1, Arg2, Arg3> {
 private:
  Result (*f_)(Arg1, Arg2, Arg3);

 public:
  explicit FunctionPointer3(Result (*f)(Arg1, Arg2, Arg3)) : f_(f) {}

  Result operator()(const Tuple &args) const {
    this->CheckArgs(args);
    return f_(args[0], args[1], args[2]);
  }
};

template <typename Arg1, typename Arg2, typename Arg3, typename Result>
FunctionPointer3<Arg1, Arg2, Arg3, Result>
  FunctionPointer(Result (*f)(Arg1, Arg2, Arg3)) {
  return FunctionPointer3<Arg1, Arg2, Arg3, Result>(f);
}

template <typename Arg1, typename Arg2,
  typename Arg3, typename Arg4, typename Result>
class FunctionPointer4 : public FunctionWithTypes<Arg1, Arg2, Arg3, Arg4> {
 private:
  Result (*f_)(Arg1, Arg2, Arg3, Arg4);

 public:
  explicit FunctionPointer4(Result (*f)(Arg1, Arg2, Arg3, Arg4)) : f_(f) {}

  Result operator()(const Tuple &args) const {
    this->CheckArgs(args);
    return f_(args[0], args[1], args[2], args[3]);
  }
};

template <typename Arg1, typename Arg2,
  typename Arg3, typename Arg4, typename Result>
FunctionPointer4<Arg1, Arg2, Arg3, Arg4, Result>
  FunctionPointer(Result (*f)(Arg1, Arg2, Arg3, Arg4)) {
  return FunctionPointer4<Arg1, Arg2, Arg3, Arg4, Result>(f);
}

// A functor class with one argument bound.
template <typename F, typename Arg>
class OneBinder {
 private:
  F f_;
  Arg value_;
  unsigned bound_arg_index_;

 public:
  OneBinder(F f, Arg value, unsigned bound_arg_index);

  unsigned GetNumArgs() const { return f_.GetNumArgs() - 1; }

  Type GetArgType(unsigned arg_index) const {
    if (arg_index >= bound_arg_index_)
      ++arg_index;
    return f_.GetArgType(arg_index);
  }

  double operator()(const Tuple &args) const;
};

template <typename F, typename Arg>
OneBinder<F, Arg>::OneBinder(F f, Arg value, unsigned bound_arg_index)
: f_(f), value_(value), bound_arg_index_(bound_arg_index) {
  if (bound_arg_index >= f.GetNumArgs())
    throw std::out_of_range("argument index is out of range");
}

template <typename F, typename Arg>
double OneBinder<F, Arg>::operator()(const Tuple &args) const {
  unsigned num_args = args.size();
  Tuple part_args(Tuple::GetTupleWithSize(num_args + 1));
  for (unsigned i = 0; i < bound_arg_index_; ++i)
    part_args[i] = args[i];
  for (unsigned i = bound_arg_index_; i < num_args; ++i)
    part_args[i + 1] = args[i];
  part_args[bound_arg_index_] = value_;
  return f_(part_args);
}

// Binds one arguments.
template <typename F, typename Arg>
OneBinder<F, Arg> BindOne(F f, Arg value, unsigned bound_arg_index) {
  return OneBinder<F, Arg>(f, value, bound_arg_index);
}

// A functor class with all but one arguments bound.
template <typename F>
class AllButOneBinder {
 private:
  F f_;
  mutable Tuple args_;
  unsigned unbound_arg_index_;

 public:
  AllButOneBinder(F f, const Tuple &args, unsigned unbound_arg_index);

  double operator()(double x) const {
    args_[unbound_arg_index_] = x;
    return f_(args_);
  }
};

template <typename F>
AllButOneBinder<F>::AllButOneBinder(
    F f, const Tuple &args, unsigned unbound_arg_index)
: f_(f), args_(args), unbound_arg_index_(unbound_arg_index) {
  if (unbound_arg_index >= args.size())
    throw std::out_of_range("argument index is out of range");
}

// Binds all but one arguments.
template <typename F>
AllButOneBinder<F> BindAllButOne(
    F f, const Tuple &args, unsigned unbound_arg_index) {
  return AllButOneBinder<F>(f, args, unbound_arg_index);
}

// A utility class for computing the derivative by Ridders' method
// of polynomial extrapolation. The implementation is taken from
// "Numerical Recipes in C", Chapter 5.7.
class Differentiator {
 private:
  enum {NTAB = 200};

  // Successive columns in the Neville tableau will go to smaller
  // step sizes and higher orders of extrapolation.
  std::vector<double> table_;

  // Returns the element at row i and column j of the Neville tableau.
  double &at(unsigned i, unsigned j) {
    return table_[i * NTAB + j];
  }

  template <typename F>
  static double SymmetricDifference(F f, double x, double h) {
    return (f(x + h) - f(x - h)) / (2 * h);
  }

  template <typename F>
  static double LeftDifference(F f, double x, double h) {
    return (f(x) - f(x - h)) / h;
  }

  template <typename F>
  static double RightDifference(F f, double x, double h) {
    return (f(x + h) - f(x)) / h;
  }

  // Returns the derivative of a function f at a point x by Ridders'
  // method of polynomial extrapolation.
  template <typename F, typename D>
  double Differentiate(F f, double x, D d, double *error = 0);

 public:
  Differentiator() : table_(NTAB * NTAB) {}

  // Returns the derivative of a function f at a point x by Ridders'
  // method of polynomial extrapolation trying to detect indeterminate case.
  template <typename F>
  double operator()(F f, double x, double *error = 0, bool *detected_nan = 0);
};

template <typename F, typename D>
double Differentiator::Differentiate(F f, double x, D d, double *error) {
  const double CON = 1.4, CON2 = CON * CON;
  double safe = 20;
  double hh = 0.125, ans = std::numeric_limits<double>::quiet_NaN(), diff = 0;
  at(0, 0) = d(f, x, hh);

  // If at(0, 0) is NaN try reducing hh a couple of times.
  for (unsigned i = 0; i < 10 && std::isnan(at(0, 0)); i++) {
    hh /= CON;
    at(0, 0) = d(f, x, hh);
  }

  double err = std::numeric_limits<double>::max();
  unsigned i = 1;
  for (; i < NTAB; i++, safe *= 0.95) {
    // Try new, smaller step size.
    hh /= CON;
    at(0, i) = d(f, x, hh);
    double fac = CON2;
    for (unsigned j = 1; j <= i; j++) {
      // Compute extrapolations of various orders, requiring no new function
      // evaluations.
      at(j, i) = (at(j - 1, i) * fac - at(j - 1, i - 1)) / (fac - 1);
      fac = CON2 * fac;
      double errt = std::max(
          std::fabs(at(j, i) - at(j - 1, i)),
          std::fabs(at(j, i) - at(j - 1, i - 1)));
      // The error strategy is to compare each new extrapolation to one order
      // lower, both at the present step size and the previous one.
      if (errt <= err) {
        // If error is decreased, save the improved answer.
        err = errt;
        ans = at(j, i);
      }
    }
    // If higher order is worse by a significant factor 'safe', then quit early.
    diff = std::fabs(at(i, i) - at(i - 1, i - 1));
    if (diff >= safe * err || std::isnan(diff))
      break;
  }

#ifdef DEBUG_DIFFERENTIATOR
  std::cout << "deriv=" << ans << " err=" << err << " iter=" << i
      << " diff=" << diff << std::endl;
#endif

  if (error)
    *error = err;
  return ans;
}

template <typename F>
double Differentiator::operator()(
    F f, double x, double *error, bool *detected_nan) {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  if (detected_nan)
    *detected_nan = false;
  if (std::isnan(f(x)))
    return nan;
  double dummy_error = 0;
  if (!error)
    error = &dummy_error;
  double deriv = Differentiate(f, x, SymmetricDifference<F>, error);
  double right_error = 0;
  double right_deriv = Differentiate(f, x, RightDifference<F>, &right_error);
  if (std::isnan(deriv)) {
    *error = right_error;
    return right_deriv;
  }
  double left_error = nan;
  double left_deriv = Differentiate(f, x, LeftDifference<F>, &left_error);
  if ((!(std::fabs(left_deriv - right_deriv) <= 1e-2) &&
      left_error / (std::fabs(left_deriv) + 1) < 0.05) ||
          (std::isnan(left_deriv) && std::isnan(right_deriv))) {
    if (detected_nan)
      *detected_nan = true;
    return nan;
  }
  return deriv;
}

class Function;

// Function information that can't be obtained automatically, in particular
// due to limitations of numerical differentiation.
class FunctionInfo {
 private:
  std::vector<std::string> arg_names_;

 public:
  explicit FunctionInfo(const char *arg_names = "") {
    SetArgNames(arg_names);
  }
  virtual ~FunctionInfo();

  std::string GetArgName(unsigned index) const {
    return arg_names_.at(index);
  }

  // Sets argument names separated by spaces.
  FunctionInfo &SetArgNames(const char *arg_names);

  class Result {
   private:
    double value_;
    std::string error_;

   public:
    explicit Result(double value) : value_(value) {}
    explicit Result(const char *error = "") :
      value_(std::numeric_limits<double>::quiet_NaN()), error_(error) {}

    double value() const { return value_; }
    const char *error() const { return error_.empty() ? 0 : error_.c_str(); }
  };

  // Returns the value of the derivative at the point specified by args.
  // In most cases this is not needed as the derivative is computed using
  // numerical differentiation. This function can also return an error message.
  virtual Result GetDerivative(
      const Function &f, unsigned arg_index, const Tuple &args) const;

  // Returns the value of the second derivative at the point specified by args.
  // In most cases this is not needed as the derivative is computed using
  // numerical differentiation. This function can also return an error message.
  virtual Result GetSecondDerivative(const Function &f,
      unsigned arg1_index, unsigned arg2_index, const Tuple &args) const;
};

// Flags for an AMPL function call.
enum {
  DERIVS = 1,  // Get first partial derivatives.
  HES    = 3   // Get both first and second partial derivatives.
};

// An AMPL function.
class Function {
 private:
  ASL *asl_;
  func_info *fi_;
  const FunctionInfo *info_;

 public:
  Function(ASL *asl, func_info *fi, const FunctionInfo *info) :
    asl_(asl), fi_(fi), info_(info) {}

  const char *name() const;
  int nargs() const;

  const FunctionInfo *info() const { return info_; }

  // A result of an AMPL function call.
  class Result {
   private:
    double value_;
    std::vector<double> derivs_;
    std::vector<double> hes_;
    const char *error_;

    void CheckError() const {
      if (error_)
        throw std::runtime_error(error_);
    }

   public:
    Result(double value, const std::vector<double> &derivs,
        const std::vector<double> &hes, const char *error) :
      value_(value), derivs_(derivs), hes_(hes), error_(error) {}

    operator double() const {
      CheckError();
      return value_;
    }

    double deriv(size_t index = 0) const {
      CheckError();
      return derivs_.at(index);
    }
    double hes(size_t index = 0) const {
      CheckError();
      return hes_.at(index);
    }

    const char *error() const { return error_; }
  };

  // Calls a function.
  Result operator()(const Tuple &args, int flags = 0,
      const BitSet &use_deriv = BitSet(), void *info = 0) const;

  std::string GetArgName(unsigned index) const {
    return info_->GetArgName(index);
  }

  FunctionInfo::Result GetDerivative(
      unsigned arg_index, const Tuple &args) const {
    return info_->GetDerivative(*this, arg_index, args);
  }

  FunctionInfo::Result GetSecondDerivative(
      unsigned arg1_index, unsigned arg2_index, const Tuple &args) const {
    return info_->GetSecondDerivative(*this, arg1_index, arg2_index, args);
  }
};

// A helper class that wraps a Function's derivative and binds
// all but one argument to the given values.
class DerivativeBinder {
 private:
  Function f_;
  unsigned deriv_arg_;
  unsigned eval_arg_;
  Tuple args_;
  BitSet use_deriv_;

 public:
  // Creates a Derivative object.
  // deriv_arg: index of an argument with respect to which
  //            the derivative is taken
  // eval_arg:  index of an argument which is not bound
  DerivativeBinder(Function f, unsigned deriv_arg,
      unsigned eval_arg, const Tuple &args);

  double operator()(double x);
};
}

#endif  // TESTS_FUNCTION_H_
