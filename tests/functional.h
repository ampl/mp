// Function adapters, binders, numerical differentiator and other
// function-related functionality.

#ifndef TESTS_FUNCTIONAL_H
#define TESTS_FUNCTIONAL_H

#include <algorithm>
#include <limits>
#include <vector>
#include <cmath>

#ifdef DEBUG_DIFFERENTIATOR
# include <iostream>
#endif

namespace fun {

// A base class for ternary function objects.
template <typename Arg1, typename Arg2, typename Arg3, typename Result>
struct ternary_function {
  typedef Arg1 first_argument_type;
  typedef Arg2 second_argument_type;
  typedef Arg3 third_argument_type;
  typedef Result result_type;
};

// Adapter of a ternary function pointer to the ternary_function interface.
template <typename Arg1, typename Arg2, typename Arg3, typename Result>
class pointer_to_ternary_function :
  public ternary_function<Arg1, Arg2, Arg3, Result>
{
private:
  Result (*f_)(Arg1, Arg2, Arg3);

public:
  explicit pointer_to_ternary_function(Result (*f)(Arg1, Arg2, Arg3)) :
    f_(f) {}
  Result operator()(Arg1 x, Arg2 y, Arg3 z) const { return f_(x, y, z); }
};

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
  double operator()(F f, double x, double *error, bool *detected_nan);
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

}

#endif // TESTS_FUNCTIONAL_H
