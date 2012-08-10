/*
 Tests of the AMPL bindings for GNU Scientific Library.

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

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>

#include <functional>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <cstring>

#include "gtest/gtest.h"
#include "tests/config.h"
// #define DEBUG_DIFFERENTIATOR
#include "tests/function.h"
#include "solvers/asl.h"

using std::string;
using std::vector;

using fun::BitSet;
using fun::DERIVS;
using fun::Differentiator;
using fun::Function;
using fun::FunctionInfo;
using fun::HES;
using fun::Tuple;
using fun::pointer_to_ternary_function;
using fun::ternary_function;

namespace {

// Converts error estimate returned by Diff into an absolute tolerance to
// be used in EXPECT_NEAR.
double ConvertErrorToTolerance(double error) {
  return error != 0 ? error * 1000 : 1e-10;
}

class Error {
 private:
  string str;

 public:
  explicit Error(const string &s) : str(s) {}
  operator const char*() const { return str.c_str(); }
  const char* c_str() const { return str.c_str(); }
  operator FunctionInfo::Result() const {
    return FunctionInfo::Result(str.c_str());
  }
};

Error EvalError(const Function &f, const Tuple &args, const char *suffix = "") {
  std::ostringstream os;
  os << "can't evaluate " << f.name() << suffix << args;
  return Error(os.str());
}

Error NotIntError(const string &arg_name, double value = 0.5) {
  std::ostringstream os;
  os << "argument '" << arg_name
      << "' can't be represented as int, " << arg_name << " = " << value;
  return Error(os.str());
}

#define EXPECT_ERROR(expected_message, result) \
  EXPECT_STREQ(expected_message, (result).error())

// Check if the value returned by af is correct.
void CheckFunction(double value, const Function &f, const Tuple &args) {
  std::ostringstream os;
  os << "Checking if " << f.name() << args << " = " << value;
  SCOPED_TRACE(os.str());
  if (gsl_isnan(value))
    EXPECT_ERROR(EvalError(f, args), f(args));
  else
    EXPECT_EQ(value, f(args)) << f.name() << args;
}

// A helper class that wraps a Function's derivative and binds
// all but one argument to the given values.
class DerivativeBinder {
 private:
  Function f_;
  unsigned deriv_var_;
  unsigned eval_var_;
  Tuple args_;
  BitSet use_deriv_;

 public:
  // Creates a Derivative object.
  // deriv_var: index of a variable with respect to which
  //            the derivative is taken
  // eval_var:  index of a variable which is not bound
  DerivativeBinder(Function f, unsigned deriv_var,
      unsigned eval_var, const Tuple &args);

  double operator()(double x);
};

DerivativeBinder::DerivativeBinder(Function f, unsigned deriv_var,
    unsigned eval_var, const Tuple &args)
: f_(f), deriv_var_(deriv_var), eval_var_(eval_var),
  args_(args), use_deriv_(args.size(), false) {
  unsigned num_vars = args_.size();
  if (deriv_var >= num_vars || eval_var >= num_vars)
    throw std::out_of_range("variable index is out of range");
  use_deriv_[deriv_var] = true;
}

double DerivativeBinder::operator()(double x) {
  args_[eval_var_] = x;
  Function::Result r = f_(args_, DERIVS, use_deriv_);
  return r.error() ? GSL_NAN : r.deriv(deriv_var_);
}

typedef double (*FuncU)(unsigned);
typedef double (*Func3Mode)(double, double, double, gsl_mode_t);
typedef double (*FuncND)(int, double);

const double POINTS[] = {-5, -2, -1.23, -1, 0, 1, 1.23, 2, 5};
const size_t NUM_POINTS = sizeof(POINTS) / sizeof(*POINTS);

const double POINTS_FOR_N[] = {-2, -1, 0, 1, 2};
const size_t NUM_POINTS_FOR_N = sizeof(POINTS_FOR_N) / sizeof(*POINTS_FOR_N);

class GSLTest : public ::testing::Test {
 protected:
  ASL *asl;
  FunctionInfo info;  // Default function info.

  Differentiator diff;

  // Differentiator statistics.
  struct Stats {
    unsigned num_calls_;
    unsigned num_nans_;
    double min_error_;
    double max_error_;

    Stats() :
      num_calls_(0), num_nans_(0),
      min_error_(std::numeric_limits<double>::max()),
      max_error_(std::numeric_limits<double>::min()) {}

    ~Stats() {
      std::cout << "Called numerical differentiation "
          << num_calls_ << " times with " << num_nans_ << " "
          << (num_nans_ == 1 ? "NaN" : "NaNs")
          << " detected" << std::endl;
      std::cout << "Error min=" << min_error_;
      std::cout << "  max=" << max_error_ << std::endl;
    }
  };
  static Stats stats_;

  template <typename F>
  double Diff(F f, double x, double *error = 0) {
    ++stats_.num_calls_;
    double err = GSL_NAN;
    bool detected_nan = false;
    double deriv = diff(f, x, &err, &detected_nan);
    if (!gsl_isnan(deriv)) {
      stats_.min_error_ = std::min(stats_.min_error_, err);
      stats_.max_error_ = std::max(stats_.max_error_, err);
    }
    if (detected_nan)
      ++stats_.num_nans_;
    if (error)
      *error = err;
    return deriv;
  }

  void SetUp() {
    asl = ASL_alloc(ASL_read_f);
    i_option_ASL = "../solvers/gsl/libamplgsl.so";
    func_add(asl);
  }

  void TearDown() {
    ASL_free(&asl);
  }

  // Returns an AMPL function by name.
  Function GetFunction(const char *name, FunctionInfo *info = 0) const {
    func_info *fi = func_lookup(asl, name, 0);
    if (!fi)
      throw std::runtime_error(string("function not found: ") + name);
    return Function(asl, fi, info);
  }

  template <typename F>
  bool CheckDerivative(F f, const Function &af,
      unsigned var_index, const Tuple &args);

  static const unsigned NO_VAR = ~0u;

  void CheckSecondDerivatives(const Function &f,
      const Tuple &args, unsigned skip_var = NO_VAR);

  // Tests a function taking a single argument.
  template <typename F>
  void TestUnaryFunc(const Function &af, F f);
  void TestFunc(const Function &af, double (*f)(double x)) {
    TestUnaryFunc(af, f);
  }
  void TestFunc(const Function &af, double (*f)(double, gsl_mode_t)) {
    TestUnaryFunc(af, std::bind2nd(std::ptr_fun(f), GSL_PREC_DOUBLE));
  }

  void TestZeroFunc(const Function &af,
      double value, const Tuple &args, unsigned s_index);

  // Tests a zero function.
  void TestFunc(const Function &af, FuncU f) {
    for (size_t i = 0; i != NUM_POINTS; ++i) {
      double s = POINTS[i];
      TestZeroFunc(af, f(s), Tuple(s), 0);
    }
  }

  // Tests a function taking an integer and a double parameter.
  // test_x is a value of x where the function can be computed for very large
  // and very small n. If there is no such x or it is not known, then test_x
  // should be GSL_NAN.
  void TestFuncND(const Function &af, FuncND f,
      double test_x, const string &arg_name);

  template <typename F>
  void TestBinaryFunc(const Function &af, F f);

  void TestFunc(const Function &af, double (*f)(double, double)) {
    TestBinaryFunc(af, std::ptr_fun(f));
  }

  typedef double (*Func2Mode)(double, double, gsl_mode_t);

  // Binds the mode argument of Func2Mode to GSL_PREC_DOUBLE.
  class Func2DoubleMode : public std::binary_function<double, double, double> {
   private:
    Func2Mode f_;

   public:
    explicit Func2DoubleMode(Func2Mode f) : f_(f) {}

    double operator()(double x, double y) const {
      return f_(x, y, GSL_PREC_DOUBLE);
    }
  };

  void TestFunc(const Function &af, Func2Mode f) {
    TestBinaryFunc(af, Func2DoubleMode(f));
  }

  // Binds the mode argument of Func2Mode to GSL_PREC_DOUBLE.
  class Func3DoubleMode :
    public ternary_function<double, double, double, double> {
   private:
    Func3Mode f_;

   public:
    explicit Func3DoubleMode(Func3Mode f) : f_(f) {}

    double operator()(double x, double y, double z) const {
      return f_(x, y, z, GSL_PREC_DOUBLE);
    }
  };

  template <typename F>
  void TestTernaryFunc(const Function &af, F f);

  template <typename Arg1, typename Arg2, typename Arg3, typename Result>
  void TestFunc(const Function &af, Result (*f)(Arg1, Arg2, Arg3)) {
    TestTernaryFunc(af,
        pointer_to_ternary_function<Arg1, Arg2, Arg3, Result>(f));
  }
  void TestFunc(const Function &af, Func3Mode f) {
    TestTernaryFunc(af, Func3DoubleMode(f));
  }
};

GSLTest::Stats GSLTest::stats_;

// Checks if the value of the derivative returned by af agrees with the
// value returned by numerical differentiation of function f.
// var_index: index of the variable with respect to which to differentiate
// args: point at which the derivative is computed
template <typename F>
bool GSLTest::CheckDerivative(F f, const Function &af,
    unsigned var_index, const Tuple &args) {
  std::ostringstream os;
  os << "Checking d/dx" << var_index << " " << af.name() << " at " << args;
  SCOPED_TRACE(os.str());

  BitSet use_deriv(args.size(), false);
  use_deriv[var_index] = true;

  double error = 0;
  double x = args[var_index];
  FunctionInfo::Result deriv_result = af.GetDerivative(var_index, args);
  if (!deriv_result.error()) {
    double numerical_deriv = Diff(f, x, &error);
    double overridden_deriv = deriv_result.value();
    if (!gsl_isnan(overridden_deriv) && overridden_deriv != numerical_deriv) {
      std::cout << "Overriding d/dx" << var_index << " " << af.name()
        << " at " << args << ", computed = " << numerical_deriv
        << ", overridden = " << overridden_deriv << std::endl;
      numerical_deriv = overridden_deriv;
    }
    if (!gsl_isnan(numerical_deriv)) {
      double deriv = af(args, DERIVS, use_deriv).deriv(var_index);
      if (numerical_deriv != deriv)
        EXPECT_NEAR(numerical_deriv, deriv, ConvertErrorToTolerance(error));
      return true;
    }
  }

  Function::Result r = af(args, DERIVS, use_deriv);
  if (!gsl_isnan(f(x))) {
    if (deriv_result.error())
      EXPECT_ERROR(deriv_result.error(), r);
    else
      EXPECT_ERROR(EvalError(af, args, "'"), r);
  } else {
    EXPECT_TRUE(r.error() != nullptr);
  }
  return false;
}

// Checks if the values of the second partial derivatives returned by af
// agree with the values returned by numerical differentiation of the first
// partial derivatives.
// args: point at which the derivatives are computed
// skip_var: index of the variable with respect to which not to differentiate
void GSLTest::CheckSecondDerivatives(const Function &f,
    const Tuple &args, unsigned skip_var) {
  unsigned num_args = args.size();
  if (skip_var == NO_VAR) {
    for (unsigned i = 0; i < num_args; ++i) {
      if (f.GetDerivative(i, args).error()) {
        skip_var = i;
        break;
      }
    }
  }
  for (unsigned i = 0; i < num_args; ++i) {
    if (i == skip_var) continue;
    for (unsigned j = 0; j < num_args; ++j) {
      if (j == skip_var) continue;
      BitSet use_deriv(num_args, false);
      use_deriv[i] = true;
      use_deriv[j] = true;
      double error = 0;
      FunctionInfo::Result deriv_result = f.GetSecondDerivative(i, j, args);
      if (!deriv_result.error()) {
        double d = Diff(DerivativeBinder(f, j, i, args), args[i], &error);
        double overridden_deriv = deriv_result.value();
        if (!gsl_isnan(overridden_deriv) && overridden_deriv != d) {
          std::cout << "Overriding d/dx" << i << " d/dx" << j << " "
            << f.name() << " at " << args << ", computed = " << d
            << ", overridden = " << overridden_deriv << std::endl;
          d = overridden_deriv;
        }
        std::ostringstream os;
        os << "Checking if d/dx" << i << " d/dx" << j
          << " " << f.name() << " at " << args << " is " << d;
        SCOPED_TRACE(os.str());
        if (!gsl_isnan(d)) {
          unsigned ii = i, jj = j;
          if (ii > jj) std::swap(ii, jj);
          unsigned hes_index = ii * (2 * num_args - ii - 1) / 2 + jj;
          double actual_deriv = f(args, HES, use_deriv).hes(hes_index);
          if (d != actual_deriv)
            EXPECT_NEAR(d, actual_deriv, ConvertErrorToTolerance(error));
          continue;
        }
      }
      Function::Result r = f(args, HES, use_deriv);
      if (f(args, DERIVS, use_deriv).error())
        EXPECT_TRUE(r.error() != nullptr);
      else if (deriv_result.error())
        EXPECT_ERROR(deriv_result.error(), r);
      else
        EXPECT_ERROR(EvalError(f, args, "''"), r);
    }
  }
}

template <typename F>
void GSLTest::TestUnaryFunc(const Function &af, F f) {
  for (size_t i = 0; i != NUM_POINTS; ++i) {
    double x = POINTS[i];
    Tuple arg(x);
    CheckFunction(f(x), af, arg);
    CheckDerivative(f, af, 0, arg);
    CheckSecondDerivatives(af, arg);
  }
}

void GSLTest::TestZeroFunc(const Function &af,
    double value, const Tuple &args, unsigned s_index) {
  double s = args[s_index];
  if (static_cast<unsigned>(s) == s) {
    if (gsl_isnan(value)) {
      EXPECT_ERROR(EvalError(af, args), af(args));
      EXPECT_ERROR(EvalError(af, args), af(args, DERIVS));
      EXPECT_ERROR(EvalError(af, args), af(args, HES));
    } else {
      EXPECT_EQ(value, af(args)) << af.name() << " at " << args;
      EXPECT_ERROR(EvalError(af, args, "'"), af(args, DERIVS));
      EXPECT_ERROR(EvalError(af, args, "'"), af(args, HES));
    }
  } else {
    std::ostringstream os;
    os << "argument 's' can't be represented as unsigned int, s = " << s;
    EXPECT_ERROR(os.str().c_str(), af(args));
    EXPECT_ERROR(os.str().c_str(), af(args, DERIVS));
    EXPECT_ERROR(os.str().c_str(), af(args, HES));
  }
}

void GSLTest::TestFuncND(const Function &af, FuncND f,
    double test_x, const string &arg_name) {
  for (size_t i = 0; i != NUM_POINTS_FOR_N; ++i) {
    int n = POINTS_FOR_N[i];
    for (size_t j = 0; j != NUM_POINTS; ++j) {
      double x = POINTS[j];
      Tuple args(n, x);
      CheckFunction(f(n, x), af, args);
      string error("argument '" + arg_name + "' is not constant");
      EXPECT_ERROR(error.c_str(), af(args, DERIVS));
      EXPECT_ERROR(error.c_str(), af(args, HES));
      CheckDerivative(std::bind1st(std::ptr_fun(f), n), af, 1, args);
      CheckSecondDerivatives(af, args, 0);
    }
  }
  EXPECT_ERROR(("argument '" + arg_name + "' can't be represented as int, " +
      arg_name + " = 0.5").c_str(), af(Tuple(0.5, 0)));

  if (gsl_isnan(test_x))
    return;

  // These points are tested separately because of various problems, e.g.
  // gsl_sf_bessel_Jn(n, x) and gsl_sf_bessel_In(n, x) take too much time
  // (hang?) for n = INT_MIN and gsl_sf_bessel_Yn(n, x) returns different
  // values close to 0 when called different times for n = INT_MIN.

  BitSet use_deriv("01");
  EXPECT_ERROR(
      ("can't compute derivative: argument '" + arg_name + "' too small, " +
      arg_name + " = -2147483648").c_str(),
      af(Tuple(INT_MIN, test_x), DERIVS, use_deriv));
  EXPECT_ERROR(
      ("can't compute derivative: argument '" + arg_name + "' too large, " +
      arg_name + " = 2147483647").c_str(),
      af(Tuple(INT_MAX, test_x), DERIVS, use_deriv));
  EXPECT_TRUE(!gsl_isnan(af(Tuple(INT_MIN + 1, test_x), DERIVS,
      use_deriv).deriv(1)));
  EXPECT_TRUE(!gsl_isnan(af(Tuple(INT_MAX - 1, test_x), DERIVS,
      use_deriv).deriv(1)));

  EXPECT_ERROR(
      ("can't compute derivative: argument '" + arg_name + "' too small, " +
      arg_name + " = -2147483647").c_str(),
      af(Tuple(INT_MIN + 1, test_x), HES, use_deriv));
  EXPECT_ERROR(
      ("can't compute derivative: argument '" + arg_name + "' too large, " +
      arg_name + " = 2147483646").c_str(),
      af(Tuple(INT_MAX - 1, test_x), HES, use_deriv));
  EXPECT_TRUE(!gsl_isnan(af(Tuple(INT_MIN + 2, test_x), HES,
      use_deriv).hes(2)));
  EXPECT_TRUE(!gsl_isnan(af(Tuple(INT_MAX - 2, test_x), HES,
      use_deriv).hes(2)));
}

template <typename F>
void GSLTest::TestBinaryFunc(const Function &af, F f) {
  for (size_t i = 0; i != NUM_POINTS; ++i) {
    for (size_t j = 0; j != NUM_POINTS; ++j) {
      double x = POINTS[i], y = POINTS[j];
      Tuple args(x, y);
      CheckFunction(f(x, y), af, args);
      CheckDerivative(std::bind2nd(f, y), af, 0, args);
      CheckDerivative(std::bind1st(f, x), af, 1, args);
      CheckSecondDerivatives(af, args);
    }
  }
}

// Checks that the error is returned when trying to get a derivative
// with respect to an integer argument.
// Returns true if the argument is integer, false otherwise.
template <typename Arg>
bool CheckArg(const Function &, const Tuple &, unsigned);

template <>
bool CheckArg<double>(const Function &, const Tuple &, unsigned) {
  return false;
}

template <>
bool CheckArg<int>(const Function &f, const Tuple &args, unsigned arg_index) {
  std::ostringstream os;
  BitSet use_deriv(args.size(), false);
  use_deriv[arg_index] = true;
  os << "argument '" << f.GetArgName(arg_index) << "' is not constant";
  EXPECT_STREQ(os.str().c_str(), f(args, DERIVS, use_deriv).error());
  return true;
}

template <typename F>
void GSLTest::TestTernaryFunc(const Function &af, F f) {
  for (size_t i = 0; i != NUM_POINTS; ++i) {
    for (size_t j = 0; j != NUM_POINTS; ++j) {
      for (size_t k = 0; k != NUM_POINTS; ++k) {
        double x = POINTS[i], y = POINTS[j], z = POINTS[k];
        Tuple args(x, y, z);
        if (static_cast<typename F::first_argument_type>(x) != x) {
          EXPECT_STREQ(NotIntError(af.GetArgName(0), x), af(args).error());
          continue;
        }
        if (static_cast<typename F::second_argument_type>(y) != y) {
          EXPECT_STREQ(NotIntError(af.GetArgName(1), y), af(args).error());
          continue;
        }
        CheckFunction(f(x, y, z), af, args);
        if (!CheckArg<typename F::first_argument_type>(af, args, 0))
          CheckDerivative(Bind2Of3(f, args, 0), af, 0, args);
        if (!CheckArg<typename F::second_argument_type>(af, args, 1))
          CheckDerivative(Bind2Of3(f, args, 1), af, 1, args);
        if (!CheckArg<typename F::third_argument_type>(af, args, 2))
          CheckDerivative(Bind2Of3(f, args, 2), af, 2, args);
        CheckSecondDerivatives(af, args);
      }
    }
  }
}

#define TEST_FUNC(name) TestFunc(GetFunction("gsl_" #name, &info), gsl_##name)

#define TEST_FUNC_ND(name, test_x, arg) \
  TestFuncND(GetFunction("gsl_" #name, &info), gsl_##name, test_x, #arg)

double ellint_E(double x) {
  return gsl_sf_ellint_E(-1.23, x, GSL_PREC_DOUBLE);
}

TEST_F(GSLTest, DerivativeBinder) {
  DerivativeBinder d(GetFunction("gsl_hypot"), 0, 1, Tuple(1, 0));
  ASSERT_EQ(1, d(0));
  ASSERT_EQ(1 / sqrt(2), d(1));
  d = DerivativeBinder(GetFunction("gsl_hypot"), 1, 1, Tuple(1, 0));
  ASSERT_EQ(0, d(0));
  ASSERT_EQ(1 / sqrt(2), d(1));
  EXPECT_THROW(
      DerivativeBinder(GetFunction("gsl_hypot"), 2, 0, Tuple(0, 0)),
      std::out_of_range);
  EXPECT_THROW(
      DerivativeBinder(GetFunction("gsl_hypot"), 0, 2, Tuple(0, 0)),
      std::out_of_range);
}

TEST_F(GSLTest, Elementary) {
  TEST_FUNC(log1p);
  TEST_FUNC(expm1);
  TEST_FUNC(hypot);
  TEST_FUNC(hypot3);
}

TEST_F(GSLTest, AiryA) {
  TEST_FUNC(sf_airy_Ai);
  TEST_FUNC(sf_airy_Ai_scaled);
}

TEST_F(GSLTest, AiryB) {
  TEST_FUNC(sf_airy_Bi);
  TEST_FUNC(sf_airy_Bi_scaled);
}

TEST_F(GSLTest, AiryZero) {
  TEST_FUNC(sf_airy_zero_Ai);
  TEST_FUNC(sf_airy_zero_Bi);
  TEST_FUNC(sf_airy_zero_Ai_deriv);
  TEST_FUNC(sf_airy_zero_Bi_deriv);
}

TEST_F(GSLTest, BesselJ) {
  TEST_FUNC(sf_bessel_J0);
  TEST_FUNC(sf_bessel_J1);
  TEST_FUNC_ND(sf_bessel_Jn, 0, n);
}

TEST_F(GSLTest, BesselY) {
  TEST_FUNC(sf_bessel_Y0);
  TEST_FUNC(sf_bessel_Y1);
  TEST_FUNC_ND(sf_bessel_Yn, 1, n);
}

TEST_F(GSLTest, BesselI) {
  TEST_FUNC(sf_bessel_I0);
  TEST_FUNC(sf_bessel_I1);
  TEST_FUNC_ND(sf_bessel_In, 0, n);
  TEST_FUNC(sf_bessel_I0_scaled);
  TEST_FUNC(sf_bessel_I1_scaled);
  TEST_FUNC_ND(sf_bessel_In_scaled, 0, n);
}

TEST_F(GSLTest, BesselK) {
  TEST_FUNC(sf_bessel_K0);
  TEST_FUNC(sf_bessel_K1);
  TEST_FUNC_ND(sf_bessel_Kn, 1, n);
  TEST_FUNC(sf_bessel_K0_scaled);
  TEST_FUNC(sf_bessel_K1_scaled);
  TEST_FUNC_ND(sf_bessel_Kn_scaled, 1, n);
}

TEST_F(GSLTest, Besselj) {
  TEST_FUNC(sf_bessel_j0);
  TEST_FUNC(sf_bessel_j1);
  TEST_FUNC(sf_bessel_j2);
  TEST_FUNC_ND(sf_bessel_jl, GSL_NAN, l);
}

TEST_F(GSLTest, Bessely) {
  TEST_FUNC(sf_bessel_y0);
  TEST_FUNC(sf_bessel_y1);
  TEST_FUNC(sf_bessel_y2);
  TEST_FUNC_ND(sf_bessel_yl, GSL_NAN, l);
}

TEST_F(GSLTest, Besseli) {
  TEST_FUNC(sf_bessel_i0_scaled);
  TEST_FUNC(sf_bessel_i1_scaled);
  TEST_FUNC(sf_bessel_i2_scaled);
  TEST_FUNC_ND(sf_bessel_il_scaled, GSL_NAN, l);
}

TEST_F(GSLTest, Besselk) {
  TEST_FUNC(sf_bessel_k0_scaled);
  TEST_FUNC(sf_bessel_k1_scaled);
  TEST_FUNC(sf_bessel_k2_scaled);
  TEST_FUNC_ND(sf_bessel_kl_scaled, GSL_NAN, l);
}

struct BesselFractionalOrderInfo : FunctionInfo {
  Result GetDerivative(const Function &f,
      unsigned var_index, const Tuple &args) {
    // Partial derivative with respect to nu is not provided.
    if (var_index == 0)
      return Result("argument 'nu' is not constant");
    // Computing gsl_sf_bessel_*nu'(nu, x) requires
    // gsl_sf_bessel_*nu(nu - 1, x) which doesn't work when the
    // first argument is non-negative, so nu should be >= 1.
    if (args[0] < 1)
      return EvalError(f, args, "'");
    return Result();
  }

  Result GetSecondDerivative(
      const Function &f, unsigned, unsigned, const Tuple &args) {
    // Computing gsl_sf_bessel_*nu''(nu, x) requires
    // gsl_sf_bessel_*nu(nu - 2, x) which doesn't work when the
    // first argument is non-negative, so nu should be >= 2.
    return args[0] < 2 ? EvalError(f, args, "''") : Result();
  }
};

TEST_F(GSLTest, BesselFractionalOrder) {
  BesselFractionalOrderInfo info;
  TEST_FUNC(sf_bessel_Jnu);
  TEST_FUNC(sf_bessel_Ynu);
  TEST_FUNC(sf_bessel_Inu);
  TEST_FUNC(sf_bessel_Inu_scaled);
  TEST_FUNC(sf_bessel_Knu);
  TEST_FUNC(sf_bessel_lnKnu);
  TEST_FUNC(sf_bessel_Knu_scaled);
}

TEST_F(GSLTest, BesselZero) {
  TEST_FUNC(sf_bessel_zero_J0);
  TEST_FUNC(sf_bessel_zero_J1);

  const char *name = "gsl_sf_bessel_zero_Jnu";
  Function af = GetFunction(name);
  for (size_t i = 0; i != NUM_POINTS; ++i) {
    for (size_t j = 0; j != NUM_POINTS; ++j) {
      double nu = POINTS[i], s = POINTS[j];
      TestZeroFunc(af, gsl_sf_bessel_zero_Jnu(nu, s), Tuple(nu, s), 1);
    }
  }
}

struct ClausenFunctionInfo : FunctionInfo {
  Result GetDerivative(const Function &, unsigned, const Tuple &args) {
    return Result(args[0] == 0 ? GSL_POSINF : GSL_NAN);
  }
};

TEST_F(GSLTest, Clausen) {
  ClausenFunctionInfo info;
  TEST_FUNC(sf_clausen);
}

TEST_F(GSLTest, Hydrogenic) {
  TEST_FUNC(sf_hydrogenicR_1);

  Function f = GetFunction("gsl_sf_hydrogenicR");
  EXPECT_EQ(2, f(Tuple(1, 0, 1, 0)));
  EXPECT_ERROR("argument 'n' can't be represented as int, n = 1.1",
      f(Tuple(1.1, 0, 1, 0)));
  EXPECT_ERROR("argument 'l' can't be represented as int, l = 0.1",
      f(Tuple(1, 0.1, 1, 0)));
  for (size_t in = 0; in != NUM_POINTS_FOR_N; ++in) {
    int n = POINTS_FOR_N[in];
    for (size_t il = 0; il != NUM_POINTS_FOR_N; ++il) {
      int el = POINTS_FOR_N[il];
      for (size_t iz = 0; iz != NUM_POINTS; ++iz) {
        for (size_t ir = 0; ir != NUM_POINTS; ++ir) {
          double z = POINTS[iz], r = POINTS[ir];
          Tuple args(n, el, z, r);
          CheckFunction(gsl_sf_hydrogenicR(n, el, z, r), f, args);
          const char *error = "argument 'n' is not constant";
          EXPECT_ERROR(error, f(args, DERIVS));
          EXPECT_ERROR(error, f(args, HES));
          error = "argument 'l' is not constant";
          EXPECT_ERROR(error, f(args, DERIVS, BitSet("0111")));
          EXPECT_ERROR(error, f(args, HES, BitSet("0111")));
          error = "derivatives are not provided";
          EXPECT_ERROR(error, f(args, DERIVS, BitSet("0011")));
          EXPECT_ERROR(error, f(args, HES, BitSet("0011")));
        }
      }
    }
  }
}

TEST_F(GSLTest, Coulomb) {
  Function f = GetFunction("gsl_sf_coulomb_CL");
  for (size_t i = 0; i != NUM_POINTS; ++i) {
    for (size_t j = 0; j != NUM_POINTS; ++j) {
      double x = POINTS[i], y = POINTS[j];
      Tuple args(x, y);
      gsl_sf_result result = {};
      double value = gsl_sf_coulomb_CL_e(x, y, &result) ? GSL_NAN : result.val;
      CheckFunction(value, f, args);
      const char *error = "derivatives are not provided";
      EXPECT_ERROR(error, f(args, DERIVS));
      EXPECT_ERROR(error, f(args, HES));
    }
  }
}

TEST_F(GSLTest, Coupling3j) {
  double value = gsl_sf_coupling_3j(8, 20, 12, -2, 12, -10);
  EXPECT_NEAR(0.0812695955, value, 1e-5);
  Function f = GetFunction("gsl_sf_coupling_3j");
  Tuple args(8, 20, 12, -2, 12, -10);
  EXPECT_EQ(value, f(args));
  f(Tuple(0, 0, 0, 0, 0, 0));
  EXPECT_ERROR(NotIntError("two_ja"), f(Tuple(0.5, 0, 0, 0, 0, 0)));
  EXPECT_ERROR(NotIntError("two_jb"), f(Tuple(0, 0.5, 0, 0, 0, 0)));
  EXPECT_ERROR(NotIntError("two_jc"), f(Tuple(0, 0, 0.5, 0, 0, 0)));
  EXPECT_ERROR(NotIntError("two_ma"), f(Tuple(0, 0, 0, 0.5, 0, 0)));
  EXPECT_ERROR(NotIntError("two_mb"), f(Tuple(0, 0, 0, 0, 0.5, 0)));
  EXPECT_ERROR(NotIntError("two_mc"), f(Tuple(0, 0, 0, 0, 0, 0.5)));
  const char *error = "argument 'two_ja' is not constant";
  EXPECT_ERROR(error, f(args, DERIVS));
  EXPECT_ERROR(error, f(args, HES));
}

TEST_F(GSLTest, Coupling6j) {
  double value = gsl_sf_coupling_6j(2, 4, 6, 8, 10, 12);
  EXPECT_NEAR(0.0176295295, value, 1e-7);
  Function f = GetFunction("gsl_sf_coupling_6j");
  Tuple args(2, 4, 6, 8, 10, 12);
  EXPECT_EQ(value, f(args));
  EXPECT_TRUE(f(Tuple(0, 0, 0, 0, 0, 0)).error() == nullptr);
  EXPECT_ERROR(NotIntError("two_ja"), f(Tuple(0.5, 0, 0, 0, 0, 0)));
  EXPECT_ERROR(NotIntError("two_jb"), f(Tuple(0, 0.5, 0, 0, 0, 0)));
  EXPECT_ERROR(NotIntError("two_jc"), f(Tuple(0, 0, 0.5, 0, 0, 0)));
  EXPECT_ERROR(NotIntError("two_jd"), f(Tuple(0, 0, 0, 0.5, 0, 0)));
  EXPECT_ERROR(NotIntError("two_je"), f(Tuple(0, 0, 0, 0, 0.5, 0)));
  EXPECT_ERROR(NotIntError("two_jf"), f(Tuple(0, 0, 0, 0, 0, 0.5)));
  const char *error = "argument 'two_ja' is not constant";
  EXPECT_ERROR(error, f(args, DERIVS));
  EXPECT_ERROR(error, f(args, HES));
}

TEST_F(GSLTest, Coupling9j) {
  double value = gsl_sf_coupling_9j(6, 16, 18, 8, 20, 14, 12, 10, 4);
  EXPECT_NEAR(-0.000775648399, value, 1e-9);
  Function f = GetFunction("gsl_sf_coupling_9j");
  Tuple args(6, 16, 18, 8, 20, 14, 12, 10, 4);
  EXPECT_EQ(value, f(args));
  EXPECT_TRUE(f(Tuple(0, 0, 0, 0, 0, 0, 0, 0, 0)).error() == nullptr);
  EXPECT_ERROR(NotIntError("two_ja"), f(Tuple(0.5, 0, 0, 0, 0, 0, 0, 0, 0)));
  EXPECT_ERROR(NotIntError("two_jb"), f(Tuple(0, 0.5, 0, 0, 0, 0, 0, 0, 0)));
  EXPECT_ERROR(NotIntError("two_jc"), f(Tuple(0, 0, 0.5, 0, 0, 0, 0, 0, 0)));
  EXPECT_ERROR(NotIntError("two_jd"), f(Tuple(0, 0, 0, 0.5, 0, 0, 0, 0, 0)));
  EXPECT_ERROR(NotIntError("two_je"), f(Tuple(0, 0, 0, 0, 0.5, 0, 0, 0, 0)));
  EXPECT_ERROR(NotIntError("two_jf"), f(Tuple(0, 0, 0, 0, 0, 0.5, 0, 0, 0)));
  EXPECT_ERROR(NotIntError("two_jg"), f(Tuple(0, 0, 0, 0, 0, 0, 0.5, 0, 0)));
  EXPECT_ERROR(NotIntError("two_jh"), f(Tuple(0, 0, 0, 0, 0, 0, 0, 0.5, 0)));
  EXPECT_ERROR(NotIntError("two_ji"), f(Tuple(0, 0, 0, 0, 0, 0, 0, 0, 0.5)));
  const char *error = "argument 'two_ja' is not constant";
  EXPECT_ERROR(error, f(args, DERIVS));
  EXPECT_ERROR(error, f(args, HES));
}

TEST_F(GSLTest, Dawson) {
  TEST_FUNC(sf_dawson);
}

TEST_F(GSLTest, Debye) {
  TEST_FUNC(sf_debye_1);
  TEST_FUNC(sf_debye_2);
  TEST_FUNC(sf_debye_3);
  TEST_FUNC(sf_debye_4);
  TEST_FUNC(sf_debye_5);
  TEST_FUNC(sf_debye_6);
}

struct DilogFunctionInfo : FunctionInfo {
  Result GetDerivative(const Function &, unsigned, const Tuple &args) {
    return Result(args[0] == 1 ? GSL_POSINF : GSL_NAN);
  }
};

TEST_F(GSLTest, Dilog) {
  DilogFunctionInfo info;
  TEST_FUNC(sf_dilog);
}

struct NoDerivativeInfo : FunctionInfo {
  Result GetDerivative(const Function &, unsigned, const Tuple &) {
    return Result("derivatives are not provided");
  }
};

TEST_F(GSLTest, EllInt) {
  TEST_FUNC(sf_ellint_Kcomp);
  TEST_FUNC(sf_ellint_Ecomp);
  TEST_FUNC(sf_ellint_Pcomp);
  TEST_FUNC(sf_ellint_F);
  TEST_FUNC(sf_ellint_E);
  {
    NoDerivativeInfo info;
    TEST_FUNC(sf_ellint_P);
    TEST_FUNC(sf_ellint_D);
    TEST_FUNC(sf_ellint_RC);
    TEST_FUNC(sf_ellint_RD);
    TEST_FUNC(sf_ellint_RF);
  }
  Function f = GetFunction("gsl_sf_ellint_RJ");
  for (size_t ix = 0; ix != NUM_POINTS; ++ix) {
    for (size_t iy = 0; iy != NUM_POINTS; ++iy) {
      for (size_t iz = 0; iz != NUM_POINTS; ++iz) {
        for (size_t ip = 0; ip != NUM_POINTS; ++ip) {
          double x = POINTS[ix], y = POINTS[iy];
          double z = POINTS[iz], p = POINTS[ip];
          Tuple args(x, y, z, p);
          CheckFunction(gsl_sf_ellint_RJ(x, y, z, p, GSL_PREC_DOUBLE), f, args);
          const char *error = "derivatives are not provided";
          EXPECT_ERROR(error, f(args, DERIVS));
          EXPECT_ERROR(error, f(args, HES));
        }
      }
    }
  }
}

TEST_F(GSLTest, Erf) {
  TEST_FUNC(sf_erf);
  TEST_FUNC(sf_erfc);
  TEST_FUNC(sf_log_erfc);
  TEST_FUNC(sf_erf_Z);
  TEST_FUNC(sf_erf_Q);
  TEST_FUNC(sf_hazard);
}

TEST_F(GSLTest, ExpInt) {
  TEST_FUNC(sf_expint_E1);
  TEST_FUNC(sf_expint_E2);
  TEST_FUNC_ND(sf_expint_En, GSL_NAN, n);
  TEST_FUNC(sf_expint_Ei);
  TEST_FUNC(sf_Shi);
  TEST_FUNC(sf_Chi);
  TEST_FUNC(sf_expint_3);
  TEST_FUNC(sf_Si);
  TEST_FUNC(sf_Ci);
  TEST_FUNC(sf_atanint);
}

TEST_F(GSLTest, FermiDirac) {
  TEST_FUNC(sf_fermi_dirac_m1);
  TEST_FUNC(sf_fermi_dirac_0);
  TEST_FUNC(sf_fermi_dirac_1);
  TEST_FUNC(sf_fermi_dirac_2);
  TEST_FUNC_ND(sf_fermi_dirac_int, GSL_NAN, j);
  {
    NoDerivativeInfo info;
    TEST_FUNC(sf_fermi_dirac_mhalf);
    TEST_FUNC(sf_fermi_dirac_half);
  }
  TEST_FUNC(sf_fermi_dirac_3half);
  TEST_FUNC(sf_fermi_dirac_inc_0);
}

struct LnGammaInfo : FunctionInfo {
  Result GetDerivative(const Function &af, unsigned , const Tuple &args) {
    double x = args[0];
    return x == -1 || x == -2 ? EvalError(af, args, "'") : Result();
  }
};

TEST_F(GSLTest, Gamma) {
  Function gamma = GetFunction("gsl_sf_gamma", &info);
  EXPECT_NEAR(-0.129354, gamma(Tuple(-0.5), DERIVS).deriv(), 1e-6);
  EXPECT_NEAR(-31.6778, gamma(Tuple(-0.5), HES).hes(), 1e-4);
  EXPECT_NEAR(1.19786e100, gamma(Tuple(71)), 1e95);
  EXPECT_TRUE(gsl_isinf(gamma(Tuple(1000))));
  TEST_FUNC(sf_gamma);
  {
    LnGammaInfo info;
    TEST_FUNC(sf_lngamma);
  }
  TEST_FUNC(sf_gammastar);
  TEST_FUNC(sf_gammainv);
}

TEST_F(GSLTest, Poch) {
  NoDerivativeInfo info;
  TEST_FUNC(sf_poch);
  TEST_FUNC(sf_lnpoch);
  TEST_FUNC(sf_pochrel);
}

struct GammaIncInfo : FunctionInfo {
  Result GetDerivative(
      const Function &af, unsigned var_index, const Tuple &args) {
    // Partial derivative with respect to a is not provided.
    if (var_index == 0)
      return Result("argument 'a' is not constant");
    if (args[1] == 0)
      return EvalError(af, args, "'");
    return Result();
  }
};

TEST_F(GSLTest, GammaInc) {
  {
    GammaIncInfo info;
    TEST_FUNC(sf_gamma_inc);
  }
  {
    NoDerivativeInfo info;
    TEST_FUNC(sf_gamma_inc_Q);
    TEST_FUNC(sf_gamma_inc_P);
  }
}

struct BetaInfo : FunctionInfo {
  Result GetDerivative(
      const Function &af, unsigned var_index, const Tuple &args) {
    if (gsl_isnan(gsl_sf_psi(args[0] + args[1])) || args[var_index] == 0)
      return EvalError(af, args, "'");
    return Result();
  }
};

TEST_F(GSLTest, Beta) {
  {
    BetaInfo info;
    TEST_FUNC(sf_beta);
  }
  TEST_FUNC(sf_lnbeta);
  {
    NoDerivativeInfo info;
    TEST_FUNC(sf_beta_inc);
  }
}

TEST_F(GSLTest, GegenPoly) {
  TEST_FUNC(sf_gegenpoly_1);
  TEST_FUNC(sf_gegenpoly_2);
  TEST_FUNC(sf_gegenpoly_3);

  NoDerivativeInfo info;
  info.SetArgNames("n");
  TEST_FUNC(sf_gegenpoly_n);
}

struct Hyperg0F1Info : FunctionInfo {
  Result GetDerivative(const Function &, unsigned var_index, const Tuple &) {
    // Partial derivative with respect to c is not provided.
    return Result(var_index == 0 ? "argument 'c' is not constant" : "");
  }
};

struct Hyperg1F1Info : FunctionInfo {
  Result GetDerivative(const Function &f, unsigned, const Tuple &args) {
    return args[1] <= 0 ? EvalError(f, args, "'") : Result();
  }
};

TEST_F(GSLTest, Hyperg) {
  {
    Hyperg0F1Info info;
    TEST_FUNC(sf_hyperg_0F1);
  }
  {
    Hyperg1F1Info info;
    info.SetArgNames("m n x");
    TEST_FUNC(sf_hyperg_1F1_int);
  }
}
}
