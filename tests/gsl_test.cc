// GSL wrapper test.

#include <functional>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>

#include "gtest/gtest.h"
#include "solvers/asl.h"
#include "tests/config.h"
//#define DEBUG_DIFFERENTIATOR
#include "tests/functional.h"

using std::string;
using std::vector;
using fun::Tuple;

namespace {

// An immutable result of an AMPL function call.
class Result {
 private:
  real value_;
  vector<real> derivs_;
  vector<real> hes_;
  const char *error_;

  void CheckError() const {
    if (error_)
      throw std::runtime_error(error_);
  }

 public:
  Result(real value, const vector<real> &derivs,
      const vector<real> &hes, const char *error) :
    value_(value), derivs_(derivs), hes_(hes), error_(error) {}

  operator real() const {
    CheckError();
    return value_;
  }

  real deriv(size_t index = 0) const {
    CheckError();
    return derivs_.at(index);
  }
  real hes(size_t index = 0) const {
    CheckError();
    return hes_.at(index);
  }

  const char *error() const { return error_; }
};

class Function;

// Function information that can't be obtained automatically, in particular
// due to limitations of numerical differentiation.
class FunctionInfo {
 private:
  vector<string> arg_names_;

 public:
  virtual ~FunctionInfo() {}

  string arg_name(unsigned index) const {
    return index < arg_names_.size() ? arg_names_[index] : string();
  }

  void set_arg_names(const char *arg_names) {
    std::istringstream is(arg_names);
    copy(std::istream_iterator<string>(is), std::istream_iterator<string>(),
        std::back_inserter<vector<string>>(arg_names_));
  }

  virtual double GetDerivative(unsigned, const Tuple &) { return GSL_NAN; }
  virtual double GetSecondDerivative(unsigned, unsigned, const Tuple &) {
    return GSL_NAN;
  }

  virtual string DerivativeError(const Function &, unsigned, const Tuple &) {
    return "";
  }
  virtual string Derivative2Error(const Function &, const Tuple &) {
    return "";
  }
};

// Flags for an AMPL function call.
enum {
  DERIVS = 1, // Get first partial derivatives.
  HES    = 3  // Get both first and second partial derivatives.
};

// An AMPL function.
class Function {
 private:
  ASL *asl_;
  func_info *fi_;
  FunctionInfo *info_;

 public:
  Function(ASL *asl, func_info *fi, FunctionInfo *info) :
    asl_(asl), fi_(fi), info_(info) {}

  const char *name() const { return fi_->name; }

  FunctionInfo *info() const { return info_; }

  // Calls a function.
  // Argument vector is passed by value intentionally to avoid
  // rogue functions accidentally overwriting arguments.
  Result operator()(const Tuple &args,
      int flags = 0, char *dig = 0, void *info = 0) const;

  string DerivativeError(unsigned var_index, const Tuple &args) const {
    return info_->DerivativeError(*this, var_index, args);
  }
  string Derivative2Error(const Tuple &args) const {
    return info_->Derivative2Error(*this, args);
  }

  string arg_name(unsigned index) const {
    return info_->arg_name(index);
  }
};

Result Function::operator()(
    const Tuple &args, int flags, char *dig, void *info) const {
  if (fi_->nargs != static_cast<int>(args.size()))
    throw std::runtime_error("Invalid number of arguments in function call");

  // Initialize the argument list.
  vector<real> ra(args.size());
  for (unsigned i = 0; i < args.size(); ++i)
    ra[i] = args[i];
  arglist al = {};
  TMInfo tmi = {};
  al.ra = &ra[0];
  al.nr = al.n = args.size();
  al.TMI = &tmi;
  al.AE = asl_->i.ae;
  al.dig = dig;
  al.funcinfo = info;

  // Allocate storage for the derivatives if needed.
  vector<real> derivs, hes;
  if ((flags & DERIVS) != 0) {
    derivs.resize(al.n);
    al.derivs = &derivs[0];
  }
  if ((flags & HES) == HES) {
    hes.resize(al.n * (al.n + 1) / 2);
    al.hes = &hes[0];
  }

  // Call the function and return the result.
  real value = fi_->funcp(&al);
  return Result(value, derivs, hes, al.Errmsg);
}

// Converts error estimate returned by Diff into an absolute tolerance to
// be used in EXPECT_NEAR.
double ConvertErrorToTolerance(double error) {
  return error != 0 ? error * 1000 : 1e-10;
}

class Error {
 private:
  string str;

 public:
  Error(const string &s) : str(s) {}
  operator const char*() const { return str.c_str(); }
  const char* c_str() const { return str.c_str(); }
};

Error EvalError(const Function &f, const Tuple &args, const char *suffix = "") {
  std::ostringstream os;
  os << "can't evaluate " << f.name() << suffix << args;
  return os.str();
}

Error NotIntError(const string &arg_name, double value = 0.5) {
  std::ostringstream os;
  os << "argument '" << arg_name
      << "' can't be represented as int, " << arg_name << " = " << value;
  return os.str();
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

// A helper class that converts a variable index into the dig array.
// See the dig member of the arglist struct for details.
class Dig {
 private:
  char *dig_;
  vector<char> store_;

 public:
  static const unsigned NO_VAR = ~0u;

  Dig(const Tuple &args, unsigned skip_var) : dig_(0) {
    if (skip_var == NO_VAR) return;
    store_.resize(args.size());
    store_.at(skip_var) = 1;
    dig_ = &store_[0];
  }

  operator char*() { return dig_; }
};

// A helper class that wraps a Function's derivative and binds
// all but one argument to the given values.
class DerivativeBinder {
 private:
  Function f_;
  unsigned deriv_var_;
  unsigned eval_var_;
  Tuple args_;
  vector<char> dig_;

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
  args_(args), dig_(args.size(), 1) {
  unsigned num_vars = args_.size();
  if (deriv_var >= num_vars || eval_var >= num_vars)
    throw std::out_of_range("variable index is out of range");
  dig_[deriv_var] = 0;
}

double DerivativeBinder::operator()(double x) {
  args_[eval_var_] = x;
  Result r = f_(args_, DERIVS, &dig_[0]);
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
  FunctionInfo info; // Default function info.

  fun::Differentiator diff;

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
      throw std::runtime_error(string("Function not found: ") + name);
    return Function(asl, fi, info);
  }

  template <typename F>
  bool CheckDerivative(F f, const Function &af,
      unsigned var_index, const Tuple &args);

  void CheckSecondDerivatives(const Function &f,
      const Tuple &args, unsigned skip_var = Dig::NO_VAR);

  // Tests a function taking a single argument.
  template <typename F>
  void TestUnaryFunc(const Function &af, F f);
  void TestFunc(const Function &af, double (*f)(double)) {
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
      TestZeroFunc(af, f(s), s, 0);
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
    Func2DoubleMode(Func2Mode f) : f_(f) {}

    double operator()(double x, double y) const {
      return f_(x, y, GSL_PREC_DOUBLE);
    }
  };

  void TestFunc(const Function &af, Func2Mode f) {
    TestBinaryFunc(af, Func2DoubleMode(f));
  }

  // Binds the mode argument of Func2Mode to GSL_PREC_DOUBLE.
  class Func3DoubleMode :
    public fun::ternary_function<double, double, double, double> {
   private:
    Func3Mode f_;

   public:
    Func3DoubleMode(Func3Mode f) : f_(f) {}

    double operator()(double x, double y, double z) const {
      return f_(x, y, z, GSL_PREC_DOUBLE);
    }
  };

  template <typename F>
  void TestTernaryFunc(const Function &af, F f);

  template <typename Arg1, typename Arg2, typename Arg3, typename Result>
  void TestFunc(const Function &af, Result (*f)(Arg1, Arg2, Arg3)) {
    TestTernaryFunc(af,
        fun::pointer_to_ternary_function<Arg1, Arg2, Arg3, Result>(f));
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
  vector<char> dig(args.size(), 1);
  dig[var_index] = 0;
  double error = 0;
  double x = args[var_index];
  string error_message = af.DerivativeError(var_index, args);
  if (error_message.empty()) {
    double numerical_deriv = Diff(f, x, &error);
    double overridden_deriv = af.info()->GetDerivative(var_index, args);
    if (!gsl_isnan(overridden_deriv) && overridden_deriv != numerical_deriv) {
      std::cout << "Overriding d/dx" << var_index << " " << af.name()
        << " at " << args << ", computed = " << numerical_deriv
        << ", overridden = " << overridden_deriv << std::endl;
      numerical_deriv = overridden_deriv;
    }
    if (!gsl_isnan(numerical_deriv)) {
      double deriv = af(args, DERIVS, &dig[0]).deriv(var_index);
      if (numerical_deriv != deriv)
        EXPECT_NEAR(numerical_deriv, deriv, ConvertErrorToTolerance(error));
      return true;
    }
  }
  Result r = af(args, DERIVS, &dig[0]);
  if (!gsl_isnan(f(x))){
    if (error_message.empty())
      error_message = EvalError(af, args, "'");
    EXPECT_ERROR(error_message.c_str(), r);
  } else
    EXPECT_TRUE(r.error() != nullptr);
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
  if (skip_var == Dig::NO_VAR) {
    for (unsigned i = 0; i < num_args; ++i) {
      if (!f.DerivativeError(i, args).empty()) {
        skip_var = i;
        break;
      }
    }
  }
  for (unsigned i = 0; i < num_args; ++i) {
    if (i == skip_var) continue;
    for (unsigned j = 0; j < num_args; ++j) {
      if (j == skip_var) continue;
      vector<char> dig(num_args, 1);
      dig[i] = 0;
      dig[j] = 0;
      double error = 0;
      string error_message = f.Derivative2Error(args);
      if (error_message.empty()) {
        double d = Diff(DerivativeBinder(f, j, i, args), args[i], &error);
        double overridden_deriv = f.info()->GetSecondDerivative(j, i, args);
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
          double actual_deriv = f(args, HES, &dig[0]).hes(hes_index);
          if (d != actual_deriv)
            EXPECT_NEAR(d, actual_deriv, ConvertErrorToTolerance(error));
          continue;
        }
      }
      Result r = f(args, HES, &dig[0]);
      if (f(args, DERIVS, &dig[0]).error())
        EXPECT_TRUE(r.error() != nullptr);
      else if (!error_message.empty())
        EXPECT_ERROR(error_message.c_str(), r);
      else
        EXPECT_ERROR(EvalError(f, args, "''"), r);
    }
  }
}

template <typename F>
void GSLTest::TestUnaryFunc(const Function &af, F f) {
  for (size_t i = 0; i != NUM_POINTS; ++i) {
    double x = POINTS[i];
    CheckFunction(f(x), af, x);
    CheckDerivative(f, af, 0, x);
    CheckSecondDerivatives(af, x);
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

  char dig[2] = {1, 0};
  EXPECT_ERROR(
      ("can't compute derivative: argument '" + arg_name + "' too small, " +
      arg_name + " = -2147483648").c_str(),
      af(Tuple(INT_MIN, test_x), DERIVS, dig));
  EXPECT_ERROR(
      ("can't compute derivative: argument '" + arg_name + "' too large, " +
      arg_name + " = 2147483647").c_str(),
      af(Tuple(INT_MAX, test_x), DERIVS, dig));
  EXPECT_TRUE(!gsl_isnan(af(Tuple(INT_MIN + 1, test_x), DERIVS, dig).deriv(1)));
  EXPECT_TRUE(!gsl_isnan(af(Tuple(INT_MAX - 1, test_x), DERIVS, dig).deriv(1)));

  EXPECT_ERROR(
      ("can't compute derivative: argument '" + arg_name + "' too small, " +
      arg_name + " = -2147483647").c_str(),
      af(Tuple(INT_MIN + 1, test_x), HES, dig));
  EXPECT_ERROR(
      ("can't compute derivative: argument '" + arg_name + "' too large, " +
      arg_name + " = 2147483646").c_str(),
      af(Tuple(INT_MAX - 1, test_x), HES, dig));
  EXPECT_TRUE(!gsl_isnan(af(Tuple(INT_MIN + 2, test_x), HES, dig).hes(2)));
  EXPECT_TRUE(!gsl_isnan(af(Tuple(INT_MAX - 2, test_x), HES, dig).hes(2)));
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
  vector<char> dig(args.size(), 1);
  dig.at(arg_index) = 0;
  os << "argument '" << f.arg_name(arg_index) << "' is not constant";
  EXPECT_STREQ(os.str().c_str(), f(args, DERIVS, &dig[0]).error());
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
          EXPECT_STREQ(NotIntError(af.arg_name(0), x), af(args).error());
          continue;
        }
        if (static_cast<typename F::second_argument_type>(y) != y) {
          EXPECT_STREQ(NotIntError(af.arg_name(1), y), af(args).error());
          continue;
        }
        CheckFunction(f(x, y, z), af, args);
        if (!CheckArg<typename F::first_argument_type>(af, args, 0))
          CheckDerivative(fun::Bind2Of3(f, args, 0), af, 0, args);
        if (!CheckArg<typename F::second_argument_type>(af, args, 1))
          CheckDerivative(fun::Bind2Of3(f, args, 1), af, 1, args);
        if (!CheckArg<typename F::third_argument_type>(af, args, 2))
          CheckDerivative(fun::Bind2Of3(f, args, 2), af, 2, args);
        CheckSecondDerivatives(af, args);
      }
    }
  }
}

#define TEST_FUNC(name) TestFunc(GetFunction("gsl_" #name, &info), gsl_##name)

#define TEST_FUNC_ND(name, test_x, arg) \
  TestFuncND(GetFunction("gsl_" #name, &info), gsl_##name, test_x, #arg)

TEST_F(GSLTest, Result) {
  static const real ARGS[] = {5, 7, 11, 13, 17};
  Result r(42, vector<real>(ARGS, ARGS + 2),
      vector<real>(ARGS + 2, ARGS + 5), nullptr);
  EXPECT_EQ(42, r);
  EXPECT_EQ(5, r.deriv());
  EXPECT_EQ(5, r.deriv(0));
  EXPECT_EQ(7, r.deriv(1));
  EXPECT_THROW(r.deriv(2), std::out_of_range);
  EXPECT_EQ(11, r.hes());
  EXPECT_EQ(11, r.hes(0));
  EXPECT_EQ(13, r.hes(1));
  EXPECT_EQ(17, r.hes(2));
  EXPECT_THROW(r.hes(3), std::out_of_range);
  EXPECT_TRUE(r.error() == nullptr);
}

TEST_F(GSLTest, ErrorResult) {
  static const real ARGS[] = {5, 7, 11, 13, 17};
  const char *error = "brain overflow";
  Result r(42, vector<real>(ARGS, ARGS + 2),
      vector<real>(ARGS + 2, ARGS + 5), error);
  EXPECT_THROW((double)r, std::runtime_error);
  EXPECT_THROW(r.deriv(), std::runtime_error);
  EXPECT_THROW(r.deriv(0), std::runtime_error);
  EXPECT_THROW(r.hes(), std::runtime_error);
  EXPECT_THROW(r.hes(0), std::runtime_error);
  EXPECT_ERROR(error, r);
}

struct CallData {
  AmplExports *ae;
  int n;
  int nr;
  vector<real> ra;
  real *derivs;
  real *hes;
  char *dig;
  char *error;
};

real Test(arglist *args) {
  CallData *data = reinterpret_cast<CallData*>(args->funcinfo);
  data->ae = args->AE;
  data->n = args->n;
  data->nr = args->nr;
  data->ra = vector<real>(args->ra, args->ra + args->n);
  data->derivs = args->derivs;
  data->hes = args->hes;
  data->dig = args->dig;
  data->error = args->Errmsg;
  if (args->ra[0] < 0)
    args->Errmsg = const_cast<char*>("oops");
  if (args->derivs) {
    args->derivs[0] = 123;
    args->derivs[1] = 456;
    if (args->n > 2)
      args->derivs[2] = 789;
  }
  if (args->hes) {
    args->hes[0] = 12;
    args->hes[1] = 34;
    args->hes[2] = 56;
  }
  return 42;
}

class TestFunction {
 private:
  ASL testASL_;
  AmplExports ae_;
  func_info fi_;
  Function f_;

 public:
  TestFunction(int nargs)
  : testASL_(), ae_(), fi_(), f_(&testASL_, &fi_, 0) {
    testASL_.i.ae = &ae_;
    fi_.nargs = nargs;
    fi_.funcp = Test;
  }

  const AmplExports* ae() const { return &ae_; }
  const Function& get() const { return f_; }
};

TEST_F(GSLTest, FunctionCall) {
  TestFunction f(1);
  CallData data = {};
  EXPECT_EQ(42, f.get()(777, 0, 0, &data));
  EXPECT_EQ(f.ae(), data.ae);
  ASSERT_EQ(1, data.n);
  EXPECT_EQ(1, data.nr);
  EXPECT_EQ(777, data.ra[0]);
  EXPECT_TRUE(data.derivs == nullptr);
  EXPECT_TRUE(data.hes == nullptr);
  EXPECT_TRUE(data.dig == nullptr);
  EXPECT_TRUE(data.error == nullptr);
}

TEST_F(GSLTest, FunctionReturnsError) {
  TestFunction f(1);
  CallData data = {};
  EXPECT_ERROR("oops", f.get()(-1, 0, 0, &data));
}

TEST_F(GSLTest, FunctionReturnsDerivs) {
  TestFunction f(3);
  CallData data = {};
  Result res = f.get()(Tuple(11, 22, 33), DERIVS, 0, &data);
  EXPECT_EQ(42, res);
  EXPECT_EQ(f.ae(), data.ae);
  ASSERT_EQ(3, data.n);
  EXPECT_EQ(3, data.nr);
  EXPECT_EQ(11, data.ra[0]);
  EXPECT_EQ(22, data.ra[1]);
  EXPECT_EQ(33, data.ra[2]);
  EXPECT_EQ(123, res.deriv(0));
  EXPECT_EQ(456, res.deriv(1));
  EXPECT_EQ(789, res.deriv(2));
  EXPECT_THROW(res.deriv(3), std::out_of_range);
  EXPECT_TRUE(data.hes == nullptr);
  EXPECT_TRUE(data.dig == nullptr);
  EXPECT_TRUE(data.error == nullptr);
}

TEST_F(GSLTest, FunctionReturnsHes) {
  TestFunction f(2);
  CallData data = {};
  Result res = f.get()(Tuple(111, 222), HES, 0, &data);
  EXPECT_EQ(42, res);
  EXPECT_EQ(f.ae(), data.ae);
  ASSERT_EQ(2, data.n);
  EXPECT_EQ(2, data.nr);
  EXPECT_EQ(111, data.ra[0]);
  EXPECT_EQ(222, data.ra[1]);
  EXPECT_EQ(123, res.deriv(0));
  EXPECT_EQ(456, res.deriv(1));
  EXPECT_THROW(res.deriv(2), std::out_of_range);
  EXPECT_EQ(12, res.hes(0));
  EXPECT_EQ(34, res.hes(1));
  EXPECT_EQ(56, res.hes(2));
  EXPECT_THROW(res.hes(3), std::out_of_range);
  EXPECT_TRUE(data.dig == nullptr);
  EXPECT_TRUE(data.error == nullptr);
}

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

TEST_F(GSLTest, Dig) {
  EXPECT_TRUE(Dig(0, Dig::NO_VAR) == nullptr);
  EXPECT_EQ(1, Dig(0, 0)[0]);
  EXPECT_THROW(Dig(0, 1), std::out_of_range);
  EXPECT_TRUE(Dig(Tuple(0, 0, 0), Dig::NO_VAR) == nullptr);
  EXPECT_EQ(string("\1\0\0", 3), string(Dig(Tuple(0, 0, 0), 0), 3));
  EXPECT_EQ(string("\0\1\0", 3), string(Dig(Tuple(0, 0, 0), 1), 3));
  EXPECT_EQ(string("\0\0\1", 3), string(Dig(Tuple(0, 0, 0), 2), 3));
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
  string DerivativeError(const Function &f,
      unsigned var_index, const Tuple &args) {
    // Partial derivative with respect to nu is not provided.
    if (var_index == 0)
      return "argument 'nu' is not constant";
    // Computing gsl_sf_bessel_*nu'(nu, x) requires
    // gsl_sf_bessel_*nu(nu - 1, x) which doesn't work when the
    // first argument is non-negative, so nu should be >= 1.
    if (args[0] < 1)
      return EvalError(f, args, "'").c_str();
    return "";
  }

  string Derivative2Error(const Function &f, const Tuple &args) {
    // Computing gsl_sf_bessel_*nu''(nu, x) requires
    // gsl_sf_bessel_*nu(nu - 2, x) which doesn't work when the
    // first argument is non-negative, so nu should be >= 2.
    return args[0] < 2 ? EvalError(f, args, "''").c_str() : "";
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
  double GetDerivative(unsigned, const Tuple &args) {
    return args[0] == 0 ? GSL_POSINF : GSL_NAN;
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
          char dig1[] = {1, 0, 0, 0};
          error = "argument 'l' is not constant";
          EXPECT_ERROR(error, f(args, DERIVS, dig1));
          EXPECT_ERROR(error, f(args, HES, dig1));
          error = "derivatives are not provided";
          char dig2[] = {1, 1, 0, 0};
          EXPECT_ERROR(error, f(args, DERIVS, dig2));
          EXPECT_ERROR(error, f(args, HES, dig2));
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
  double GetDerivative(unsigned, const Tuple &args) {
    return args[0] == 1 ? GSL_POSINF : GSL_NAN;
  }
};

TEST_F(GSLTest, Dilog) {
  DilogFunctionInfo info;
  TEST_FUNC(sf_dilog);
}

struct NoDerivativeInfo : FunctionInfo {
  string DerivativeError(const Function &, unsigned, const Tuple &) {
    return "derivatives are not provided";
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
  string DerivativeError(const Function &af,
        unsigned , const Tuple &args) {
    double x = args[0];
    return x == -1 || x == -2 ? EvalError(af, args, "'") : "";
  }
};

TEST_F(GSLTest, Gamma) {
  Function gamma = GetFunction("gsl_sf_gamma", &info);
  EXPECT_NEAR(-0.129354, gamma(-0.5, DERIVS).deriv(), 1e-6);
  EXPECT_NEAR(-31.6778, gamma(-0.5, HES).hes(), 1e-4);
  EXPECT_NEAR(1.19786e100, gamma(71), 1e95);
  EXPECT_TRUE(gsl_isinf(gamma(1000)));
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
  string DerivativeError(const Function &af,
      unsigned var_index, const Tuple &args) {
    // Partial derivative with respect to a is not provided.
    if (var_index == 0)
      return "argument 'a' is not constant";
    if (args[1] == 0)
      return EvalError(af, args, "'").c_str();
    return "";
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
  string DerivativeError(const Function &af,
      unsigned var_index, const Tuple &args) {
    if (gsl_isnan(gsl_sf_psi(args[0] + args[1])) || args[var_index] == 0)
      return EvalError(af, args, "'").c_str();
    return "";
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
  info.set_arg_names("n");
  TEST_FUNC(sf_gegenpoly_n);
}

struct Hyperg0F1Info : FunctionInfo {
  string DerivativeError(const Function &,
      unsigned var_index, const Tuple &) {
    // Partial derivative with respect to c is not provided.
    return var_index == 0 ? "argument 'c' is not constant" : "";
  }
};

struct Hyperg1F1Info : FunctionInfo {
  string DerivativeError(const Function &f, unsigned, const Tuple &args) {
    return args[1] <= 0 ? EvalError(f, args, "'").c_str() : "";
  }
};

TEST_F(GSLTest, Hyperg) {
  {
    Hyperg0F1Info info;
    TEST_FUNC(sf_hyperg_0F1);
  }
  {
    Hyperg1F1Info info;
    info.set_arg_names("m n x");
    TEST_FUNC(sf_hyperg_1F1_int);
  }
}
}
