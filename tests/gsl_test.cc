// GSL wrapper test.

#include <functional>
#include <stdexcept>
#include <vector>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>

#include "gtest/gtest.h"
#include "solvers/asl.h"
#include "tests/config.h"

using std::vector;

namespace {

// An immutable list of arguments for an AMPL function.
class ArgList {
 private:
  vector<real> ra;

  ArgList &operator<<(real arg) {
    ra.push_back(arg);
    return *this;
  }

 public:
  ArgList(real a0) { *this << a0; }
  ArgList(real a0, real a1) { *this << a0 << a1; }
  ArgList(real a0, real a1, real a2) { *this << a0 << a1 << a2; }
  ArgList(real a0, real a1, real a2, real a3) {
    *this << a0 << a1 << a2 << a3;
  }
  ArgList(real a0, real a1, real a2, real a3, real a4) {
    *this << a0 << a1 << a2 << a3 << a4;
  }
  ArgList(real a0, real a1, real a2, real a3, real a4, real a5) {
    *this << a0 << a1 << a2 << a3 << a4 << a5;
  }
  ArgList(real a0, real a1, real a2, real a3,
      real a4, real a5, real a6, real a7, real a8) {
    *this << a0 << a1 << a2 << a3 << a4 << a5 << a6 << a7 << a8;
  }

  const vector<real> &get() const { return ra; }

  friend std::ostream &operator<<(std::ostream &os, const ArgList &args);
};

std::ostream &operator<<(std::ostream &os, const ArgList &args) {
  os << "(";
  if (!args.ra.empty()) {
    os << args.ra.front();
    for (size_t i = 1, n = args.ra.size(); i < n; ++i)
      os << ", " << args.ra[i];
  }
  os << ")";
  return os;
}

// An immutable result of an AMPL function call.
class Result {
 private:
  real value_;
  vector<real> derivs_;
  vector<real> hes_;

 public:
  Result(real value, const vector<real> &derivs, const vector<real> &hes) :
    value_(value), derivs_(derivs), hes_(hes) {}

  operator real() const { return value_; }

  real deriv(size_t index = 0) const { return derivs_.at(index); }
  real hes(size_t index = 0) const { return hes_.at(index); }
};

// Options for an AMPL function call.
enum {
  ERROR        = 1, // Function call is expected to produce an error.
  DERIVS       = 2, // Get first partial derivatives.
  HES          = 6, // Get both first and second partial derivatives.
  ERROR_TO_NAN = 8  // Convert error to NaN
};

// An AMPL function.
class Function {
 private:
  ASL *asl_;
  func_info *fi_;

 public:
  Function(ASL *asl, func_info *fi) : asl_(asl), fi_(fi) {}

  const char *name() const { return fi_->name; }

  Result operator()(vector<real> args,
      int options = 0, char *dig = 0, void *info = 0) const;

  Result operator()(const ArgList &args,
      int options = 0, char *dig = 0, void *info = 0) const {
    return (*this)(args.get(), options, dig, info);
  }

  Result operator()(real x,
      int options = 0, char *dig = 0, void *info = 0) const {
    return (*this)((ArgList(x)), options, dig, info);
  }
};

Result Function::operator()(
    vector<real> args, int options, char *dig, void *info) const {
  // Initialize the argument list.
  arglist al = {};
  TMInfo tmi = {};
  al.ra = &args[0];
  al.nr = al.n = args.size();
  al.TMI = &tmi;
  al.AE = asl_->i.ae;
  al.dig = dig;
  al.funcinfo = info;

  // Allocate storage for the derivatives if needed.
  vector<real> derivs, hes;
  if ((options & DERIVS) != 0) {
    derivs.resize(args.size());
    al.derivs = &derivs[0];
  }
  if ((options & HES) == HES) {
    hes.resize(args.size() * (args.size() + 1) / 2);
    al.hes = &hes[0];
  }

  // Call the function.
  real value = fi_->funcp(&al);

  // Check the error message.
  if (al.Errmsg) {
    if ((options & ERROR_TO_NAN) != 0)
      value = GSL_NAN;
    else if ((options & ERROR) == 0)
      ADD_FAILURE() << al.Errmsg;
  } else if ((options & ERROR) != 0)
    ADD_FAILURE() << "Expected error in " << fi_->name;

  return Result(value, derivs, hes);
}

template <typename F>
inline double SymmetricDifference(F f, double x, double h) {
  return (f(x + h) - f(x - h)) / (2 * h);
}

template <typename F>
inline double LeftDifference(F f, double x, double h) {
  return (f(x) - f(x - h)) / h;
}

template <typename F>
inline double RightDifference(F f, double x, double h) {
  return (f(x + h) - f(x)) / h;
}

// Returns the derivative of a function f at a point x by Ridders'
// method of polynomial extrapolation. The implementation is taken from
// "Numerical Recipes in C", Chapter 5.7.
template <typename F, typename D>
double Diff(F f, double x, D d, double *error = 0, double *final_h = 0) {
  const int NTAB = 200;
  const double CON = 1.4, CON2 = CON * CON;
  const double BIG = std::numeric_limits<double>::max();
  const double SAFE = 2;
  double h = 0.125, ans = GSL_NAN;

  // Successive columns in the Neville tableau will go to smaller
  // step sizes and higher orders of extrapolation.
  vector<double> a(NTAB * NTAB);

  a[0] = d(f, x, h);
  double err = BIG;
  for (int i = 1; i < NTAB; i++) {
    // Try new, smaller step size.
    h /= CON;
    a[i] = d(f, x, h);
    double fac = CON2;
    for (int j = 1; j <= i; j++) {
      // Compute extrapolations of various orders, requiring no new function
      // evaluations.
      a[j * NTAB + i] =
          (a[(j - 1) * NTAB + i] * fac - a[(j - 1) * NTAB + i - 1]) /
          (fac - 1);
      fac = CON2 * fac;
      double errt = std::max(
          fabs(a[j * NTAB + i] - a[(j - 1) * NTAB + i]),
          fabs(a[j * NTAB + i] - a[(j - 1) * NTAB + i - 1]));
      // The error strategy is to compare each new extrapolation to one order
      // lower, both at the present step size and the previous one.
      if (errt <= err) {
        // If error is decreased, save the improved answer.
        err = errt;
        ans = a[j * NTAB + i];
      }
    }
    // If higher order is worse by a significant factor SAFE, then quit early.
    if (fabs(a[i * NTAB + i] - a[(i - 1) * NTAB + i - 1]) >= SAFE * err)
      break;
  }

  if (error)
    *error = err;
  if (final_h)
    *final_h = h;
  return ans;
}

// Returns the derivative of a function f at a point x by Ridders'
// method of polynomial extrapolation. The implementation is taken from
// "Numerical Recipes in C", Chapter 5.7.
template <typename F>
double Diff(F f, double x, double *err = 0) {
  if (gsl_isnan(f(x)))
    return GSL_NAN;
  double h = 0;
  double deriv = Diff(f, x, SymmetricDifference<F>, err, &h);
  double right_err = 0;
  double right_deriv = Diff(f, x, RightDifference<F>, &right_err);
  if (gsl_isnan(deriv)) {
    if (err) *err = right_err;
    return right_deriv;
  }
  double left_deriv = Diff(f, x, LeftDifference<F>);
  if (!(fabs(left_deriv - right_deriv) <= 1e-2))
    return GSL_NAN;
  if (deriv > 1) {
    // Choose h so that x + h and x differ by a number exactly representable
    // as double. See "Numerical Recipes in C", Chapter 5.7.
    double small_h = h / 100;
    volatile double temp = x + h;
    small_h = temp - x;
    // A heuristic to detect infinity.
    double check_deriv = (f(x + small_h) - f(x - small_h)) / (2 * small_h);
    if (check_deriv > deriv * 1.1)
      return GSL_POSINF;
  }
  return deriv;
}

template <typename F>
void CheckDerivative(F f, double x,
    const Function &af, size_t var_index, const ArgList &args) {
  std::ostringstream os;
  os << "Checking d/dx" << var_index << " " << af.name() << " at " << args;
  SCOPED_TRACE(os.str());
  double dx = Diff(f, x);
  if (gsl_isnan(dx)) {
    af(args, ERROR | DERIVS);
    return;
  }
  EXPECT_NEAR(dx, af(args, DERIVS).deriv(var_index), 1e-5);
}

class Deriv {
 private:
  Function af_;
  unsigned deriv_var_index_;
  unsigned eval_var_index_;
  vector<real> args_;
  char *dig_;

 public:
  Deriv(Function af, unsigned eval_var_index,
      unsigned deriv_var_index, const vector<real> &args, char *dig = 0)
  : af_(af), deriv_var_index_(deriv_var_index),
    eval_var_index_(eval_var_index),  args_(args), dig_(dig) {}

  double operator()(double x) {
    args_[eval_var_index_] = x;
    Result result = af_(args_, DERIVS | ERROR_TO_NAN, dig_);
    return gsl_isnan(result) ? GSL_NAN : result.deriv(deriv_var_index_);
  }
};

const unsigned NO_VAR = ~0u;

void CheckSecondDerivatives(const Function &af,
    const ArgList &args, unsigned skip_var = NO_VAR) {
  vector<real> ra(args.get());
  char dig = skip_var != NO_VAR ? skip_var + 1 : 0;
  for (unsigned i = 0, n = ra.size(); i < n; ++i) {
    if (i == skip_var) continue;
    for (unsigned j = 0; j < n; ++j) {
      if (j == skip_var) continue;
      double err = 0;
      double d2 = Diff(Deriv(af, i, j, ra, &dig), ra[i], &err);
      std::ostringstream os;
      os << "Checking if d/dx" << i << " d/dx" << j
          << " " << af.name() << " at " << args << " is " << d2;
      SCOPED_TRACE(os.str());
      if (gsl_isnan(d2)) {
        af(args, ERROR | HES);
        return;
      }
      unsigned ii = i, jj = j;
      if (ii > jj) std::swap(ii, jj);
      unsigned hes_index = ii * (2 * n - ii - 1) / 2 + jj;
      EXPECT_NEAR(d2, af(args, HES, &dig).hes(hes_index),
          err != 0 ? err * 1000 : 1e-10);
    }
  }
}

typedef double (*FuncU)(unsigned);
typedef double (*Func2)(double, double);
typedef double (*Func3)(double, double, double);

class Options {
 public:
  virtual ~Options() {}

  virtual bool HasDerivative(unsigned, const vector<real>&) const {
    return true;
  }

  virtual bool HasDerivative2(const vector<real>&) const {
    return true;
  }
};

class GSLTest : public ::testing::Test {
 protected:
  ASL *asl;

  void SetUp() {
    asl = ASL_alloc(ASL_read_f);
    i_option_ASL = "../solvers/gsl/libamplgsl.so";
    func_add(asl);
  }

  void TearDown() {
    ASL_free(&asl);
  }

  // Get an AMPL function by name.
  Function GetFunction(const char *name) const {
    func_info *fi = func_lookup(asl, name, 0);
    if (!fi)
      throw std::runtime_error(std::string("Function not found: ") + name);
    return Function(asl, fi);
  }

  // Returns true iff expected is at most kMaxUlps ULP's away from
  // actual or both are NaN.  In particular, this function:
  //
  //   - returns false if either number is (but not both) NAN.
  //   - treats really large numbers as almost equal to infinity.
  //   - thinks +0.0 and -0.0 are 0 DLP's apart.
  static bool AlmostEqualOrNaN(double expected, double actual) {
    testing::internal::Double lhs(expected), rhs(actual);
    return (lhs.is_nan() && rhs.is_nan()) || lhs.AlmostEquals(rhs);
  }

  // Test a function taking a single argument.
  template <typename F>
  void TestUnaryFunc(const char *name, F f);
  void TestFunc(const char *name, double (*f)(double)) {
    TestUnaryFunc(name, f);
  }
  void TestFunc(const char *name, double (*f)(double, gsl_mode_t)) {
    TestUnaryFunc(name, std::bind2nd(std::ptr_fun(f), GSL_PREC_DOUBLE));
  }

  // Test a function taking a single argument of type unsigned int.
  void TestFunc(const char *name, FuncU f);

  void TestFunc(const char *name, double (*f)(int, double));
  void TestFunc(const char *name, Func2 f, const Options& opt = Options());
  void TestFunc(const char *name, Func3 f);
};

const double POINTS[] = {-5, -1.23, 0, 1.23, 5};
const size_t NUM_POINTS = sizeof(POINTS) / sizeof(*POINTS);

const double POINTS_FOR_N[] = {
    INT_MIN, INT_MIN + 1, -2, -1, 0, 1, 2, INT_MAX - 1, INT_MAX};
const size_t NUM_POINTS_FOR_N = sizeof(POINTS_FOR_N) / sizeof(*POINTS_FOR_N);

// TODO: remove because if f(x) returns NaN, af(x) should return error
#define EXPECT_ALMOST_EQUAL_OR_NAN(expected, actual) \
    EXPECT_PRED2(AlmostEqualOrNaN, expected, actual)

template <typename F>
void GSLTest::TestUnaryFunc(const char *name, F f) {
  Function af = GetFunction(name);
  for (size_t i = 0; i != NUM_POINTS; ++i) {
    double x = POINTS[i];
    double value = f(x);
    if (gsl_isnan(value)) {
      af(x, ERROR);
      continue;
    }
    EXPECT_ALMOST_EQUAL_OR_NAN(f(x), af(x)) << name << " at " << x;
    double dx = Diff(f, x);
    if (gsl_isnan(dx)) {
      af(x, ERROR | DERIVS);
      continue;
    }
    double actual_deriv = af(x, DERIVS).deriv();
    if (dx != actual_deriv)
      EXPECT_NEAR(dx, actual_deriv, 1e-5) << name << " at " << x;
    CheckSecondDerivatives(af, ArgList(x));
  }
}

void GSLTest::TestFunc(const char *name, FuncU f) {
  Function af = GetFunction(name);
  for (size_t i = 0; i != NUM_POINTS; ++i) {
    double x = POINTS[i];
    if (static_cast<unsigned>(x) != x) {
      af(x, ERROR);
      continue;
    }
    EXPECT_ALMOST_EQUAL_OR_NAN(f(x), af(x)) << name << " at " << x;
    af(x, DERIVS | ERROR);
  }
}

void GSLTest::TestFunc(const char *name, double (*f)(int, double)) {
  Function af = GetFunction(name);
  for (size_t i = 0; i != NUM_POINTS_FOR_N; ++i) {
    for (size_t j = 0; j != NUM_POINTS; ++j) {
      int n = POINTS_FOR_N[i];
      if (n < -100 || n > 100) continue;
      double x = POINTS[j];
      ArgList args(n, x);
      double value = f(n, x);
      if (gsl_isnan(value)) {
        af(args, ERROR);
        continue;
      }
      EXPECT_ALMOST_EQUAL_OR_NAN(value, af(args)) << name << " at " << x;

      af(args, DERIVS | ERROR);
      char dig = 1;
      double deriv = Diff(std::bind1st(std::ptr_fun(f), n), x);
      if (gsl_isnan(deriv)) {
        af(args, ERROR | DERIVS, &dig);
        continue;
      }
      EXPECT_NEAR(deriv, af(args, DERIVS, &dig).deriv(1), 1e-5)
        << name << " at " << n << ", " << x;

      CheckSecondDerivatives(af, args, 0);
    }
  }
}

void GSLTest::TestFunc(const char *name, Func2 f, const Options& opt) {
  Function af = GetFunction(name);
  for (size_t i = 0; i != NUM_POINTS; ++i) {
    for (size_t j = 0; j != NUM_POINTS; ++j) {
      double x = POINTS[i], y = POINTS[j];
      ArgList args(x, y);
      double value = f(x, y);
      if (gsl_isnan(value)) {
        af(args, ERROR);
        continue;
      }
      EXPECT_ALMOST_EQUAL_OR_NAN(value, af(args));
      char dig[2] = {0, 0};
      double dx = opt.HasDerivative(0, args.get()) ?
          Diff(std::bind2nd(std::ptr_fun(f), y), x) : GSL_NAN;
      if (gsl_isnan(dx)) {
        af(args, ERROR | DERIVS);
        dig[0] = 1;
      } else {
        EXPECT_NEAR(dx, af(args, DERIVS, dig).deriv(0), 1e-5)
                << name << " at " << x << ", " << y;
      }

      double dy = opt.HasDerivative(1, args.get()) ?
          Diff(std::bind1st(std::ptr_fun(f), x), y) : GSL_NAN;
      if (gsl_isnan(dy)) {
        af(args, ERROR | DERIVS, dig);
      } else {
        EXPECT_NEAR(dy, af(args, DERIVS, dig).deriv(1), 1e-5)
          << name << " at " << x << ", " << y;
      }

      if (!opt.HasDerivative2(args.get()) || gsl_isnan(dx)) {
        af(args, ERROR | HES);
        continue;
      }
      // TODO
      CheckSecondDerivatives(af, args);
    }
  }
}

class Bind12 {
 private:
  Func3 f_;
  double arg1_;
  double arg2_;

 public:
  Bind12(Func3 f, double arg1, double arg2)
  : f_(f), arg1_(arg1), arg2_(arg2) {}

  double operator()(double arg3) const { return f_(arg1_, arg2_, arg3); }
};

class Bind13 {
 private:
  Func3 f_;
  double arg1_;
  double arg3_;

 public:
  Bind13(Func3 f, double arg1, double arg3)
  : f_(f), arg1_(arg1), arg3_(arg3) {}

  double operator()(double arg2) const { return f_(arg1_, arg2, arg3_); }
};

class Bind23 {
 private:
  Func3 f_;
  double arg2_;
  double arg3_;

 public:
  Bind23(Func3 f, double arg2, double arg3)
  : f_(f), arg2_(arg2), arg3_(arg3) {}

  double operator()(double arg1) const { return f_(arg1, arg2_, arg3_); }
};

void GSLTest::TestFunc(const char *name, Func3 f) {
  Function af = GetFunction(name);
  for (size_t i = 0; i != NUM_POINTS; ++i) {
    for (size_t j = 0; j != NUM_POINTS; ++j) {
      for (size_t k = 0; k != NUM_POINTS; ++k) {
        double x = POINTS[i], y = POINTS[j], z = POINTS[k];
        ArgList args(x, y, z);
        EXPECT_EQ(f(x, y, z), af(args)) << name;
        CheckDerivative(Bind23(f, y, z), x, af, 0, args);
        CheckDerivative(Bind13(f, x, z), y, af, 1, args);
        CheckDerivative(Bind12(f, x, y), z, af, 2, args);
        CheckSecondDerivatives(af, args);
      }
    }
  }
}

#define TEST_FUNC(name) TestFunc("gsl_" #name, gsl_##name);
#define TEST_FUNC_OPT(name) TestFunc("gsl_" #name, gsl_##name, opt);

TEST_F(GSLTest, TestArgList) {
  static const real ARGS[] = {5, 7, 11, 13, 17, 19, 23, 29, 31};
  EXPECT_EQ(vector<real>(ARGS, ARGS + 1), ArgList(5).get());
  EXPECT_EQ(vector<real>(ARGS, ARGS + 2), ArgList(5, 7).get());
  EXPECT_EQ(vector<real>(ARGS, ARGS + 3), ArgList(5, 7, 11).get());
  EXPECT_EQ(vector<real>(ARGS, ARGS + 4), ArgList(5, 7, 11, 13).get());
  EXPECT_EQ(vector<real>(ARGS, ARGS + 5), ArgList(5, 7, 11, 13, 17).get());
  EXPECT_EQ(vector<real>(ARGS, ARGS + 6),
      ArgList(5, 7, 11, 13, 17, 19).get());
  EXPECT_EQ(vector<real>(ARGS, ARGS + 9),
      ArgList(5, 7, 11, 13, 17, 19, 23, 29, 31).get());

  std::ostringstream oss;
  oss << ArgList(3, 5, 7);
  EXPECT_EQ("(3, 5, 7)", oss.str());
}

TEST_F(GSLTest, TestResult) {
  static const real ARGS[] = {5, 7, 11, 13, 17};
  Result r(42, vector<real>(ARGS, ARGS + 2), vector<real>(ARGS + 2, ARGS + 5));
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
}

struct CheckData {
  AmplExports *ae;
  int n;
  int nr;
  vector<real> ra;
  real *derivs;
  real *hes;
  char *dig;
  char *error;
};

real Check(arglist *args) {
  CheckData *data = reinterpret_cast<CheckData*>(args->funcinfo);
  data->ae = args->AE;
  data->n = args->n;
  data->nr = args->nr;
  data->ra = vector<real>(args->ra, args->ra + args->n);
  data->derivs = args->derivs;
  data->hes = args->hes;
  data->dig = args->dig;
  data->error = args->Errmsg;
  if (args->derivs)
    *args->derivs = 123;
  return 42;
}

TEST_F(GSLTest, TestAMPLFunction) {
  ASL testASL = {};
  AmplExports ae = {};
  testASL.i.ae = &ae;
  func_info fi = {};
  fi.funcp = Check;
  Function f(&testASL, &fi);
  CheckData data = {};
  EXPECT_EQ(42, f(777, 0, 0, &data));
  EXPECT_EQ(&ae, data.ae);
  ASSERT_EQ(1, data.n);
  EXPECT_EQ(1, data.nr);
  EXPECT_EQ(777, data.ra[0]);
  EXPECT_TRUE(data.derivs == nullptr);
  EXPECT_TRUE(data.hes == nullptr);
  EXPECT_TRUE(data.dig == nullptr);
  EXPECT_TRUE(data.error == nullptr);
}

TEST_F(GSLTest, TestAMPLFunctionDerivs) {
  ASL testASL = {};
  AmplExports ae = {};
  testASL.i.ae = &ae;
  func_info fi = {};
  fi.funcp = Check;
  Function f(&testASL, &fi);
  CheckData data = {};
  Result res = f(777, DERIVS, 0, &data);
  EXPECT_EQ(42, res);
  EXPECT_EQ(&ae, data.ae);
  ASSERT_EQ(1, data.n);
  EXPECT_EQ(1, data.nr);
  EXPECT_EQ(777, data.ra[0]);
  EXPECT_EQ(123, res.deriv(0));
  EXPECT_TRUE(data.hes == nullptr);
  EXPECT_TRUE(data.dig == nullptr);
  EXPECT_TRUE(data.error == nullptr);
}

TEST_F(GSLTest, Diff) {
  EXPECT_NEAR(1, Diff(sin, 0), 1e-7);
  EXPECT_NEAR(0.25, Diff(sqrt, 4), 1e-7);
  EXPECT_TRUE(gsl_isnan(Diff(
      std::bind2nd(std::ptr_fun(gsl_hypot), 0), 0)));
  EXPECT_NEAR(0, Diff(std::bind2nd(std::ptr_fun(gsl_hypot), -5), 0), 1e-7);
  EXPECT_TRUE(gsl_isnan(Diff(log, 0)));
  EXPECT_NEAR(1, Diff(gsl_log1p, 0), 1e-7);
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
  TEST_FUNC(sf_bessel_Jn);
}

TEST_F(GSLTest, BesselY) {
  TEST_FUNC(sf_bessel_Y0);
  TEST_FUNC(sf_bessel_Y1);
  TEST_FUNC(sf_bessel_Yn);
}

TEST_F(GSLTest, BesselI) {
  TEST_FUNC(sf_bessel_I0);
  TEST_FUNC(sf_bessel_I1);
  TEST_FUNC(sf_bessel_In);
  TEST_FUNC(sf_bessel_I0_scaled);
  TEST_FUNC(sf_bessel_I1_scaled);
  TEST_FUNC(sf_bessel_In_scaled);
}

TEST_F(GSLTest, BesselK) {
  TEST_FUNC(sf_bessel_K0);
  TEST_FUNC(sf_bessel_K1);
  TEST_FUNC(sf_bessel_Kn);
  TEST_FUNC(sf_bessel_K0_scaled);
  TEST_FUNC(sf_bessel_K1_scaled);
  TEST_FUNC(sf_bessel_Kn_scaled);
}

TEST_F(GSLTest, Besselj) {
  TEST_FUNC(sf_bessel_j0);
  TEST_FUNC(sf_bessel_j1);
  TEST_FUNC(sf_bessel_j2);
  TEST_FUNC(sf_bessel_jl);
}

TEST_F(GSLTest, Bessely) {
  TEST_FUNC(sf_bessel_y0);
  TEST_FUNC(sf_bessel_y1);
  TEST_FUNC(sf_bessel_y2);
  TEST_FUNC(sf_bessel_yl);
}

TEST_F(GSLTest, Besseli) {
  TEST_FUNC(sf_bessel_i0_scaled);
  EXPECT_NEAR(0.0999955, gsl_sf_bessel_i0_scaled(5), 1e-5);
  TEST_FUNC(sf_bessel_i1_scaled);
  EXPECT_NEAR(0.0800054, gsl_sf_bessel_i1_scaled(5), 1e-5);
  TEST_FUNC(sf_bessel_i2_scaled);
  EXPECT_NEAR(0.0519922, gsl_sf_bessel_i2_scaled(5), 1e-5);
  TEST_FUNC(sf_bessel_il_scaled);
  EXPECT_NEAR(0.0280133, gsl_sf_bessel_il_scaled(3, 5), 1e-5);
}

TEST_F(GSLTest, Besselk) {
  TEST_FUNC(sf_bessel_k0_scaled);
  EXPECT_NEAR(0.314159, gsl_sf_bessel_k0_scaled(5), 1e-5);
  TEST_FUNC(sf_bessel_k1_scaled);
  EXPECT_NEAR(0.376991, gsl_sf_bessel_k1_scaled(5), 1e-5);
  TEST_FUNC(sf_bessel_k2_scaled);
  EXPECT_NEAR(0.540354, gsl_sf_bessel_k2_scaled(5), 1e-5);
  TEST_FUNC(sf_bessel_kl_scaled);
  EXPECT_NEAR(0.917345, gsl_sf_bessel_kl_scaled(3, 5), 1e-5);
}

class BesselFractionalOrderOptions : public Options {
 public:
  bool HasDerivative(unsigned var_index, const vector<real>& args) const {
    // Computing gsl_sf_bessel_*nu'(nu, x) requires
    // gsl_sf_bessel_*nu(nu - 1, x) which doesn't work when the
    // first argument is non-negative, so nu should be >= 1.
    // Partial derivatives with respect to nu are not provided.
    return var_index == 1 && args[0] >= 1;
  }

  bool HasDerivative2(const vector<real>& args) const {
    // Computing gsl_sf_bessel_*nu''(nu, x) requires
    // gsl_sf_bessel_*nu(nu - 2, x) which doesn't work when the
    // first argument is non-negative, so nu should be >= 2.
    return args[0] >= 2;
  }
};

TEST_F(GSLTest, BesselFractionalOrder) {
  BesselFractionalOrderOptions opt;
  TEST_FUNC_OPT(sf_bessel_Jnu);
  TEST_FUNC_OPT(sf_bessel_Ynu);
  TEST_FUNC_OPT(sf_bessel_Inu);
  TEST_FUNC_OPT(sf_bessel_Inu_scaled);
  TEST_FUNC_OPT(sf_bessel_Knu);
  TEST_FUNC_OPT(sf_bessel_lnKnu);
  TEST_FUNC_OPT(sf_bessel_Knu_scaled);
}

TEST_F(GSLTest, BesselZero) {
  TEST_FUNC(sf_bessel_zero_J0);
  TEST_FUNC(sf_bessel_zero_J1);

  const char *name = "gsl_sf_bessel_zero_Jnu";
  Function af = GetFunction(name);
  for (size_t i = 0; i != NUM_POINTS; ++i) {
    for (size_t j = 0; j != NUM_POINTS; ++j) {
      double nu = POINTS[i], x = POINTS[j];
      ArgList args(nu, x);
      if (static_cast<unsigned>(x) != x) {
        af(args, ERROR);
        continue;
      }
      EXPECT_ALMOST_EQUAL_OR_NAN(gsl_sf_bessel_zero_Jnu(nu, x), af(args))
        << name << " at " << x;
      af(args, DERIVS | ERROR);
    }
  }
}

TEST_F(GSLTest, Clausen) {
  TEST_FUNC(sf_clausen);
}

TEST_F(GSLTest, Hydrogenic) {
  TEST_FUNC(sf_hydrogenicR_1);

  const char *name = "gsl_sf_hydrogenicR";
  Function af = GetFunction(name);
  af(ArgList(0, 0, 0, 0));
  af(ArgList(0.5, 0, 0, 0), ERROR);
  af(ArgList(0, 0.5, 0, 0), ERROR);
  for (size_t i = 0; i != NUM_POINTS_FOR_N; ++i) {
    for (size_t i = 0; i != NUM_POINTS_FOR_N; ++i) {
      for (size_t k = 0; k != NUM_POINTS; ++k) {
        for (size_t l = 0; l != NUM_POINTS; ++l) {
          int n = POINTS_FOR_N[i], ll = POINTS_FOR_N[l];
          if (n < -1000 || n > 1000) continue;
          double Z = POINTS[k], r = POINTS[l];
          ArgList args(n, ll, Z, r);
          EXPECT_ALMOST_EQUAL_OR_NAN(
              gsl_sf_hydrogenicR(n, ll, Z, r), af(args))
            << name << " at " << n;
          af(args, DERIVS | ERROR);
        }
      }
    }
  }
}

TEST_F(GSLTest, Coulomb) {
  const char *name = "gsl_sf_coulomb_CL";
  Function af = GetFunction(name);
  for (size_t i = 0; i != NUM_POINTS; ++i) {
    for (size_t j = 0; j != NUM_POINTS; ++j) {
      double x = POINTS[i], y = POINTS[j];
      ArgList args(x, y);
      gsl_sf_result result = {};
      double value = gsl_sf_coulomb_CL_e(x, y, &result) ? GSL_NAN : result.val;
      EXPECT_ALMOST_EQUAL_OR_NAN(value, af(args));
      af(args, DERIVS | ERROR);
    }
  }
}

TEST_F(GSLTest, Coupling3j) {
  double value = gsl_sf_coupling_3j(8, 20, 12, -2, 12, -10);
  EXPECT_NEAR(0.0812695955, value, 1e-5);
  Function af = GetFunction("gsl_sf_coupling_3j");
  EXPECT_ALMOST_EQUAL_OR_NAN(value, af(ArgList(8, 20, 12, -2, 12, -10)));
  af(ArgList(0, 0, 0, 0, 0, 0));
  af(ArgList(0.5, 0, 0, 0, 0, 0), ERROR);
  af(ArgList(0, 0.5, 0, 0, 0, 0), ERROR);
  af(ArgList(0, 0, 0.5, 0, 0, 0), ERROR);
  af(ArgList(0, 0, 0, 0.5, 0, 0), ERROR);
  af(ArgList(0, 0, 0, 0, 0.5, 0), ERROR);
  af(ArgList(0, 0, 0, 0, 0, 0.5), ERROR);
  af(ArgList(8, 20, 12, -2, 12, -10), ERROR | DERIVS);
}

TEST_F(GSLTest, Coupling6j) {
  double value = gsl_sf_coupling_6j(2, 4, 6, 8, 10, 12);
  EXPECT_NEAR(0.0176295295, value, 1e-5);
  Function af = GetFunction("gsl_sf_coupling_6j");
  EXPECT_ALMOST_EQUAL_OR_NAN(value, af(ArgList(2, 4, 6, 8, 10, 12)));
  af(ArgList(0, 0, 0, 0, 0, 0));
  af(ArgList(0.5, 0, 0, 0, 0, 0), ERROR);
  af(ArgList(0, 0.5, 0, 0, 0, 0), ERROR);
  af(ArgList(0, 0, 0.5, 0, 0, 0), ERROR);
  af(ArgList(0, 0, 0, 0.5, 0, 0), ERROR);
  af(ArgList(0, 0, 0, 0, 0.5, 0), ERROR);
  af(ArgList(0, 0, 0, 0, 0, 0.5), ERROR);
  af(ArgList(2, 4, 6, 8, 10, 12), ERROR | DERIVS);
}

TEST_F(GSLTest, Coupling9j) {
  double value = gsl_sf_coupling_9j(6, 16, 18, 8, 20, 14, 12, 10, 4);
  EXPECT_NEAR(-0.000775648399, value, 1e-9);
  Function af = GetFunction("gsl_sf_coupling_9j");
  return;
  EXPECT_ALMOST_EQUAL_OR_NAN(value,
      af(ArgList(6, 16, 18, 8, 20, 14, 12, 10, 4)));
  af(ArgList(0, 0, 0, 0, 0, 0, 0, 0, 0));
  af(ArgList(0.5, 0, 0, 0, 0, 0, 0, 0, 0), ERROR);
  af(ArgList(0, 0.5, 0, 0, 0, 0, 0, 0, 0), ERROR);
  af(ArgList(0, 0, 0.5, 0, 0, 0, 0, 0, 0), ERROR);
  af(ArgList(0, 0, 0, 0.5, 0, 0, 0, 0, 0), ERROR);
  af(ArgList(0, 0, 0, 0, 0.5, 0, 0, 0, 0), ERROR);
  af(ArgList(0, 0, 0, 0, 0, 0.5, 0, 0, 0), ERROR);
  af(ArgList(0, 0, 0, 0, 0, 0, 0.5, 0, 0), ERROR);
  af(ArgList(0, 0, 0, 0, 0, 0, 0, 0.5, 0), ERROR);
  af(ArgList(0, 0, 0, 0, 0, 0, 0, 0, 0.5), ERROR);
  af(ArgList(6, 16, 18, 8, 20, 14, 12, 10, 4), ERROR | DERIVS);
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

TEST_F(GSLTest, Dilog) {
  TEST_FUNC(sf_dilog);
}

TEST_F(GSLTest, EllInt) {
  TEST_FUNC(sf_ellint_Kcomp);
  TEST_FUNC(sf_ellint_Ecomp);

  ArgList Zero(0, 0);
  ArgList TestPt(0.5, 0.5);
  Function f = GetFunction("gsl_sf_ellint_Pcomp");
  EXPECT_NEAR(1.36647395300460, f(TestPt), 1e-14);
  EXPECT_ALMOST_EQUAL_OR_NAN(GSL_NAN, f(Zero, DERIVS).deriv(0));
  EXPECT_NEAR(0.393428217409760, f(TestPt, DERIVS).deriv(0), 1e-15);
  EXPECT_ALMOST_EQUAL_OR_NAN(GSL_NAN, f(Zero, DERIVS).deriv(1));
  EXPECT_NEAR(-0.471628143501985, f(TestPt, DERIVS).deriv(1), 1e-15);
  EXPECT_ALMOST_EQUAL_OR_NAN(GSL_NAN, f(Zero, HES).hes(0));
  EXPECT_NEAR(1.35115, f(TestPt, HES).hes(0), 1e-5);
  EXPECT_ALMOST_EQUAL_OR_NAN(GSL_NAN, f(Zero, HES).hes(1));
  EXPECT_NEAR(-0.210152, f(TestPt, HES).hes(1), 1e-5);
  EXPECT_ALMOST_EQUAL_OR_NAN(GSL_NAN, f(Zero, HES).hes(2));
  EXPECT_NEAR(0.477835, f(TestPt, HES).hes(2), 1e-5);
}
}
