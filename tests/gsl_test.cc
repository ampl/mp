// GSL wrapper test.

#include <functional>
#include <stdexcept>
#include <vector>

#include <gsl/gsl_math.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_airy.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_clausen.h>
#include <gsl/gsl_sf_coulomb.h>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_dawson.h>
#include <gsl/gsl_sf_debye.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_ellint.h>

#include "gtest/gtest.h"
#include "solvers/asl.h"
#include "tests/config.h"

using std::vector;
using std::runtime_error;

namespace {

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

  vector<real> copy() const { return ra; }
};

// Options for an AMPL function call.
enum {
  ERROR  = 1, // Function call is expected to produce an error.
  DERIVS = 2, // Get first partial derivatives.
  HES    = 6  // Get both first and second partial derivatives.
};

class AMPLResult {
 private:
  real value_;
  vector<real> derivs_;
  vector<real> hes_;

 public:
  AMPLResult() : value_(0) {}

  real *AllocateDerivs(size_t size) {
    derivs_.resize(size);
    return &derivs_[0];
  }

  real *AllocateHes(size_t size) {
    hes_.resize(size * (size + 1) / 2);
    return &hes_[0];
  }

  void SetValue(real value) {
    value_ = value;
  }

  operator real() const { return value_; }

  real deriv(size_t index = 0) const { return derivs_[index]; }
  real hes(size_t index = 0) const { return hes_[index]; }
};

class AMPLFunction {
 private:
  ASL *asl_;
  func_info *fi_;

 public:
  AMPLFunction(ASL *asl, func_info *fi) : asl_(asl), fi_(fi) {}

  AMPLResult operator()(const ArgList &args,
      int options = 0, char *dig = 0, void *info = 0) const {
    // Initialize the argument list.
    arglist al = {};
    TMInfo tmi = {};
    vector<real> ra(args.copy());
    al.ra = &ra[0];
    al.nr = al.n = ra.size();
    al.TMI = &tmi;
    al.AE = asl_->i.ae;
    al.dig = dig;
    al.funcinfo = info;

    // Allocate derivative storage if needed.
    AMPLResult result;
    if ((options & DERIVS) != 0)
      al.derivs = result.AllocateDerivs(ra.size());
    if ((options & HES) == HES)
      al.hes = result.AllocateHes(ra.size());

    // Call the function.
    real value = fi_->funcp(&al);

    // Check the error message.
    if (al.Errmsg) {
      if ((options & ERROR) == 0)
        throw runtime_error(al.Errmsg);
    } else if ((options & ERROR) != 0)
      throw runtime_error(std::string("Expected error in ") + fi_->name);

    result.SetValue(value);
    return result;
  }

  AMPLResult operator()(real x,
      int options = 0, char *dig = 0, void *info = 0) const {
    return (*this)((ArgList(x)), options, dig, info);
  }
};

struct Result {
  double value;
  bool error;

  Result(double val, bool err = false) : value(val), error(err) {}
};

typedef double (*Func1)(double);
typedef double (*FuncU)(unsigned);
typedef double (*FuncN1)(int, double);
typedef Result (*FuncN1Result)(int, double);
typedef double (*Func2)(double, double);
typedef double (*Func3)(double, double, double);

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
  AMPLFunction GetFunction(const char *name) const {
    func_info *fi = func_lookup(asl, name, 0);
    if (!fi)
      throw runtime_error(std::string("Function not found: ") + name);
    return AMPLFunction(asl, fi);
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

  static bool NearOrNaN(double expected, double actual, double abs_error) {
    testing::internal::Double lhs(expected), rhs(actual);
    return (lhs.is_nan() && rhs.is_nan()) ||
        fabs(expected - actual) <= abs_error;
  }

  // Test a function taking a single argument.
  void TestFunc(const char *name, Func1 f, Func1 dx2);

  // Test a function taking a single argument of type unsigned int.
  void TestFunc(const char *name, FuncU f);

  void TestFunc(const char *name, FuncN1 f, FuncN1Result dx2);
  void TestFunc(const char *name, Func2 f, Func2 dx, Func2 dy,
      Func2 dx2, Func2 dxdy, Func2 dy2);
  void TestFunc(const char *name, Func3 f, Func3 dx, Func3 dy, Func3 dz,
      Func3 dx2, Func3 dxdy, Func3 dxdz, Func3 dy2, Func3 dydz, Func3 dz2);
};

const double POINTS[] = {-5, -1.23, 0, 1.23, 5};
const size_t NUM_POINTS = sizeof(POINTS) / sizeof(*POINTS);

const double POINTS_FOR_N[] = {
    INT_MIN, INT_MIN + 1, -2, -1, 0, 1, 2, INT_MAX - 1, INT_MAX};
const size_t NUM_POINTS_FOR_N = sizeof(POINTS_FOR_N) / sizeof(*POINTS_FOR_N);

#define EXPECT_ALMOST_EQUAL_OR_NAN(expected, actual) \
    EXPECT_PRED2(AlmostEqualOrNaN, expected, actual)

#define EXPECT_NEAR_OR_NAN(expected, actual, abs_error) \
    EXPECT_PRED3(NearOrNaN, expected, actual, abs_error)

// Differentiates f numerically.
template <typename F>
double Diff(F f, double x) {
  double h = x != 0 ?
      sqrt(std::numeric_limits<double>::epsilon()) * x : 1e-10;
  double f_x = f(x), f_left = f(x - h), f_right = f(x + h);
  double left_deriv = (f_x - f_left) / h;
  double right_deriv = (f_right - f_x) / h;
  if (gsl_isnan(f_left))
    return right_deriv;
  if (fabs(left_deriv - right_deriv) > 1e-3)
    return GSL_NAN;
  double deriv = (f_right - f_left) / (2 * h);
  if (deriv > 1) {
    // A heuristic to detect infinity.
    double small_h = h / 100;
    double check_deriv = (f(x + small_h) - f(x - small_h)) / (2 * small_h);
    if (check_deriv > deriv * 1.1)
      return GSL_POSINF;
  }
  return deriv;
}

void GSLTest::TestFunc(const char *name, Func1 f, Func1 dx2) {
  AMPLFunction af = GetFunction(name);
  for (size_t i = 0; i != NUM_POINTS; ++i) {
    double x = POINTS[i];
    double value = f(x);
    if (gsl_isnan(value)) {
      af(x, ERROR);
      continue;
    }
    EXPECT_ALMOST_EQUAL_OR_NAN(f(x), af(x)) << name << " at " << x;
    double deriv = Diff(f, x);
    if (gsl_isnan(deriv)) {
      af(x, ERROR | DERIVS);
      continue;
    }
    double actual_deriv = af(x, DERIVS).deriv();
    if (deriv != actual_deriv) {
      EXPECT_NEAR(deriv, actual_deriv, 1e-5)
        << name << " at " << x;
    }
    EXPECT_ALMOST_EQUAL_OR_NAN(dx2(x), af(x, HES).hes())
      << name << " at " << x;
  }
}

void GSLTest::TestFunc(const char *name, FuncU f) {
  AMPLFunction af = GetFunction(name);
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

void GSLTest::TestFunc(const char *name, FuncN1 f, FuncN1Result dx2) {
  AMPLFunction af = GetFunction(name);
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

      Result r = dx2(n, x);
      if (r.error) {
        af(args, HES | ERROR, &dig);
      } else {
        EXPECT_ALMOST_EQUAL_OR_NAN(r.value, af(args, HES, &dig).hes(2))
          << name << " at " << x;
      }
    }
  }
}

void GSLTest::TestFunc(const char *name, Func2 f, Func2 dx, Func2 dy,
    Func2 dx2, Func2 dxdy, Func2 dy2) {
  AMPLFunction af = GetFunction(name);
  for (size_t i = 0; i != NUM_POINTS; ++i) {
    for (size_t j = 0; j != NUM_POINTS; ++j) {
      double x = POINTS[i], y = POINTS[j];
      ArgList args(x, y);
      EXPECT_ALMOST_EQUAL_OR_NAN(f(x, y), af(args));
      char dig = 0;
      AMPLResult res;
      if (dx) {
        res = af(args, DERIVS);
        EXPECT_ALMOST_EQUAL_OR_NAN(dx(x, y), res.deriv(0));
      } else {
        af(args, DERIVS | ERROR);
        dig = 1;
        res = af(args, DERIVS, &dig);
      }
      EXPECT_ALMOST_EQUAL_OR_NAN(dy(x, y), res.deriv(1));
      res = af(args, HES, &dig);
      if (dx2)
        EXPECT_ALMOST_EQUAL_OR_NAN(dx2(x, y), res.hes(0));
      if (dxdy)
        EXPECT_ALMOST_EQUAL_OR_NAN(dxdy(x, y), res.hes(1));
      EXPECT_ALMOST_EQUAL_OR_NAN(dy2(x, y), res.hes(2));
    }
  }
}

void GSLTest::TestFunc(const char *name, Func3 f, Func3 dx, Func3 dy, Func3 dz,
    Func3 dx2, Func3 dxdy, Func3 dxdz, Func3 dy2, Func3 dydz, Func3 dz2) {
  AMPLFunction af = GetFunction(name);
  for (size_t i = 0; i != NUM_POINTS; ++i) {
    for (size_t j = 0; j != NUM_POINTS; ++j) {
      for (size_t k = 0; k != NUM_POINTS; ++k) {
        double x = POINTS[i], y = POINTS[j], z = POINTS[k];
        ArgList args(x, y, z);
        EXPECT_ALMOST_EQUAL_OR_NAN(f(x, y, z), af(args)) << name;
        AMPLResult res = af(args, DERIVS);
        EXPECT_ALMOST_EQUAL_OR_NAN(dx(x, y, z), res.deriv(0)) << name;
        EXPECT_ALMOST_EQUAL_OR_NAN(dy(x, y, z), res.deriv(1)) << name;
        EXPECT_ALMOST_EQUAL_OR_NAN(dz(x, y, z), res.deriv(2)) << name;
        res = af(args, HES);
        EXPECT_ALMOST_EQUAL_OR_NAN(dx2(x, y, z),  res.hes(0)) << name;
        EXPECT_ALMOST_EQUAL_OR_NAN(dxdy(x, y, z), res.hes(1)) << name;
        EXPECT_ALMOST_EQUAL_OR_NAN(dxdz(x, y, z), res.hes(2)) << name;
        EXPECT_ALMOST_EQUAL_OR_NAN(dy2(x, y, z),  res.hes(3)) << name;
        EXPECT_ALMOST_EQUAL_OR_NAN(dydz(x, y, z), res.hes(4)) << name;
        EXPECT_ALMOST_EQUAL_OR_NAN(dz2(x, y, z),  res.hes(5)) << name;
      }
    }
  }
}

double log1p_dx2(double x) { return -1 / ((x + 1) * (x + 1)); }

double expm1_dx2(double x) { return exp(x); }

double hypot_dx(double x, double y) { return x / sqrt(x * x + y * y); }
double hypot_dy(double x, double y) { return y / sqrt(x * x + y * y); }
double hypot_dx2(double x, double y) {
  return y * y / pow(x * x + y * y, 1.5);
}
double hypot_dxdy(double x, double y) {
  return -x * y / pow(x * x + y * y, 1.5);
}
double hypot_dy2(double x, double y) {
  return x * x / pow(x * x + y * y, 1.5);
}

double hypot3_dx(double x, double y, double z) {
  return x / gsl_hypot3(x, y, z);
}
double hypot3_dy(double x, double y, double z) {
  return y / gsl_hypot3(x, y, z);
}
double hypot3_dz(double x, double y, double z) {
  return z / gsl_hypot3(x, y, z);
}
double hypot3_dx2(double x, double y, double z) {
  return (y * y + z * z) / pow(x * x + y * y + z * z, 1.5);
}
double hypot3_dxdy(double x, double y, double z) {
  return - x * y / pow(x * x + y * y + z * z, 1.5);
}
double hypot3_dxdz(double x, double y, double z) {
  return - x * z / pow(x * x + y * y + z * z, 1.5);
}
double hypot3_dy2(double x, double y, double z) {
  return (x * x + z * z) / pow(x * x + y * y + z * z, 1.5);
}
double hypot3_dydz(double x, double y, double z) {
  return - y * z / pow(x * x + y * y + z * z, 1.5);
}
double hypot3_dz2(double x, double y, double z) {
  return (x * x + y * y) / pow(x * x + y * y + z * z, 1.5);
}

double gsl_sf_airy_Ai(double x) {
  return ::gsl_sf_airy_Ai(x, GSL_PREC_DOUBLE);
}
double sf_airy_Ai_dx2(double x) {
  return x * gsl_sf_airy_Ai(x);
}

double gsl_sf_airy_Bi(double x) {
  return ::gsl_sf_airy_Bi(x, GSL_PREC_DOUBLE);
}
double sf_airy_Bi_dx2(double x) {
  return x * gsl_sf_airy_Bi(x);
}

double gsl_sf_airy_Ai_scaled(double x) {
  return ::gsl_sf_airy_Ai_scaled(x, GSL_PREC_DOUBLE);
}
double sf_airy_Ai_scaled_dx2(double x) {
  double dx = x > 0 ? (gsl_sf_airy_Ai_deriv_scaled(x, GSL_PREC_DOUBLE) +
      sqrt(x) * gsl_sf_airy_Ai_scaled(x)) :
      gsl_sf_airy_Ai_deriv(x, GSL_PREC_DOUBLE);
  return x > 0 ? (gsl_sf_airy_Ai_scaled(x) +
      4 * x * dx) / (2 * sqrt(x)) :
      x * gsl_sf_airy_Ai(x);
}

double gsl_sf_airy_Bi_scaled(double x) {
  return ::gsl_sf_airy_Bi_scaled(x, GSL_PREC_DOUBLE);
}
double sf_airy_Bi_scaled_dx2(double x) {
  double dx = x > 0 ? (gsl_sf_airy_Bi_deriv_scaled(x, GSL_PREC_DOUBLE) -
      sqrt(x) * gsl_sf_airy_Bi_scaled(x)) :
      gsl_sf_airy_Bi_deriv(x, GSL_PREC_DOUBLE);
  return x > 0 ? -(gsl_sf_airy_Bi_scaled(x) +
      4 * x * dx) / (2 * sqrt(x)) :
      x * gsl_sf_airy_Bi(x);
}

double sf_bessel_J0_dx2(double x) {
  return 0.5 * (gsl_sf_bessel_Jn(2, x) - gsl_sf_bessel_J0(x));
}

double sf_bessel_J1_dx2(double x) {
  return 0.25 * (gsl_sf_bessel_Jn(3, x) - 3 * gsl_sf_bessel_J1(x));
}

Result sf_bessel_Jn_dx2(int n, double x) {
  if (n <= INT_MIN + 1 || n >= INT_MAX - 1)
    return Result(0, true);
  return 0.25 * (gsl_sf_bessel_Jn(n - 2, x) -
      2 * gsl_sf_bessel_Jn(n, x) + gsl_sf_bessel_Jn(n + 2, x));
}

double sf_bessel_Y0_dx2(double x) {
  return 0.5 * (gsl_sf_bessel_Yn(2, x) - gsl_sf_bessel_Y0(x));
}

double sf_bessel_Y1_dx2(double x) {
  return 0.25 * (gsl_sf_bessel_Yn(3, x) - 3 * gsl_sf_bessel_Y1(x));
}

Result sf_bessel_Yn_dx2(int n, double x) {
  if (n <= INT_MIN + 1 || n >= INT_MAX - 1)
    return Result(0, true);
  return 0.25 * (gsl_sf_bessel_Yn(n - 2, x) -
      2 * gsl_sf_bessel_Yn(n, x) + gsl_sf_bessel_Yn(n + 2, x));
}

double sf_bessel_I0_dx2(double x) {
  return 0.5 * (gsl_sf_bessel_In(2, x) + gsl_sf_bessel_I0(x));
}

double sf_bessel_I0_scaled_dx2(double x) {
  return 1.5 * gsl_sf_bessel_I0_scaled(x) -
      2 * x / fabs(x) * gsl_sf_bessel_I1_scaled(x) +
      0.5 * gsl_sf_bessel_In_scaled(2, x);
}

double sf_bessel_I1_dx2(double x) {
  return 0.25 * (gsl_sf_bessel_In(3, x) + 3 * gsl_sf_bessel_I1(x));
}

Result sf_bessel_In_dx2(int n, double x) {
  if (n <= INT_MIN + 1 || n >= INT_MAX - 1)
    return Result(0, true);
  return 0.25 * (gsl_sf_bessel_In(n - 2, x) +
      2 * gsl_sf_bessel_In(n, x) + gsl_sf_bessel_In(n + 2, x));
}

double sf_bessel_I1_scaled_dx2(double x) {
  return -x / fabs(x) * gsl_sf_bessel_I0_scaled(x) +
      1.75 * gsl_sf_bessel_I1_scaled(x) -
      x / fabs(x) * gsl_sf_bessel_In_scaled(2, x) +
      0.25 * gsl_sf_bessel_In_scaled(3, x);
}

Result sf_bessel_In_scaled_dx2(int n, double x) {
  if (n <= INT_MIN + 1 || n >= INT_MAX - 1)
    return Result(0, true);
  return (fabs(x) * (gsl_sf_bessel_In_scaled(n - 2, x) +
                    6 * gsl_sf_bessel_In_scaled(n, x) +
                    gsl_sf_bessel_In_scaled(n + 2, x)) -
      4 * x * (gsl_sf_bessel_In_scaled(n - 1, x) +
               gsl_sf_bessel_In_scaled(n + 1, x))) / (4 * fabs(x));
}

double sf_bessel_K0_dx2(double x) {
  return 0.5 * (gsl_sf_bessel_Kn(2, x) + gsl_sf_bessel_K0(x));
}

double sf_bessel_K1_dx2(double x) {
  return 0.25 * (gsl_sf_bessel_Kn(3, x) + 3 * gsl_sf_bessel_K1(x));
}

Result sf_bessel_Kn_dx2(int n, double x) {
  if (n <= INT_MIN + 1 || n >= INT_MAX - 1)
    return Result(0, true);
  return 0.25 * (gsl_sf_bessel_Kn(n - 2, x) +
      2 * gsl_sf_bessel_Kn(n, x) + gsl_sf_bessel_Kn(n + 2, x));
}

double sf_bessel_K0_scaled_dx2(double x) {
  return 1.5 * gsl_sf_bessel_K0_scaled(x) - 2 * gsl_sf_bessel_K1_scaled(x) +
      0.5 * gsl_sf_bessel_Kn_scaled(2, x);
}

double sf_bessel_K1_scaled_dx2(double x) {
  return -gsl_sf_bessel_K0_scaled(x) + 1.75 * gsl_sf_bessel_K1_scaled(x) -
      gsl_sf_bessel_Kn_scaled(2, x) + 0.25 * gsl_sf_bessel_Kn_scaled(3, x);
}

Result sf_bessel_Kn_scaled_dx2(int n, double x) {
  if (n <= INT_MIN + 1 || n >= INT_MAX - 1)
    return Result(0, true);
  return 0.25 * (gsl_sf_bessel_Kn_scaled(n - 2, x) -
      4 * gsl_sf_bessel_Kn_scaled(n - 1, x) +
      6 * gsl_sf_bessel_Kn_scaled(n, x) -
      4 * gsl_sf_bessel_Kn_scaled(n + 1, x) +
      gsl_sf_bessel_Kn_scaled(n + 2, x));
}

double sf_bessel_j0_dx2(double x) {
  return -((x * x - 2) * sin(x) + 2 * x * cos(x)) / (x * x * x);
}

double sf_bessel_j1_dx2(double x) {
  return (x * (x * x - 6) * cos(x) -
      3 * (x * x - 2) * sin(x)) / (x * x * x * x);
}

double sf_bessel_j2_dx2(double x) {
  double x_squared = x * x;
  return (x * (5 * x_squared - 36) * cos(x) +
      (x_squared * x_squared - 17 * x_squared + 36) * sin(x)) /
      (x_squared * x_squared * x);
}

Result sf_bessel_jl_dx2(int n, double x) {
  if (n <= INT_MIN + 1 || n >= INT_MAX - 1)
    return Result(0, true);
  return (x * x * gsl_sf_bessel_jl(n - 2, x) -
      2 * x * x * gsl_sf_bessel_jl(n, x) +
      x * x * gsl_sf_bessel_jl(n + 2, x) -
      2 * x * gsl_sf_bessel_jl(n - 1, x) +
      2 * x * gsl_sf_bessel_jl(n + 1, x) +
      3 * gsl_sf_bessel_jl(n, x)) / (4 * x * x);
}

double sf_bessel_y0_dx2(double x) {
  return ((x * x - 2) * cos(x) - 2 * x * sin(x)) / (x * x * x);
}

double sf_bessel_y1_dx2(double x) {
  return (x * (x * x - 6) * sin(x) +
      3 * (x * x - 2) * cos(x)) / (x * x * x * x);
}

double sf_bessel_y2_dx2(double x) {
  return ((36 - 5 * x * x) * gsl_sf_bessel_y1(x) -
      (x * x - 12) * cos(x)) / (x * x * x);
}

Result sf_bessel_yl_dx2(int n, double x) {
  if (n <= INT_MIN + 1 || n >= INT_MAX - 1)
    return Result(0, true);
  return (x * x * gsl_sf_bessel_yl(n - 2, x) -
      2 * x * x * gsl_sf_bessel_yl(n, x) +
      x * x * gsl_sf_bessel_yl(n + 2, x) -
      2 * x * gsl_sf_bessel_yl(n - 1, x) +
      2 * x * gsl_sf_bessel_yl(n + 1, x) +
      3 * gsl_sf_bessel_yl(n, x)) / (4 * x * x);
}

double sf_bessel_i0_scaled_dx2(double x) {
  double coef = -2 * (1 + 2 * fabs(x)) / x;
  double hyp_coef = exp(-fabs(x)) / x;
  double i_minus1 = hyp_coef * cosh(x);
  double i_minus2 = hyp_coef * sinh(x) - i_minus1 / x;
  return 0.25 * (
      i_minus2 +
      coef * i_minus1 +
      (3 + 6 * x * x + 4 * fabs(x)) * gsl_sf_bessel_i0_scaled(x) / (x * x) +
      coef * gsl_sf_bessel_i1_scaled(x) +
      gsl_sf_bessel_i2_scaled(x));
}

double sf_bessel_i1_scaled_dx2(double x) {
  double coef = -2 * (1 + 2 * fabs(x)) / x;
  double i_minus1 = (exp(-fabs(x)) * sqrt(1 / x) * cosh(x)) / sqrt(x);
  return 0.25 * (
      i_minus1 +
      coef * gsl_sf_bessel_i0_scaled(x) +
      (3 + 6 * x * x + 4 * fabs(x)) * gsl_sf_bessel_i1_scaled(x) / (x * x) +
      coef * gsl_sf_bessel_i2_scaled(x) +
      gsl_sf_bessel_il_scaled(3, x));
}

double sf_bessel_i2_scaled_dx2(double x) {
  double coef = -2 * (1 + 2 * fabs(x)) / x;
  return 0.25 * (
      gsl_sf_bessel_i0_scaled(x) +
      coef * gsl_sf_bessel_i1_scaled(x) +
      (3 + 6 * x * x + 4 * fabs(x)) * gsl_sf_bessel_il_scaled(2, x) / (x * x) +
      coef * gsl_sf_bessel_il_scaled(3, x) +
      gsl_sf_bessel_il_scaled(4, x));
}

Result sf_bessel_il_scaled_dx2(int n, double x) {
  if (n <= INT_MIN + 1 || n >= INT_MAX - 1)
    return Result(0, true);
  double coef = -2 * (1 + 2 * fabs(x)) / x;
  return 0.25 * (
      gsl_sf_bessel_il_scaled(n - 2, x) +
      coef * gsl_sf_bessel_il_scaled(n - 1, x) +
      (3 + 4 * fabs(x) + 6 * x * x) * gsl_sf_bessel_il_scaled(n, x) / (x * x) +
      coef * gsl_sf_bessel_il_scaled(n + 1, x) +
      gsl_sf_bessel_il_scaled(n + 2, x));
}

double sf_bessel_k0_scaled_dx2(double x) {
  return M_PI * sqrt(1 / x) / pow(x, 2.5);
}

double sf_bessel_k1_scaled_dx2(double x) {
  return (M_PI * sqrt(1 / x) * (x + 3)) / pow(x, 3.5);
}

double sf_bessel_k2_scaled_dx2(double x) {
  return (M_PI * sqrt(1 / x) * (x * x + 9 * x + 18)) / pow(x, 4.5);
}

Result sf_bessel_kl_scaled_dx2(int n, double x) {
  if (n <= INT_MIN + 1 || n >= INT_MAX - 1)
    return Result(0, true);
  double coef = 2 * (1 - 2 * x) / x;
  return 0.25 * (
      gsl_sf_bessel_kl_scaled(n - 2, x) +
      coef * gsl_sf_bessel_kl_scaled(n - 1, x) +
      (3 - 4 * x + 6 * x * x) * gsl_sf_bessel_kl_scaled(n, x) / (x * x) +
      coef * gsl_sf_bessel_kl_scaled(n + 1, x) +
      gsl_sf_bessel_kl_scaled(n + 2, x));
}

Func2 sf_bessel_Jnu_dx, sf_bessel_Jnu_dx2, sf_bessel_Jnu_dxdy;
double sf_bessel_Jnu_dy(double x, double y) {
  return 0.5 * (gsl_sf_bessel_Jnu(x - 1, y) - gsl_sf_bessel_Jnu(x + 1, y));
}
double sf_bessel_Jnu_dy2(double x, double y) {
  return 0.25 * (gsl_sf_bessel_Jnu(x - 2, y) -
      2 * gsl_sf_bessel_Jnu(x, y) + gsl_sf_bessel_Jnu(x + 2, y));
}

Func2 sf_bessel_Ynu_dx, sf_bessel_Ynu_dx2, sf_bessel_Ynu_dxdy;
double sf_bessel_Ynu_dy(double x, double y) {
  return 0.5 * (gsl_sf_bessel_Ynu(x - 1, y) - gsl_sf_bessel_Ynu(x + 1, y));
}
double sf_bessel_Ynu_dy2(double x, double y) {
  return 0.25 * (gsl_sf_bessel_Ynu(x - 2, y) -
      2 * gsl_sf_bessel_Ynu(x, y) + gsl_sf_bessel_Ynu(x + 2, y));
}

Func2 sf_bessel_Inu_dx, sf_bessel_Inu_dx2, sf_bessel_Inu_dxdy;
double sf_bessel_Inu_dy(double x, double y) {
  return 0.5 * (gsl_sf_bessel_Inu(x - 1, y) + gsl_sf_bessel_Inu(x + 1, y));
}
double sf_bessel_Inu_dy2(double x, double y) {
  return 0.25 * (gsl_sf_bessel_Inu(x - 2, y) +
      2 * gsl_sf_bessel_Inu(x, y) + gsl_sf_bessel_Inu(x + 2, y));
}

Func2 sf_bessel_Inu_scaled_dx;
Func2 sf_bessel_Inu_scaled_dx2, sf_bessel_Inu_scaled_dxdy;
double sf_bessel_Inu_scaled_dy(double x, double y) {
  return 0.5 * (-(2 * y * gsl_sf_bessel_Inu_scaled(x, y)) / fabs(y) +
      gsl_sf_bessel_Inu_scaled(x - 1, y) + gsl_sf_bessel_Inu_scaled(x + 1, y));
}
double sf_bessel_Inu_scaled_dy2(double x, double y) {
  return (fabs(y) * (gsl_sf_bessel_Inu_scaled(x - 2, y) +
                    6 * gsl_sf_bessel_Inu_scaled(x, y) +
                    gsl_sf_bessel_Inu_scaled(x + 2, y)) -
      4 * y * (gsl_sf_bessel_Inu_scaled(x - 1, y) +
               gsl_sf_bessel_Inu_scaled(x + 1, y))) / (4 * fabs(y));
}

Func2 sf_bessel_Knu_dx, sf_bessel_Knu_dx2, sf_bessel_Knu_dxdy;
double sf_bessel_Knu_dy(double x, double y) {
  return -0.5 * (gsl_sf_bessel_Knu(x - 1, y) + gsl_sf_bessel_Knu(x + 1, y));
}
double sf_bessel_Knu_dy2(double x, double y) {
  return 0.25 * (gsl_sf_bessel_Knu(x - 2, y) +
      2 * gsl_sf_bessel_Knu(x, y) + gsl_sf_bessel_Knu(x + 2, y));
}

Func2 sf_bessel_lnKnu_dx, sf_bessel_lnKnu_dx2, sf_bessel_lnKnu_dxdy;
double sf_bessel_lnKnu_dy(double x, double y) {
  return -0.5 * (gsl_sf_bessel_Knu(x - 1, y) + gsl_sf_bessel_Knu(x + 1, y)) /
      gsl_sf_bessel_Knu(x, y);
}
double sf_bessel_lnKnu_dy2(double x, double y) {
  double kn = gsl_sf_bessel_Knu(x, y);
  double kn_plus_minus =
      gsl_sf_bessel_Knu(x - 1, y) + gsl_sf_bessel_Knu(x + 1, y);
  return 0.25 * (kn * (gsl_sf_bessel_Knu(x - 2, y) + 2 * kn +
      gsl_sf_bessel_Knu(x + 2, y)) - kn_plus_minus * kn_plus_minus) / (kn * kn);
}

Func2 sf_bessel_Knu_scaled_dx;
Func2 sf_bessel_Knu_scaled_dx2, sf_bessel_Knu_scaled_dxdy;
double sf_bessel_Knu_scaled_dy(double x, double y) {
  return -0.5 * (gsl_sf_bessel_Knu_scaled(x - 1, y) -
      2 * gsl_sf_bessel_Knu_scaled(x, y) + gsl_sf_bessel_Knu_scaled(x + 1, y));
}
double sf_bessel_Knu_scaled_dy2(double x, double y) {
  return 0.25 * (gsl_sf_bessel_Knu_scaled(x - 2, y) -
      4 * gsl_sf_bessel_Knu_scaled(x - 1, y) +
      6 * gsl_sf_bessel_Knu_scaled(x, y) -
      4 * gsl_sf_bessel_Knu_scaled(x + 1, y) +
      gsl_sf_bessel_Knu_scaled(x + 2, y));
}

double sf_clausen_dx2(double x) {
  return -0.5 * tan(0.5 * M_PI - x);
}

double sf_hydrogenicR_1_dx(double x, double y) {
  double Z = x, r = y;
  return sqrt(Z) *  exp(-Z * r) * (3 - 2 * r * Z);
}
double sf_hydrogenicR_1_dy(double x, double y) {
  double Z = x, r = y;
  return -2 * pow(Z, 2.5) * exp(-Z * r);
}
double sf_hydrogenicR_1_dx2(double x, double y) {
  double Z = x, r = y;
  return (exp(-Z * r) * (4 * r * r * Z * Z - 12 * r * Z + 3)) / (2 *sqrt(Z));
}
double sf_hydrogenicR_1_dxdy(double x, double y) {
  double Z = x, r = y;
  return pow(Z, 1.5) * exp(-Z * r) * (2 * r * Z - 5);
}
double sf_hydrogenicR_1_dy2(double x, double y) {
  double Z = x, r = y;
  return 2 * pow(Z, 3.5) * exp(-Z * r);
}

double sf_dawson_dx2(double x) {
  double dx = 1 - 2 * x * gsl_sf_dawson(x);
  return - 2 * (gsl_sf_dawson(x) + x * dx);
}

double debye_dx(int n, double x, double (*func)(double)) {
  return n * (1 / (exp(x) - 1) - func(x) / x);
}
double debye_dx2(int n, double x, double (*func)(double)) {
  return n * (-exp(x) / ((exp(x) - 1) * (exp(x) - 1)) +
      func(x) / (x * x) - debye_dx(n, x, func) / x);
}

double sf_debye_1_dx2(double x) { return debye_dx2(1, x, gsl_sf_debye_1); }
double sf_debye_2_dx2(double x) { return debye_dx2(2, x, gsl_sf_debye_2); }
double sf_debye_3_dx2(double x) { return debye_dx2(3, x, gsl_sf_debye_3); }
double sf_debye_4_dx2(double x) { return debye_dx2(4, x, gsl_sf_debye_4); }
double sf_debye_5_dx2(double x) { return debye_dx2(5, x, gsl_sf_debye_5); }
double sf_debye_6_dx2(double x) { return debye_dx2(6, x, gsl_sf_debye_6); }

double sf_dilog_dx2(double x) {
  return x != 0 ? (1 / (1 - x) +
      GSL_REAL(gsl_complex_log(gsl_complex_rect(1 - x, 0))) / x) / x : 0.5;
}

double gsl_sf_ellint_Kcomp(double x) {
  return ::gsl_sf_ellint_Kcomp(x, GSL_PREC_DOUBLE);
}
double gsl_sf_ellint_Ecomp(double x) {
  return ::gsl_sf_ellint_Ecomp(x, GSL_PREC_DOUBLE);
}

double sf_ellint_Kcomp_dx2(double x) {
  double x_squared = x * x;
  double div = x * (1 - x_squared);
  double K = gsl_sf_ellint_Kcomp(x);
  double E = gsl_sf_ellint_Ecomp(x);
  return x != 0 ? ((2 * x_squared * x_squared - 3 * x_squared + 1) * K +
      (3 * x_squared - 1) * E) / (div * div) : M_PI_4;
}

double sf_ellint_Ecomp_dx2(double x) {
  double K = gsl_sf_ellint_Kcomp(x);
  double E = gsl_sf_ellint_Ecomp(x);
  return x != 0 ?
      ((x * x - 1) * K + E) / (x * x * (x * x - 1)) : -M_PI_4;
}

#define TEST_FUNC(name) \
  TestFunc("gsl_" #name, gsl_##name, name##_dx2);

#define TEST_FUNC_U(name) \
  TestFunc("gsl_" #name, gsl_##name);

#define TEST_FUNC_N(name) \
  TestFunc("gsl_" #name, gsl_##name, name##_dx2);

#define TEST_FUNC2(name) \
    TestFunc("gsl_" #name, gsl_##name, name##_dx, name##_dy, \
        name##_dx2, name##_dxdy, name##_dy2);

#define TEST_FUNC3(name) \
    TestFunc("gsl_" #name, gsl_##name, name##_dx, name##_dy, name##_dz, \
        name##_dx2, name##_dxdy, name##_dxdz, \
        name##_dy2, name##_dydz, name##_dz2);

TEST_F(GSLTest, TestArgList) {
  static const real ARGS[] = {5, 7, 11, 13, 17, 19, 23, 29, 31};
  EXPECT_EQ(vector<real>(ARGS, ARGS + 1), ArgList(5).copy());
  EXPECT_EQ(vector<real>(ARGS, ARGS + 2), ArgList(5, 7).copy());
  EXPECT_EQ(vector<real>(ARGS, ARGS + 3), ArgList(5, 7, 11).copy());
  EXPECT_EQ(vector<real>(ARGS, ARGS + 4), ArgList(5, 7, 11, 13).copy());
  EXPECT_EQ(vector<real>(ARGS, ARGS + 5), ArgList(5, 7, 11, 13, 17).copy());
  EXPECT_EQ(vector<real>(ARGS, ARGS + 6),
      ArgList(5, 7, 11, 13, 17, 19).copy());
  EXPECT_EQ(vector<real>(ARGS, ARGS + 9),
      ArgList(5, 7, 11, 13, 17, 19, 23, 29, 31).copy());
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
  AMPLFunction f(&testASL, &fi);
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
  AMPLFunction f(&testASL, &fi);
  CheckData data = {};
  AMPLResult res = f(777, DERIVS, 0, &data);
  EXPECT_EQ(42, res);
  EXPECT_EQ(&ae, data.ae);
  ASSERT_EQ(1, data.n);
  EXPECT_EQ(1, data.nr);
  EXPECT_EQ(777, data.ra[0]);
  EXPECT_EQ(123, *data.derivs);
  EXPECT_EQ(123, res.deriv(0));
  EXPECT_TRUE(data.hes == nullptr);
  EXPECT_TRUE(data.dig == nullptr);
  EXPECT_TRUE(data.error == nullptr);
}

TEST_F(GSLTest, Elementary) {
  TEST_FUNC(log1p);
  TEST_FUNC(expm1);
  TEST_FUNC2(hypot);
  TEST_FUNC3(hypot3);
}

TEST_F(GSLTest, AiryA) {
  TEST_FUNC(sf_airy_Ai);
  TEST_FUNC(sf_airy_Ai_scaled);
  EXPECT_NEAR(0.00207512, sf_airy_Ai_scaled_dx2(5), 1e-5);
}

TEST_F(GSLTest, AiryB) {
  TEST_FUNC(sf_airy_Bi);
  TEST_FUNC(sf_airy_Bi_scaled);
  EXPECT_NEAR(0.00559418, sf_airy_Bi_scaled_dx2(5), 1e-5);
}

TEST_F(GSLTest, AiryZero) {
  TEST_FUNC_U(sf_airy_zero_Ai);
  TEST_FUNC_U(sf_airy_zero_Bi);
  TEST_FUNC_U(sf_airy_zero_Ai_deriv);
  TEST_FUNC_U(sf_airy_zero_Bi_deriv);
}

TEST_F(GSLTest, BesselJ) {
  TEST_FUNC(sf_bessel_J0);
  TEST_FUNC(sf_bessel_J1);
  TEST_FUNC_N(sf_bessel_Jn);
  EXPECT_NEAR(-0.199025, sf_bessel_Jn_dx2(3, 5).value, 1e-5);
}

TEST_F(GSLTest, BesselY) {
  TEST_FUNC(sf_bessel_Y0);
  TEST_FUNC(sf_bessel_Y1);
  TEST_FUNC_N(sf_bessel_Yn);
  EXPECT_NEAR(-0.149592, sf_bessel_Yn_dx2(3, 5).value, 1e-5);
}

TEST_F(GSLTest, BesselI) {
  TEST_FUNC(sf_bessel_I0);
  TEST_FUNC(sf_bessel_I1);
  TEST_FUNC_N(sf_bessel_In);
  EXPECT_NEAR(11.78897, sf_bessel_In_dx2(3, 5).value, 1e-5);

  TEST_FUNC(sf_bessel_I0_scaled);
  EXPECT_NEAR(0.00634264, sf_bessel_I0_scaled_dx2(5), 1e-5);

  TEST_FUNC(sf_bessel_I1_scaled);
  EXPECT_NEAR(0.00286143, sf_bessel_I1_scaled_dx2(5), 1e-5);

  TEST_FUNC_N(sf_bessel_In_scaled);
  EXPECT_NEAR(-0.00332666, sf_bessel_In_scaled_dx2(3, 5).value, 1e-5);
}

TEST_F(GSLTest, BesselK) {
  TEST_FUNC(sf_bessel_K0);
  TEST_FUNC(sf_bessel_K1);
  TEST_FUNC_N(sf_bessel_Kn);
  EXPECT_NEAR(0.0133336, sf_bessel_Kn_dx2(3, 5).value, 1e-5);

  TEST_FUNC(sf_bessel_K0_scaled);
  EXPECT_NEAR(0.0151222, sf_bessel_K0_scaled_dx2(5), 1e-5);

  TEST_FUNC(sf_bessel_K1_scaled);
  EXPECT_NEAR(0.0224065, sf_bessel_K1_scaled_dx2(5), 1e-5);

  TEST_FUNC_N(sf_bessel_Kn_scaled);
  EXPECT_NEAR(0.156927, sf_bessel_Kn_scaled_dx2(3, 5).value, 1e-5);
}

TEST_F(GSLTest, Besselj) {
  TEST_FUNC(sf_bessel_j0);
  EXPECT_NEAR(0.153749, sf_bessel_j0_dx2(5), 1e-5);

  TEST_FUNC(sf_bessel_j1);
  EXPECT_NEAR(0.148982, sf_bessel_j1_dx2(5), 1e-5);

  TEST_FUNC(sf_bessel_j2);
  EXPECT_NEAR(-0.0320245, sf_bessel_j2_dx2(5), 1e-5);

  TEST_FUNC_N(sf_bessel_jl);
  EXPECT_NEAR(-0.0998566, sf_bessel_jl_dx2(3, 5).value, 1e-5);
}

TEST_F(GSLTest, Bessely) {
  TEST_FUNC(sf_bessel_y0);
  EXPECT_NEAR(0.128908, sf_bessel_y0_dx2(5), 1e-5);

  TEST_FUNC(sf_bessel_y1);
  EXPECT_NEAR(-0.11444, sf_bessel_y1_dx2(5), 1e-5);

  TEST_FUNC(sf_bessel_y2);
  EXPECT_NEAR(-0.157973, sf_bessel_y2_dx2(5), 1e-5);

  TEST_FUNC_N(sf_bessel_yl);
  EXPECT_NEAR(-0.0629096, sf_bessel_yl_dx2(3, 5).value, 1e-5);
}

TEST_F(GSLTest, Besseli) {
  TEST_FUNC(sf_bessel_i0_scaled);
  EXPECT_NEAR(0.0999955, gsl_sf_bessel_i0_scaled(5), 1e-5);
  EXPECT_NEAR(0.00797784, sf_bessel_i0_scaled_dx2(5), 1e-5);

  TEST_FUNC(sf_bessel_i1_scaled);
  EXPECT_NEAR(0.0800054, gsl_sf_bessel_i1_scaled(5), 1e-5);
  EXPECT_NEAR(0.00322746, sf_bessel_i1_scaled_dx2(5), 1e-5);

  TEST_FUNC(sf_bessel_i2_scaled);
  EXPECT_NEAR(0.0519922, gsl_sf_bessel_i2_scaled(5), 1e-5);
  EXPECT_NEAR(-0.000681812, sf_bessel_i2_scaled_dx2(5), 1e-5);

  TEST_FUNC_N(sf_bessel_il_scaled);
  EXPECT_NEAR(0.0280133, gsl_sf_bessel_il_scaled(3, 5), 1e-5);
  EXPECT_NEAR(-0.00152293, sf_bessel_il_scaled_dx2(3, 5).value, 1e-5);
}

TEST_F(GSLTest, Besselk) {
  TEST_FUNC(sf_bessel_k0_scaled);
  EXPECT_NEAR(0.314159, gsl_sf_bessel_k0_scaled(5), 1e-5);
  EXPECT_NEAR(0.0251327, sf_bessel_k0_scaled_dx2(5), 1e-5);

  TEST_FUNC(sf_bessel_k1_scaled);
  EXPECT_NEAR(0.376991, gsl_sf_bessel_k1_scaled(5), 1e-5);
  EXPECT_NEAR(0.0402124, sf_bessel_k1_scaled_dx2(5), 1e-5);

  TEST_FUNC(sf_bessel_k2_scaled);
  EXPECT_NEAR(0.540354, gsl_sf_bessel_k2_scaled(5), 1e-5);
  EXPECT_NEAR(0.0884672, sf_bessel_k2_scaled_dx2(5), 1e-5);

  TEST_FUNC_N(sf_bessel_kl_scaled);
  EXPECT_NEAR(0.917345, gsl_sf_bessel_kl_scaled(3, 5), 1e-5);
  EXPECT_NEAR(0.236248, sf_bessel_kl_scaled_dx2(3, 5).value, 1e-5);
}

TEST_F(GSLTest, BesselFractionalOrder) {
  TEST_FUNC2(sf_bessel_Jnu);
  EXPECT_NEAR(0.410029, gsl_sf_bessel_Jnu(3.5, 5), 1e-5);
  EXPECT_NEAR(-0.0466428, sf_bessel_Jnu_dy(3.5, 5), 1e-5);
  EXPECT_NEAR(-0.199786, sf_bessel_Jnu_dy2(3.5, 5), 1e-5);

  TEST_FUNC2(sf_bessel_Ynu);
  EXPECT_NEAR(-0.0275521, gsl_sf_bessel_Ynu(3.5, 5), 1e-5);
  EXPECT_NEAR(0.313659, sf_bessel_Ynu_dy(3.5, 5), 1e-5);
  EXPECT_NEAR(-0.0486802, sf_bessel_Ynu_dy2(3.5, 5), 1e-5);

  TEST_FUNC2(sf_bessel_Inu);
  EXPECT_NEAR(7.417560126111555, gsl_sf_bessel_Inu(3.5, 5), 1e-5);
  EXPECT_NEAR(8.574590050404494, sf_bessel_Inu_dy(3.5, 5), 1e-5);
  EXPECT_NEAR(9.337246577825318, sf_bessel_Inu_dy2(3.5, 5), 1e-5);

  TEST_FUNC2(sf_bessel_Inu_scaled);
  EXPECT_NEAR(0.04997912699226937, gsl_sf_bessel_Inu_scaled(3.5, 5), 1e-5);
  EXPECT_NEAR(0.00779600630624169, sf_bessel_Inu_scaled_dy(3.5, 5), 1e-5);
  EXPECT_NEAR(-0.0026572670459736, sf_bessel_Inu_scaled_dy2(3.5, 5), 1e-5);

  TEST_FUNC2(sf_bessel_Knu);
  EXPECT_NEAR(0.011027711053957217, gsl_sf_bessel_Knu(3.5, 5), 1e-5);
  EXPECT_NEAR(-0.014215172742155810, sf_bessel_Knu_dy(3.5, 5), 1e-5);
  EXPECT_NEAR(0.01927432401882742, sf_bessel_Knu_dy2(3.5, 5), 1e-5);

  TEST_FUNC2(sf_bessel_lnKnu);
  EXPECT_NEAR(-4.5073439872921324, gsl_sf_bessel_lnKnu(3.5, 5), 1e-5);
  EXPECT_NEAR(-1.289041095890411, sf_bessel_lnKnu_dy(3.5, 5), 1e-5);
  EXPECT_NEAR(0.08618127228373, sf_bessel_lnKnu_dy2(3.5, 5), 1e-5);

  TEST_FUNC2(sf_bessel_Knu_scaled);
  EXPECT_NEAR(1.6366574351881952, gsl_sf_bessel_Knu_scaled(3.5, 5), 1e-5);
  EXPECT_NEAR(-0.4730612586639852, sf_bessel_Knu_scaled_dy(3.5, 5), 1e-5);
  EXPECT_NEAR(0.2777833646846813, sf_bessel_Knu_scaled_dy2(3.5, 5), 1e-5);
}

TEST_F(GSLTest, BesselZero) {
  TEST_FUNC_U(sf_bessel_zero_J0);
  TEST_FUNC_U(sf_bessel_zero_J1);

  const char *name = "gsl_sf_bessel_zero_Jnu";
  AMPLFunction af = GetFunction(name);
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
  TEST_FUNC2(sf_hydrogenicR_1);
  EXPECT_NEAR(0.0633592, gsl_sf_hydrogenicR_1(3, 1.7), 1e-5);
  EXPECT_NEAR(-0.0760311, sf_hydrogenicR_1_dx(3, 1.7), 1e-5);
  EXPECT_NEAR(-0.190078, sf_hydrogenicR_1_dy(3, 1.7), 1e-5);
  EXPECT_NEAR(0.0806774, sf_hydrogenicR_1_dx2(3, 1.7), 1e-5);
  EXPECT_NEAR(0.164734, sf_hydrogenicR_1_dxdy(3, 1.7), 1e-5);
  EXPECT_NEAR(0.570233, sf_hydrogenicR_1_dy2(3, 1.7), 1e-5);

  const char *name = "gsl_sf_hydrogenicR";
  AMPLFunction af = GetFunction(name);
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
  AMPLFunction af = GetFunction(name);
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
  AMPLFunction af = GetFunction("gsl_sf_coupling_3j");
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
  AMPLFunction af = GetFunction("gsl_sf_coupling_6j");
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
  AMPLFunction af = GetFunction("gsl_sf_coupling_9j");
  return;
  EXPECT_ALMOST_EQUAL_OR_NAN(value, af(ArgList(6, 16, 18, 8, 20, 14, 12, 10, 4)));
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
  EXPECT_EQ(0.5, sf_dilog_dx2(0));
  EXPECT_NEAR(0.00545177, sf_dilog_dx2(5), 1e-5);
}

TEST_F(GSLTest, EllInt) {
  TEST_FUNC(sf_ellint_Kcomp);
  EXPECT_EQ(M_PI_4, sf_ellint_Kcomp_dx2(0));
  EXPECT_NEAR(1.88651, sf_ellint_Kcomp_dx2(0.5), 1e-5);

  TEST_FUNC(sf_ellint_Ecomp);
  EXPECT_EQ(-M_PI_4, sf_ellint_Ecomp_dx2(0));
  EXPECT_NEAR(-1.08346, sf_ellint_Ecomp_dx2(0.5), 1e-5);

  ArgList Zero(0, 0);
  ArgList TestPt(0.5, 0.5);
  AMPLFunction f = GetFunction("gsl_sf_ellint_Pcomp");
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
