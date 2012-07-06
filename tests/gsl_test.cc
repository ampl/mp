// GSL wrapper test.

#include <stdexcept>
#include <vector>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_airy.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_clausen.h>
#include <gsl/gsl_sf_coulomb.h>
#include <gsl/gsl_sf_coupling.h>

#include "gtest/gtest.h"
#include "solvers/asl.h"
#include "tests/config.h"

using std::vector;

namespace {

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

  func_info *GetFunction(const char *name) const {
    func_info *fi = func_lookup(asl, name, 0);
    if (fi == nullptr)
      throw std::runtime_error(std::string("Function not found: ") + name);
    return fi;
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

  void TestFunc(const char *name, Func1 f, Func1 dx, Func1 dx2);
  void TestFunc(const char *name, FuncU f);
  void TestFunc(const char *name, FuncN1 f, FuncN1Result dx, FuncN1Result dx2);
  void TestFunc(const char *name, Func2 f, Func2 dx, Func2 dy,
      Func2 dx2, Func2 dxdy, Func2 dy2);
  void TestFunc(const char *name, Func3 f, Func3 dx, Func3 dy, Func3 dz,
      Func3 dx2, Func3 dxdy, Func3 dxdz, Func3 dy2, Func3 dydz, Func3 dz2);
};

class ArgList {
 private:
  TMInfo tmi;
  vector<real> ra;
  vector<real> derivs_;
  vector<real> hes_;
  arglist args;

  void init(ASL *asl, size_t size) {
    args.TMI = &tmi;
    args.AE = asl->i.ae;
    args.nr = args.n = size;
    ra.resize(size);
    args.ra = &ra[0];
  }

 public:
  ArgList(ASL *asl, real x) : tmi(), args() {
    init(asl, 1);
    ra[0] = x;
  }

  ArgList(ASL *asl, real x, real y) : tmi(), args() {
    init(asl, 2);
    ra[0] = x;
    ra[1] = y;
  }

  ArgList(ASL *asl, real x, real y, real z) : tmi(), args() {
    init(asl, 3);
    ra[0] = x;
    ra[1] = y;
    ra[2] = z;
  }

  ArgList &operator,(real arg) {
    ra.push_back(arg);
    args.nr = args.n = ra.size();
    args.ra = &ra[0];
    return *this;
  }

  ArgList &allocateDerivs() {
    derivs_.resize(ra.size());
    args.derivs = &derivs_[0];
    return *this;
  }

  void allocateHes() {
    hes_.resize(ra.size() * (ra.size() + 1) / 2);
    args.hes = &hes_[0];
  }

  void set(size_t index, real value) { ra[index] = value; }

  arglist *get() { return &args; }
  arglist *operator->() { return &args; }

  real deriv(size_t index = 0) const { return derivs_[index]; }
  real hes(size_t index = 0) const { return hes_[index]; }

  real Call(func_info *fi);
  void CallError(func_info *fi);
};

const double POINTS[] = {-5, -1.23, 0, 1.23, 5};
const size_t NUM_POINTS = sizeof(POINTS) / sizeof(*POINTS);

const double POINTS_FOR_N[] = {
    INT_MIN, INT_MIN + 1, -2, -1, 0, 1, 2, INT_MAX - 1, INT_MAX};
const size_t NUM_POINTS_FOR_N = sizeof(POINTS_FOR_N) / sizeof(*POINTS_FOR_N);

#define EXPECT_ALMOST_EQUAL_OR_NAN(expected, actual) \
    EXPECT_PRED2(AlmostEqualOrNaN, expected, actual)

real Call(func_info *fi, ArgList &args) {
  args->Errmsg = nullptr;
  real value = fi->funcp(args.get());
  EXPECT_TRUE(args->Errmsg == nullptr) << args->Errmsg;
  return value;
}

void CallError(func_info *fi, ArgList &args) {
  args->Errmsg = nullptr;
  fi->funcp(args.get());
  EXPECT_TRUE(args->Errmsg != nullptr);
  args->Errmsg = nullptr;
}

real ArgList::Call(func_info *fi) {
  return ::Call(fi, *this);
}

void ArgList::CallError(func_info *fi) {
  ::CallError(fi, *this);
}

void GSLTest::TestFunc(const char *name, Func1 f, Func1 dx, Func1 dx2) {
  func_info *fi = GetFunction(name);
  for (size_t i = 0; i != NUM_POINTS; ++i) {
    double x = POINTS[i];
    ArgList args(asl, x);
    EXPECT_ALMOST_EQUAL_OR_NAN(f(x), Call(fi, args)) << name << " at " << x;
    args.allocateDerivs();
    Call(fi, args);
    EXPECT_ALMOST_EQUAL_OR_NAN(dx(x), args.deriv()) << name << " at " << x;
    args.allocateHes();
    Call(fi, args);
    EXPECT_ALMOST_EQUAL_OR_NAN(dx2(x), args.hes()) << name << " at " << x;
  }
}

void GSLTest::TestFunc(const char *name, FuncU f) {
  func_info *fi = GetFunction(name);
  for (size_t i = 0; i != NUM_POINTS; ++i) {
    double x = POINTS[i];
    ArgList args(asl, x);
    if (static_cast<unsigned>(x) != x) {
      CallError(fi, args);
      continue;
    }
    EXPECT_ALMOST_EQUAL_OR_NAN(f(x), Call(fi, args))
      << name << " at " << x;
    args.allocateDerivs();
    CallError(fi, args);
  }
}

void GSLTest::TestFunc(const char *name, FuncN1 f,
    FuncN1Result dx, FuncN1Result dx2) {
  func_info *fi = GetFunction(name);
  for (size_t i = 0; i != NUM_POINTS_FOR_N; ++i) {
    for (size_t j = 0; j != NUM_POINTS; ++j) {
      int n = POINTS_FOR_N[i];
      bool skip = n < -100 || n > 100;
      double x = POINTS[j];
      ArgList args(asl, static_cast<real>(n), x);
      if (!skip)
        EXPECT_ALMOST_EQUAL_OR_NAN(f(n, x), Call(fi, args)) << name << " at " << x;

      args.allocateDerivs();
      CallError(fi, args);
      char dig = 0;
      args->dig = &dig;
      CallError(fi, args);
      dig = 1;
      Result r = dx(n, x);
      if (r.error) {
        CallError(fi, args);
      } else if (!skip) {
        Call(fi, args);
        EXPECT_ALMOST_EQUAL_OR_NAN(r.value, args.deriv(1))
          << name << " at " << n << ", " << x;
      }

      args.allocateHes();
      r = dx2(n, x);
      if (r.error) {
        CallError(fi, args);
      } else if (!skip) {
        Call(fi, args);
        EXPECT_ALMOST_EQUAL_OR_NAN(r.value, args.hes(2)) << name << " at " << x;
      }
    }
  }
}

void GSLTest::TestFunc(const char *name, Func2 f, Func2 dx, Func2 dy,
    Func2 dx2, Func2 dxdy, Func2 dy2) {
  func_info *fi = GetFunction(name);
  for (size_t i = 0; i != NUM_POINTS; ++i) {
    for (size_t j = 0; j != NUM_POINTS; ++j) {
      double x = POINTS[i], y = POINTS[j];
      ArgList args(asl, x, y);
      EXPECT_ALMOST_EQUAL_OR_NAN(f(x, y), Call(fi, args));
      args.allocateDerivs();
      char dig = 1;
      if (dx) {
        Call(fi, args);
        EXPECT_ALMOST_EQUAL_OR_NAN(dx(x, y), args.deriv(0));
      } else {
        CallError(fi, args);
        args->dig = &dig;
        Call(fi, args);
      }
      EXPECT_ALMOST_EQUAL_OR_NAN(dy(x, y), args.deriv(1));
      args.allocateHes();
      Call(fi, args);
      if (dx2)
        EXPECT_ALMOST_EQUAL_OR_NAN(dx2(x, y),  args.hes(0));
      if (dxdy)
        EXPECT_ALMOST_EQUAL_OR_NAN(dxdy(x, y), args.hes(1));
      EXPECT_ALMOST_EQUAL_OR_NAN(dy2(x, y),  args.hes(2));
    }
  }
}

void GSLTest::TestFunc(const char *name, Func3 f, Func3 dx, Func3 dy, Func3 dz,
    Func3 dx2, Func3 dxdy, Func3 dxdz, Func3 dy2, Func3 dydz, Func3 dz2) {
  func_info *fi = GetFunction(name);
  for (size_t i = 0; i != NUM_POINTS; ++i) {
    for (size_t j = 0; j != NUM_POINTS; ++j) {
      for (size_t k = 0; k != NUM_POINTS; ++k) {
        double x = POINTS[i], y = POINTS[j], z = POINTS[k];
        ArgList args(asl, x, y, z);
        EXPECT_ALMOST_EQUAL_OR_NAN(f(x, y, z), Call(fi, args)) << name;
        args.allocateDerivs();
        Call(fi, args);
        EXPECT_ALMOST_EQUAL_OR_NAN(dx(x, y, z), args.deriv(0)) << name;
        EXPECT_ALMOST_EQUAL_OR_NAN(dy(x, y, z), args.deriv(1)) << name;
        EXPECT_ALMOST_EQUAL_OR_NAN(dz(x, y, z), args.deriv(2)) << name;
        args.allocateHes();
        Call(fi, args);
        EXPECT_ALMOST_EQUAL_OR_NAN(dx2(x, y, z),  args.hes(0)) << name;
        EXPECT_ALMOST_EQUAL_OR_NAN(dxdy(x, y, z), args.hes(1)) << name;
        EXPECT_ALMOST_EQUAL_OR_NAN(dxdz(x, y, z), args.hes(2)) << name;
        EXPECT_ALMOST_EQUAL_OR_NAN(dy2(x, y, z),  args.hes(3)) << name;
        EXPECT_ALMOST_EQUAL_OR_NAN(dydz(x, y, z), args.hes(4)) << name;
        EXPECT_ALMOST_EQUAL_OR_NAN(dz2(x, y, z),  args.hes(5)) << name;
      }
    }
  }
}

double log1p_dx(double x) { return 1 / (x + 1); }
double log1p_dx2(double x) { return -1 / ((x + 1) * (x + 1)); }

double expm1_dx(double x) { return exp(x); }
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
double sf_airy_Ai_dx(double x) {
  return gsl_sf_airy_Ai_deriv(x, GSL_PREC_DOUBLE);
}
double sf_airy_Ai_dx2(double x) {
  return x * gsl_sf_airy_Ai(x);
}

double gsl_sf_airy_Bi(double x) {
  return ::gsl_sf_airy_Bi(x, GSL_PREC_DOUBLE);
}
double sf_airy_Bi_dx(double x) {
  return gsl_sf_airy_Bi_deriv(x, GSL_PREC_DOUBLE);
}
double sf_airy_Bi_dx2(double x) {
  return x * gsl_sf_airy_Bi(x);
}

double gsl_sf_airy_Ai_scaled(double x) {
  return ::gsl_sf_airy_Ai_scaled(x, GSL_PREC_DOUBLE);
}
double sf_airy_Ai_scaled_dx(double x) {
  return x > 0 ? (gsl_sf_airy_Ai_deriv_scaled(x, GSL_PREC_DOUBLE) +
      sqrt(x) * gsl_sf_airy_Ai_scaled(x)) :
      gsl_sf_airy_Ai_deriv(x, GSL_PREC_DOUBLE);
}
double sf_airy_Ai_scaled_dx2(double x) {
  return x > 0 ? (gsl_sf_airy_Ai_scaled(x) +
      4 * x * sf_airy_Ai_scaled_dx(x)) / (2 * sqrt(x)) :
      x * gsl_sf_airy_Ai(x);
}

double gsl_sf_airy_Bi_scaled(double x) {
  return ::gsl_sf_airy_Bi_scaled(x, GSL_PREC_DOUBLE);
}
double sf_airy_Bi_scaled_dx(double x) {
  return x > 0 ? (gsl_sf_airy_Bi_deriv_scaled(x, GSL_PREC_DOUBLE) -
      sqrt(x) * gsl_sf_airy_Bi_scaled(x)) :
      gsl_sf_airy_Bi_deriv(x, GSL_PREC_DOUBLE);
}
double sf_airy_Bi_scaled_dx2(double x) {
  return x > 0 ? -(gsl_sf_airy_Bi_scaled(x) +
      4 * x * sf_airy_Bi_scaled_dx(x)) / (2 * sqrt(x)) :
      x * gsl_sf_airy_Bi(x);
}

double sf_bessel_J0_dx(double x) { return -gsl_sf_bessel_J1(x); }
double sf_bessel_J0_dx2(double x) {
  return 0.5 * (gsl_sf_bessel_Jn(2, x) - gsl_sf_bessel_J0(x));
}

double sf_bessel_J1_dx(double x) {
  return 0.5 * (gsl_sf_bessel_J0(x) - gsl_sf_bessel_Jn(2, x));
}
double sf_bessel_J1_dx2(double x) {
  return 0.25 * (gsl_sf_bessel_Jn(3, x) - 3 * gsl_sf_bessel_J1(x));
}

Result sf_bessel_Jn_dx(int n, double x) {
  if (n == INT_MIN || n == INT_MAX)
    return Result(0, true);
  if (n <= INT_MIN + 1 || n >= INT_MAX - 1)
    return 42;
  return 0.5 * (gsl_sf_bessel_Jn(n - 1, x) - gsl_sf_bessel_Jn(n + 1, x));
}
Result sf_bessel_Jn_dx2(int n, double x) {
  if (n <= INT_MIN + 1 || n >= INT_MAX - 1)
    return Result(0, true);
  return 0.25 * (gsl_sf_bessel_Jn(n - 2, x) -
      2 * gsl_sf_bessel_Jn(n, x) + gsl_sf_bessel_Jn(n + 2, x));
}

double sf_bessel_Y0_dx(double x) { return -gsl_sf_bessel_Y1(x); }
double sf_bessel_Y0_dx2(double x) {
  return 0.5 * (gsl_sf_bessel_Yn(2, x) - gsl_sf_bessel_Y0(x));
}

double sf_bessel_Y1_dx(double x) {
  return 0.5 * (gsl_sf_bessel_Y0(x) - gsl_sf_bessel_Yn(2, x));
}
double sf_bessel_Y1_dx2(double x) {
  return 0.25 * (gsl_sf_bessel_Yn(3, x) - 3 * gsl_sf_bessel_Y1(x));
}

Result sf_bessel_Yn_dx(int n, double x) {
  if (n == INT_MIN || n == INT_MAX)
    return Result(0, true);
  if (n <= INT_MIN + 1 || n >= INT_MAX - 1)
    return 42;
  return 0.5 * (gsl_sf_bessel_Yn(n - 1, x) - gsl_sf_bessel_Yn(n + 1, x));
}
Result sf_bessel_Yn_dx2(int n, double x) {
  if (n <= INT_MIN + 1 || n >= INT_MAX - 1)
    return Result(0, true);
  return 0.25 * (gsl_sf_bessel_Yn(n - 2, x) -
      2 * gsl_sf_bessel_Yn(n, x) + gsl_sf_bessel_Yn(n + 2, x));
}

double sf_bessel_I0_dx(double x) { return gsl_sf_bessel_I1(x); }
double sf_bessel_I0_dx2(double x) {
  return 0.5 * (gsl_sf_bessel_In(2, x) + gsl_sf_bessel_I0(x));
}

double sf_bessel_I0_scaled_dx(double x) {
  return gsl_sf_bessel_I1_scaled(x) - x / abs(x) * gsl_sf_bessel_I0_scaled(x);
}
double sf_bessel_I0_scaled_dx2(double x) {
  return 1.5 * gsl_sf_bessel_I0_scaled(x) -
      2 * x / abs(x) * gsl_sf_bessel_I1_scaled(x) +
      0.5 * gsl_sf_bessel_In_scaled(2, x);
}

double sf_bessel_I1_dx(double x) {
  return 0.5 * (gsl_sf_bessel_I0(x) + gsl_sf_bessel_In(2, x));
}
double sf_bessel_I1_dx2(double x) {
  return 0.25 * (gsl_sf_bessel_In(3, x) + 3 * gsl_sf_bessel_I1(x));
}

Result sf_bessel_In_dx(int n, double x) {
  if (n == INT_MIN || n == INT_MAX)
    return Result(0, true);
  if (n <= INT_MIN + 1 || n >= INT_MAX - 1)
    return 42;
  return 0.5 * (gsl_sf_bessel_In(n - 1, x) + gsl_sf_bessel_In(n + 1, x));
}
Result sf_bessel_In_dx2(int n, double x) {
  if (n <= INT_MIN + 1 || n >= INT_MAX - 1)
    return Result(0, true);
  return 0.25 * (gsl_sf_bessel_In(n - 2, x) +
      2 * gsl_sf_bessel_In(n, x) + gsl_sf_bessel_In(n + 2, x));
}

double sf_bessel_I1_scaled_dx(double x) {
  return 0.5 * gsl_sf_bessel_I0_scaled(x) -
      x / abs(x) * gsl_sf_bessel_I1_scaled(x) +
      0.5 * gsl_sf_bessel_In_scaled(2, x);
}
double sf_bessel_I1_scaled_dx2(double x) {
  return -x / abs(x) * gsl_sf_bessel_I0_scaled(x) +
      1.75 * gsl_sf_bessel_I1_scaled(x) -
      x / abs(x) * gsl_sf_bessel_In_scaled(2, x) +
      0.25 * gsl_sf_bessel_In_scaled(3, x);
}

Result sf_bessel_In_scaled_dx(int n, double x) {
  if (n == INT_MIN || n == INT_MAX)
    return Result(0, true);
  if (n <= INT_MIN + 1 || n >= INT_MAX - 1)
    return 42;
  return 0.5 * (-(2 * x * gsl_sf_bessel_In_scaled(n, x)) / abs(x) +
      gsl_sf_bessel_In_scaled(n - 1, x) + gsl_sf_bessel_In_scaled(n + 1, x));
}
Result sf_bessel_In_scaled_dx2(int n, double x) {
  if (n <= INT_MIN + 1 || n >= INT_MAX - 1)
    return Result(0, true);
  return (abs(x) * (gsl_sf_bessel_In_scaled(n - 2, x) +
                    6 * gsl_sf_bessel_In_scaled(n, x) +
                    gsl_sf_bessel_In_scaled(n + 2, x)) -
      4 * x * (gsl_sf_bessel_In_scaled(n - 1, x) +
               gsl_sf_bessel_In_scaled(n + 1, x))) / (4 * abs(x));
}

double sf_bessel_K0_dx(double x) { return -gsl_sf_bessel_K1(x); }
double sf_bessel_K0_dx2(double x) {
  return 0.5 * (gsl_sf_bessel_Kn(2, x) + gsl_sf_bessel_K0(x));
}

double sf_bessel_K1_dx(double x) {
  return -0.5 * (gsl_sf_bessel_K0(x) + gsl_sf_bessel_Kn(2, x));
}
double sf_bessel_K1_dx2(double x) {
  return 0.25 * (gsl_sf_bessel_Kn(3, x) + 3 * gsl_sf_bessel_K1(x));
}

Result sf_bessel_Kn_dx(int n, double x) {
  if (n == INT_MIN || n == INT_MAX)
    return Result(0, true);
  if (n <= INT_MIN + 1 || n >= INT_MAX - 1)
    return 42;
  return -0.5 * (gsl_sf_bessel_Kn(n - 1, x) + gsl_sf_bessel_Kn(n + 1, x));
}
Result sf_bessel_Kn_dx2(int n, double x) {
  if (n <= INT_MIN + 1 || n >= INT_MAX - 1)
    return Result(0, true);
  return 0.25 * (gsl_sf_bessel_Kn(n - 2, x) +
      2 * gsl_sf_bessel_Kn(n, x) + gsl_sf_bessel_Kn(n + 2, x));
}

double sf_bessel_K0_scaled_dx(double x) {
  return gsl_sf_bessel_K0_scaled(x) - gsl_sf_bessel_K1_scaled(x);
}
double sf_bessel_K0_scaled_dx2(double x) {
  return 1.5 * gsl_sf_bessel_K0_scaled(x) - 2 * gsl_sf_bessel_K1_scaled(x) +
      0.5 * gsl_sf_bessel_Kn_scaled(2, x);
}

double sf_bessel_K1_scaled_dx(double x) {
  return -0.5 * gsl_sf_bessel_K0_scaled(x) + gsl_sf_bessel_K1_scaled(x) -
      0.5 * gsl_sf_bessel_Kn_scaled(2, x);
}
double sf_bessel_K1_scaled_dx2(double x) {
  return -gsl_sf_bessel_K0_scaled(x) + 1.75 * gsl_sf_bessel_K1_scaled(x) -
      gsl_sf_bessel_Kn_scaled(2, x) + 0.25 * gsl_sf_bessel_Kn_scaled(3, x);
}

Result sf_bessel_Kn_scaled_dx(int n, double x) {
  if (n == INT_MIN || n == INT_MAX)
    return Result(0, true);
  if (n <= INT_MIN + 1 || n >= INT_MAX - 1)
    return 42;
  return -0.5 * (gsl_sf_bessel_Kn_scaled(n - 1, x) -
      2 * gsl_sf_bessel_Kn_scaled(n, x) + gsl_sf_bessel_Kn_scaled(n + 1, x));
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

double sf_bessel_j0_dx(double x) { return (x * cos(x) - sin(x)) / (x * x); }
double sf_bessel_j0_dx2(double x) {
  return -((x * x - 2) * sin(x) + 2 * x * cos(x)) / (x * x * x);
}

double sf_bessel_j1_dx(double x) {
  return (sin(x) - 2 * gsl_sf_bessel_j1(x)) / x;
}
double sf_bessel_j1_dx2(double x) {
  return (x * (x * x - 6) * cos(x) -
      3 * (x * x - 2) * sin(x)) / (x * x * x * x);
}

double sf_bessel_j2_dx(double x) {
  return gsl_sf_bessel_j1(x) - (3 * gsl_sf_bessel_j2(x)) / x;
}
double sf_bessel_j2_dx2(double x) {
  double x_squared = x * x;
  return (x * (5 * x_squared - 36) * cos(x) +
      (x_squared * x_squared - 17 * x_squared + 36) * sin(x)) /
      (x_squared * x_squared * x);
}

Result sf_bessel_jl_dx(int n, double x) {
  if (n == INT_MAX)
    return Result(0, true);
  if (n <= INT_MIN + 1 || n >= INT_MAX - 1)
    return 42;
  return (n * gsl_sf_bessel_jl(n, x)) / x - gsl_sf_bessel_jl(n + 1, x);
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

double sf_bessel_y0_dx(double x) { return (x * sin(x) + cos(x)) / (x * x); }
double sf_bessel_y0_dx2(double x) {
  return ((x * x - 2) * cos(x) - 2 * x * sin(x)) / (x * x * x);
}

double sf_bessel_y1_dx(double x) {
  return -(2 * gsl_sf_bessel_y1(x) + cos(x)) / x;
}
double sf_bessel_y1_dx2(double x) {
  return (x * (x * x - 6) * sin(x) +
      3 * (x * x - 2) * cos(x)) / (x * x * x * x);
}

double sf_bessel_y2_dx(double x) {
  return gsl_sf_bessel_y1(x) - (3 * gsl_sf_bessel_y2(x)) / x;
}
double sf_bessel_y2_dx2(double x) {
  return ((36 - 5 * x * x) * gsl_sf_bessel_y1(x) -
      (x * x - 12) * cos(x)) / (x * x * x);
}

Result sf_bessel_yl_dx(int n, double x) {
  if (n == INT_MIN || n == INT_MAX)
    return Result(0, true);
  if (n <= INT_MIN + 1 || n >= INT_MAX - 1)
    return 42;
  return 0.5 * (gsl_sf_bessel_yl(n - 1, x) -
      gsl_sf_bessel_yl(n, x) / x - gsl_sf_bessel_yl(n + 1, x));
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

double sf_bessel_i0_scaled_dx(double x) {
  double i_minus1 = (exp(-abs(x)) * sqrt(1 / x) / sqrt(x)) * cosh(x);
  return 0.5 * (i_minus1 -
      ((1 + 2 * abs(x)) / x) * gsl_sf_bessel_i0_scaled(x) +
      gsl_sf_bessel_i1_scaled(x));
}
double sf_bessel_i0_scaled_dx2(double x) {
  double coef = -2 * (1 + 2 * abs(x)) / x;
  double hyp_coef = exp(-abs(x)) * sqrt(1 / x) / sqrt(x);
  double i_minus1 = hyp_coef * cosh(x);
  double i_minus2 = hyp_coef * sinh(x) - i_minus1 / x;
  return 0.25 * (
      i_minus2 +
      coef * i_minus1 +
      (3 + 6 * x * x + 4 * abs(x)) * gsl_sf_bessel_i0_scaled(x) / (x * x) +
      coef * gsl_sf_bessel_i1_scaled(x) +
      gsl_sf_bessel_i2_scaled(x));
}

double sf_bessel_i1_scaled_dx(double x) {
  return 0.5 * (gsl_sf_bessel_i0_scaled(x) -
      ((1 + 2 * abs(x)) / x) * gsl_sf_bessel_i1_scaled(x) +
      gsl_sf_bessel_i2_scaled(x));
}
double sf_bessel_i1_scaled_dx2(double x) {
  double coef = -2 * (1 + 2 * abs(x)) / x;
  double i_minus1 = (exp(-abs(x)) * sqrt(1 / x) * cosh(x)) / sqrt(x);
  return 0.25 * (
      i_minus1 +
      coef * gsl_sf_bessel_i0_scaled(x) +
      (3 + 6 * x * x + 4 * abs(x)) * gsl_sf_bessel_i1_scaled(x) / (x * x) +
      coef * gsl_sf_bessel_i2_scaled(x) +
      gsl_sf_bessel_il_scaled(3, x));
}

double sf_bessel_i2_scaled_dx(double x) {
  return 0.5 * (gsl_sf_bessel_i1_scaled(x) -
      ((1 + 2 * abs(x)) / x) * gsl_sf_bessel_i2_scaled(x) +
      gsl_sf_bessel_il_scaled(3, x));
}
double sf_bessel_i2_scaled_dx2(double x) {
  double coef = -2 * (1 + 2 * abs(x)) / x;
  return 0.25 * (
      gsl_sf_bessel_i0_scaled(x) +
      coef * gsl_sf_bessel_i1_scaled(x) +
      (3 + 6 * x * x + 4 * abs(x)) * gsl_sf_bessel_il_scaled(2, x) / (x * x) +
      coef * gsl_sf_bessel_il_scaled(3, x) +
      gsl_sf_bessel_il_scaled(4, x));
}

Result sf_bessel_il_scaled_dx(int n, double x) {
  if (n == INT_MIN || n == INT_MAX)
    return Result(0, true);
  if (n <= INT_MIN + 1 || n >= INT_MAX - 1)
    return 42;
  return 0.5 * (gsl_sf_bessel_il_scaled(n - 1, x) -
      ((1 + 2 * abs(x)) / x) * gsl_sf_bessel_il_scaled(n, x) +
      gsl_sf_bessel_il_scaled(n + 1, x));
}
Result sf_bessel_il_scaled_dx2(int n, double x) {
  if (n <= INT_MIN + 1 || n >= INT_MAX - 1)
    return Result(0, true);
  double coef = -2 * (1 + 2 * abs(x)) / x;
  return 0.25 * (
      gsl_sf_bessel_il_scaled(n - 2, x) +
      coef * gsl_sf_bessel_il_scaled(n - 1, x) +
      (3 + 4 * abs(x) + 6 * x * x) * gsl_sf_bessel_il_scaled(n, x) / (x * x) +
      coef * gsl_sf_bessel_il_scaled(n + 1, x) +
      gsl_sf_bessel_il_scaled(n + 2, x));
}

double sf_bessel_k0_scaled_dx(double x) {
  return -M_PI * sqrt(1 / x) / (2 * pow(x, 1.5));
}
double sf_bessel_k0_scaled_dx2(double x) {
  return M_PI * sqrt(1 / x) / pow(x, 2.5);
}

double sf_bessel_k1_scaled_dx(double x) {
  return -(M_PI * sqrt(1 / x) * (x + 2)) / (2 * pow(x, 2.5));
}
double sf_bessel_k1_scaled_dx2(double x) {
  return (M_PI * sqrt(1 / x) * (x + 3)) / pow(x, 3.5);
}

double sf_bessel_k2_scaled_dx(double x) {
  return -(M_PI * sqrt(1 / x) * (x + 3) * (x + 3)) / (2 * pow(x, 3.5));
}
double sf_bessel_k2_scaled_dx2(double x) {
  return (M_PI * sqrt(1 / x) * (x * x + 9 * x + 18)) / pow(x, 4.5);
}

Result sf_bessel_kl_scaled_dx(int n, double x) {
  if (n == INT_MIN || n == INT_MAX)
    return Result(0, true);
  if (n <= INT_MIN + 1 || n >= INT_MAX - 1)
    return 42;
  return -0.5 * (gsl_sf_bessel_kl_scaled(n - 1, x) +
      ((1 - 2 * x) / x) * gsl_sf_bessel_kl_scaled(n, x) +
      gsl_sf_bessel_kl_scaled(n + 1, x));
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
  return 0.5 * (-(2 * y * gsl_sf_bessel_Inu_scaled(x, y)) / abs(y) +
      gsl_sf_bessel_Inu_scaled(x - 1, y) + gsl_sf_bessel_Inu_scaled(x + 1, y));
}
double sf_bessel_Inu_scaled_dy2(double x, double y) {
  return (abs(y) * (gsl_sf_bessel_Inu_scaled(x - 2, y) +
                    6 * gsl_sf_bessel_Inu_scaled(x, y) +
                    gsl_sf_bessel_Inu_scaled(x + 2, y)) -
      4 * y * (gsl_sf_bessel_Inu_scaled(x - 1, y) +
               gsl_sf_bessel_Inu_scaled(x + 1, y))) / (4 * abs(y));
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

double sf_clausen_dx(double x) {
  return -log(2 * sin(0.5 * x));
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

#define TEST_FUNC(name) \
  TestFunc("gsl_" #name, gsl_##name, name##_dx, name##_dx2);

#define TEST_FUNC_U(name) \
  TestFunc("gsl_" #name, gsl_##name);

#define TEST_FUNC_N(name) \
  TestFunc("gsl_" #name, gsl_##name, name##_dx, name##_dx2);

#define TEST_FUNC2(name) \
    TestFunc("gsl_" #name, gsl_##name, name##_dx, name##_dy, \
        name##_dx2, name##_dxdy, name##_dy2);

#define TEST_FUNC3(name) \
    TestFunc("gsl_" #name, gsl_##name, name##_dx, name##_dy, name##_dz, \
        name##_dx2, name##_dxdy, name##_dxdz, \
        name##_dy2, name##_dydz, name##_dz2);

TEST_F(GSLTest, Elementary) {
  TEST_FUNC(log1p);
  TEST_FUNC(expm1);
  TEST_FUNC2(hypot);
  TEST_FUNC3(hypot3);
}

TEST_F(GSLTest, AiryA) {
  TEST_FUNC(sf_airy_Ai);
  TEST_FUNC(sf_airy_Ai_scaled);
  ASSERT_NEAR(-0.00888609, sf_airy_Ai_scaled_dx(5), 1e-5);
  ASSERT_NEAR(0.00207512, sf_airy_Ai_scaled_dx2(5), 1e-5);
}

TEST_F(GSLTest, AiryB) {
  TEST_FUNC(sf_airy_Bi);
  TEST_FUNC(sf_airy_Bi_scaled);
  ASSERT_NEAR(-0.0203063, sf_airy_Bi_scaled_dx(5), 1e-5);
  ASSERT_NEAR(0.00559418, sf_airy_Bi_scaled_dx2(5), 1e-5);
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
  ASSERT_NEAR(-0.172334, sf_bessel_Jn_dx(3, 5).value, 1e-5);
  ASSERT_NEAR(-0.199025, sf_bessel_Jn_dx2(3, 5).value, 1e-5);
}

TEST_F(GSLTest, BesselY) {
  TEST_FUNC(sf_bessel_Y0);
  TEST_FUNC(sf_bessel_Y1);
  TEST_FUNC_N(sf_bessel_Yn);
  ASSERT_NEAR(0.279903, sf_bessel_Yn_dx(3, 5).value, 1e-5);
  ASSERT_NEAR(-0.149592, sf_bessel_Yn_dx2(3, 5).value, 1e-5);
}

TEST_F(GSLTest, BesselI) {
  TEST_FUNC(sf_bessel_I0);
  TEST_FUNC(sf_bessel_I1);
  TEST_FUNC_N(sf_bessel_In);
  ASSERT_NEAR(11.30692, sf_bessel_In_dx(3, 5).value, 1e-5);
  ASSERT_NEAR(11.78897, sf_bessel_In_dx2(3, 5).value, 1e-5);

  TEST_FUNC(sf_bessel_I0_scaled);
  ASSERT_NEAR(-0.0195685, sf_bessel_I0_scaled_dx(5), 1e-5);
  ASSERT_NEAR(0.00634264, sf_bessel_I0_scaled_dx2(5), 1e-5);

  TEST_FUNC(sf_bessel_I1_scaled);
  ASSERT_NEAR(-0.0132259, sf_bessel_I1_scaled_dx(5), 1e-5);
  ASSERT_NEAR(0.00286143, sf_bessel_I1_scaled_dx2(5), 1e-5);

  TEST_FUNC_N(sf_bessel_In_scaled);
  ASSERT_NEAR(0.00657472, sf_bessel_In_scaled_dx(3, 5).value, 1e-5);
  ASSERT_NEAR(-0.00332666, sf_bessel_In_scaled_dx2(3, 5).value, 1e-5);
}

TEST_F(GSLTest, BesselK) {
  TEST_FUNC(sf_bessel_K0);
  TEST_FUNC(sf_bessel_K1);
  TEST_FUNC_N(sf_bessel_Kn);
  ASSERT_NEAR(-0.010284, sf_bessel_Kn_dx(3, 5).value, 1e-5);
  ASSERT_NEAR(0.0133336, sf_bessel_Kn_dx2(3, 5).value, 1e-5);

  TEST_FUNC(sf_bessel_K0_scaled);
  ASSERT_NEAR(-0.0524663, sf_bessel_K0_scaled_dx(5), 1e-5);
  ASSERT_NEAR(0.0151222, sf_bessel_K0_scaled_dx2(5), 1e-5);

  TEST_FUNC(sf_bessel_K1_scaled);
  ASSERT_NEAR(-0.0675885, sf_bessel_K1_scaled_dx(5), 1e-5);
  ASSERT_NEAR(0.0224065, sf_bessel_K1_scaled_dx2(5), 1e-5);

  TEST_FUNC_N(sf_bessel_Kn_scaled);
  ASSERT_NEAR(-0.295674, sf_bessel_Kn_scaled_dx(3, 5).value, 1e-5);
  ASSERT_NEAR(0.156927, sf_bessel_Kn_scaled_dx2(3, 5).value, 1e-5);
}

TEST_F(GSLTest, Besselj) {
  TEST_FUNC(sf_bessel_j0);
  ASSERT_NEAR(0.0950894, sf_bessel_j0_dx(5), 1e-5);
  ASSERT_NEAR(0.153749, sf_bessel_j0_dx2(5), 1e-5);

  TEST_FUNC(sf_bessel_j1);
  ASSERT_NEAR(-0.153749, sf_bessel_j1_dx(5), 1e-5);
  ASSERT_NEAR(0.148982, sf_bessel_j1_dx2(5), 1e-5);

  TEST_FUNC(sf_bessel_j2);
  ASSERT_NEAR(-0.175928, sf_bessel_j2_dx(5), 1e-5);
  ASSERT_NEAR(-0.0320245, sf_bessel_j2_dx2(5), 1e-5);

  TEST_FUNC_N(sf_bessel_jl);
  ASSERT_NEAR(-0.0491253, sf_bessel_jl_dx(3, 5).value, 1e-5);
  ASSERT_NEAR(-0.0998566, sf_bessel_jl_dx2(3, 5).value, 1e-5);
}

TEST_F(GSLTest, Bessely) {
  TEST_FUNC(sf_bessel_y0);
  ASSERT_NEAR(-0.180438, sf_bessel_y0_dx(5), 1e-5);
  ASSERT_NEAR(0.128908, sf_bessel_y0_dx2(5), 1e-5);

  TEST_FUNC(sf_bessel_y1);
  ASSERT_NEAR(-0.128908, sf_bessel_y1_dx(5), 1e-5);
  ASSERT_NEAR(-0.11444, sf_bessel_y1_dx2(5), 1e-5);

  TEST_FUNC(sf_bessel_y2);
  ASSERT_NEAR(0.0814411, sf_bessel_y2_dx(5), 1e-5);
  ASSERT_NEAR(-0.157973, sf_bessel_y2_dx2(5), 1e-5);

  TEST_FUNC_N(sf_bessel_yl);
  ASSERT_NEAR(0.17735, sf_bessel_yl_dx(3, 5).value, 1e-5);
  ASSERT_NEAR(-0.0629096, sf_bessel_yl_dx2(3, 5).value, 1e-5);
}

TEST_F(GSLTest, Besseli) {
  TEST_FUNC(sf_bessel_i0_scaled);
  ASSERT_NEAR(0.0999955, gsl_sf_bessel_i0_scaled(5), 1e-5);
  ASSERT_NEAR(-0.01999, sf_bessel_i0_scaled_dx(5), 1e-5);
  ASSERT_NEAR(0.00797784, sf_bessel_i0_scaled_dx2(5), 1e-5);

  TEST_FUNC(sf_bessel_i1_scaled);
  ASSERT_NEAR(0.0800054, gsl_sf_bessel_i1_scaled(5), 1e-5);
  ASSERT_NEAR(-0.0120122, sf_bessel_i1_scaled_dx(5), 1e-5);
  ASSERT_NEAR(0.00322746, sf_bessel_i1_scaled_dx2(5), 1e-5);

  TEST_FUNC(sf_bessel_i2_scaled);
  ASSERT_NEAR(0.0519922, gsl_sf_bessel_i2_scaled(5), 1e-5);
  ASSERT_NEAR(-0.00318206, sf_bessel_i2_scaled_dx(5), 1e-5);
  ASSERT_NEAR(-0.000681812, sf_bessel_i2_scaled_dx2(5), 1e-5);

  TEST_FUNC_N(sf_bessel_il_scaled);
  ASSERT_NEAR(0.0280133, gsl_sf_bessel_il_scaled(3, 5), 1e-5);
  ASSERT_NEAR(0.00156833, sf_bessel_il_scaled_dx(3, 5).value, 1e-5);
  ASSERT_NEAR(-0.00152293, sf_bessel_il_scaled_dx2(3, 5).value, 1e-5);
}

TEST_F(GSLTest, Besselk) {
  TEST_FUNC(sf_bessel_k0_scaled);
  ASSERT_NEAR(0.314159, gsl_sf_bessel_k0_scaled(5), 1e-5);
  ASSERT_NEAR(-0.0628319, sf_bessel_k0_scaled_dx(5), 1e-5);
  ASSERT_NEAR(0.0251327, sf_bessel_k0_scaled_dx2(5), 1e-5);

  TEST_FUNC(sf_bessel_k1_scaled);
  ASSERT_NEAR(0.376991, gsl_sf_bessel_k1_scaled(5), 1e-5);
  ASSERT_NEAR(-0.0879646, sf_bessel_k1_scaled_dx(5), 1e-5);
  ASSERT_NEAR(0.0402124, sf_bessel_k1_scaled_dx2(5), 1e-5);

  TEST_FUNC(sf_bessel_k2_scaled);
  ASSERT_NEAR(0.540354, gsl_sf_bessel_k2_scaled(5), 1e-5);
  ASSERT_NEAR(-0.16085, sf_bessel_k2_scaled_dx(5), 1e-5);
  ASSERT_NEAR(0.0884672, sf_bessel_k2_scaled_dx2(5), 1e-5);

  TEST_FUNC_N(sf_bessel_kl_scaled);
  ASSERT_NEAR(0.917345, gsl_sf_bessel_kl_scaled(3, 5), 1e-5);
  ASSERT_NEAR(-0.356885, sf_bessel_kl_scaled_dx(3, 5).value, 1e-5);
  ASSERT_NEAR(0.236248, sf_bessel_kl_scaled_dx2(3, 5).value, 1e-5);
}

TEST_F(GSLTest, BesselFractionalOrder) {
  TEST_FUNC2(sf_bessel_Jnu);
  ASSERT_NEAR(0.410029, gsl_sf_bessel_Jnu(3.5, 5), 1e-5);
  ASSERT_NEAR(-0.0466428, sf_bessel_Jnu_dy(3.5, 5), 1e-5);
  ASSERT_NEAR(-0.199786, sf_bessel_Jnu_dy2(3.5, 5), 1e-5);

  TEST_FUNC2(sf_bessel_Ynu);
  ASSERT_NEAR(-0.0275521, gsl_sf_bessel_Ynu(3.5, 5), 1e-5);
  ASSERT_NEAR(0.313659, sf_bessel_Ynu_dy(3.5, 5), 1e-5);
  ASSERT_NEAR(-0.0486802, sf_bessel_Ynu_dy2(3.5, 5), 1e-5);

  TEST_FUNC2(sf_bessel_Inu);
  ASSERT_NEAR(7.417560126111555, gsl_sf_bessel_Inu(3.5, 5), 1e-5);
  ASSERT_NEAR(8.574590050404494, sf_bessel_Inu_dy(3.5, 5), 1e-5);
  ASSERT_NEAR(9.337246577825318, sf_bessel_Inu_dy2(3.5, 5), 1e-5);

  TEST_FUNC2(sf_bessel_Inu_scaled);
  ASSERT_NEAR(0.04997912699226937, gsl_sf_bessel_Inu_scaled(3.5, 5), 1e-5);
  ASSERT_NEAR(0.00779600630624169, sf_bessel_Inu_scaled_dy(3.5, 5), 1e-5);
  ASSERT_NEAR(-0.0026572670459736, sf_bessel_Inu_scaled_dy2(3.5, 5), 1e-5);

  TEST_FUNC2(sf_bessel_Knu);
  ASSERT_NEAR(0.011027711053957217, gsl_sf_bessel_Knu(3.5, 5), 1e-5);
  ASSERT_NEAR(-0.014215172742155810, sf_bessel_Knu_dy(3.5, 5), 1e-5);
  ASSERT_NEAR(0.01927432401882742, sf_bessel_Knu_dy2(3.5, 5), 1e-5);

  TEST_FUNC2(sf_bessel_lnKnu);
  ASSERT_NEAR(-4.5073439872921324, gsl_sf_bessel_lnKnu(3.5, 5), 1e-5);
  ASSERT_NEAR(-1.289041095890411, sf_bessel_lnKnu_dy(3.5, 5), 1e-5);
  ASSERT_NEAR(0.08618127228373, sf_bessel_lnKnu_dy2(3.5, 5), 1e-5);

  TEST_FUNC2(sf_bessel_Knu_scaled);
  ASSERT_NEAR(1.6366574351881952, gsl_sf_bessel_Knu_scaled(3.5, 5), 1e-5);
  ASSERT_NEAR(-0.4730612586639852, sf_bessel_Knu_scaled_dy(3.5, 5), 1e-5);
  ASSERT_NEAR(0.2777833646846813, sf_bessel_Knu_scaled_dy2(3.5, 5), 1e-5);
}

TEST_F(GSLTest, BesselZero) {
  TEST_FUNC_U(sf_bessel_zero_J0);
  TEST_FUNC_U(sf_bessel_zero_J1);

  const char *name = "gsl_sf_bessel_zero_Jnu";
  func_info *fi = GetFunction(name);
  for (size_t i = 0; i != NUM_POINTS; ++i) {
    for (size_t j = 0; j != NUM_POINTS; ++j) {
      double nu = POINTS[i], x = POINTS[j];
      ArgList args(asl, nu, x);
      if (static_cast<unsigned>(x) != x) {
        CallError(fi, args);
        continue;
      }
      EXPECT_ALMOST_EQUAL_OR_NAN(gsl_sf_bessel_zero_Jnu(nu, x), Call(fi, args))
        << name << " at " << x;
      args.allocateDerivs();
      CallError(fi, args);
    }
  }
}

TEST_F(GSLTest, Clausen) {
  TEST_FUNC(sf_clausen);
}

TEST_F(GSLTest, Hydrogenic) {
  TEST_FUNC2(sf_hydrogenicR_1);
  ASSERT_NEAR(0.0633592, gsl_sf_hydrogenicR_1(3, 1.7), 1e-5);
  ASSERT_NEAR(-0.0760311, sf_hydrogenicR_1_dx(3, 1.7), 1e-5);
  ASSERT_NEAR(-0.190078, sf_hydrogenicR_1_dy(3, 1.7), 1e-5);
  ASSERT_NEAR(0.0806774, sf_hydrogenicR_1_dx2(3, 1.7), 1e-5);
  ASSERT_NEAR(0.164734, sf_hydrogenicR_1_dxdy(3, 1.7), 1e-5);
  ASSERT_NEAR(0.570233, sf_hydrogenicR_1_dy2(3, 1.7), 1e-5);

  const char *name = "gsl_sf_hydrogenicR";
  func_info *fi = GetFunction(name);
  (ArgList(asl, 0, 0, 0), 0).Call(fi);
  (ArgList(asl, 0.5, 0, 0), 0).CallError(fi);
  (ArgList(asl, 0, 0.5, 0), 0).CallError(fi);
  for (size_t i = 0; i != NUM_POINTS_FOR_N; ++i) {
    for (size_t i = 0; i != NUM_POINTS_FOR_N; ++i) {
      for (size_t k = 0; k != NUM_POINTS; ++k) {
        for (size_t l = 0; l != NUM_POINTS; ++l) {
          int n = POINTS_FOR_N[i], ll = POINTS_FOR_N[l];
          if (n < -1000 || n > 1000) continue;
          double Z = POINTS[k], r = POINTS[l];
          ArgList args(asl, n, ll, Z);
          args, r;
          EXPECT_ALMOST_EQUAL_OR_NAN(
              gsl_sf_hydrogenicR(n, ll, Z, r), Call(fi, args))
            << name << " at " << n;
          args.allocateDerivs();
          CallError(fi, args);
        }
      }
    }
  }
}

TEST_F(GSLTest, Coulomb) {
  const char *name = "gsl_sf_coulomb_CL";
  func_info *fi = GetFunction(name);
  for (size_t i = 0; i != NUM_POINTS; ++i) {
    for (size_t j = 0; j != NUM_POINTS; ++j) {
      double x = POINTS[i], y = POINTS[j];
      ArgList args(asl, x, y);
      gsl_sf_result result = {};
      double value = gsl_sf_coulomb_CL_e(x, y, &result) ? GSL_NAN : result.val;
      EXPECT_ALMOST_EQUAL_OR_NAN(value, Call(fi, args));
      args.allocateDerivs();
      CallError(fi, args);
    }
  }
}

TEST_F(GSLTest, Coupling) {
  double value = gsl_sf_coupling_3j(12, 8, 4, 0, 0, 0);
  ASSERT_NEAR(0.186989, value, 1e-5);
  const char *name = "gsl_sf_coupling_3j";
  func_info *fi = GetFunction(name);
  EXPECT_ALMOST_EQUAL_OR_NAN(value,
      (ArgList(asl, 12, 8, 4), 0, 0, 0).Call(fi));
  (ArgList(asl, 0, 0, 0), 0, 0, 0).Call(fi);
  (ArgList(asl, 0.5, 0, 0), 0, 0, 0).CallError(fi);
  (ArgList(asl, 0, 0.5, 0), 0, 0, 0).CallError(fi);
  (ArgList(asl, 0, 0, 0.5), 0, 0, 0).CallError(fi);
  (ArgList(asl, 0, 0, 0), 0.5, 0, 0).CallError(fi);
  (ArgList(asl, 0, 0, 0), 0, 0.5, 0).CallError(fi);
  (ArgList(asl, 0, 0, 0), 0, 0, 0.5).CallError(fi);
  (ArgList(asl, 12, 8, 4), 0, 0, 0).allocateDerivs().CallError(fi);
}
}
