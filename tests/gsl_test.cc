// GSL wrapper test.

#include <stdexcept>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>

#include "gtest/gtest.h"
#include "solvers/asl.h"
#include "tests/config.h"

namespace {

typedef double (*Func1)(double);
typedef double (*Func2)(double, double);

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

  static bool IsNaN(double x) {
    return x != x;
  }

  static void ExpectNearOrNaN(double expected, double actual,
      const char *where) {
    if (IsNaN(expected))
      EXPECT_TRUE(IsNaN(actual)) << "in " << where;
    else
      EXPECT_NEAR(expected, actual, 1e-5) << "in " << where;
  }

  void TestFunc(const char *name, Func1 f, Func1 dx, Func1 dx2);
  void TestFunc(const char *name, Func2 f, Func2 dx, Func2 dy,
      Func2 dx2, Func2 dxdy, Func2 dy2);
};

const double POINTS[] = {-5, 0, 5};
const size_t NUM_POINTS = sizeof(POINTS) / sizeof(*POINTS);

void GSLTest::TestFunc(const char *name, Func1 f, Func1 dx, Func1 dx2) {
  func_info *fi = GetFunction(name);
  for (size_t i = 0; i != NUM_POINTS; ++i) {
    double x = POINTS[i];
    arglist args = {0};
    args.ra = &x;
    ExpectNearOrNaN(f(x), fi->funcp(&args), name);
    real deriv = 0;
    args.derivs = &deriv;
    fi->funcp(&args);
    ExpectNearOrNaN(dx(x), deriv, name);
    real hes = 0;
    args.hes = &hes;
    fi->funcp(&args);
    ExpectNearOrNaN(dx2(x), hes, name);
  }
}

void GSLTest::TestFunc(const char *name, Func2 f, Func2 dx, Func2 dy,
    Func2 dx2, Func2 dxdy, Func2 dy2) {
  func_info *fi = GetFunction(name);
  for (size_t i = 0; i != NUM_POINTS; ++i) {
    for (size_t j = 0; j != NUM_POINTS; ++j) {
      double x = POINTS[i], y = POINTS[j];
      arglist args = {0};
      real xy[] = {x, y};
      args.ra = xy;
      ExpectNearOrNaN(f(x, y), fi->funcp(&args), name);
      real deriv[2] = {};
      args.derivs = deriv;
      fi->funcp(&args);
      ExpectNearOrNaN(dx(x, y), deriv[0], name);
      ExpectNearOrNaN(dy(x, y), deriv[1], name);
      real hes[3] = {};
      args.hes = hes;
      fi->funcp(&args);
      ExpectNearOrNaN(dx2(x, y),  hes[0], name);
      ExpectNearOrNaN(dxdy(x, y), hes[1], name);
      ExpectNearOrNaN(dy2(x, y),  hes[2], name);
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

double sf_bessel_J0_dx(double x) { return -gsl_sf_bessel_J1(x); }
double sf_bessel_J0_dx2(double x) {
  return 0.5 * (gsl_sf_bessel_Jn(2, x) - gsl_sf_bessel_J0(x));
}

#define TEST_FUNC(name) \
  TestFunc("gsl_" #name, gsl_##name, name##_dx, name##_dx2);

#define TEST_FUNC2(name) \
    TestFunc("gsl_" #name, gsl_##name, name##_dx, name##_dy, \
        name##_dx2, name##_dxdy, name##_dy2);

TEST_F(GSLTest, Functions) {
  TEST_FUNC(log1p);
  TEST_FUNC(expm1);
  TEST_FUNC2(hypot);
  TEST_FUNC(sf_bessel_J0);
}
}
