// GSL wrapper test.

#include <stdexcept>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_airy.h>
#include <gsl/gsl_sf_bessel.h>

#include "gtest/gtest.h"
#include "solvers/asl.h"
#include "tests/config.h"

namespace {

typedef double (*Func1)(double);
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

  static bool IsNaN(double x) {
    return x != x;
  }

  static bool IsNearOrNaN(double expected, double actual) {
    testing::internal::FloatingPoint<double> lhs(expected), rhs(actual);
    return lhs.AlmostEquals(rhs) ||
        (IsNaN(expected) && IsNaN(actual));
  }

  void TestFunc(const char *name, Func1 f, Func1 dx, Func1 dx2);
  void TestFunc(const char *name, Func2 f, Func2 dx, Func2 dy,
      Func2 dx2, Func2 dxdy, Func2 dy2);
  void TestFunc(const char *name, Func3 f, Func3 dx, Func3 dy, Func3 dz,
      Func3 dx2, Func3 dxdy, Func3 dxdz, Func3 dy2, Func3 dydz, Func3 dz2);
};

const double POINTS[] = {-5, 0, 5};
const size_t NUM_POINTS = sizeof(POINTS) / sizeof(*POINTS);

#define EXPECT_NEAR_OR_NAN(expected, actual) \
    EXPECT_PRED2(IsNearOrNaN, expected, actual)

void GSLTest::TestFunc(const char *name, Func1 f, Func1 dx, Func1 dx2) {
  func_info *fi = GetFunction(name);
  for (size_t i = 0; i != NUM_POINTS; ++i) {
    double x = POINTS[i];
    arglist args = {0};
    args.ra = &x;
    EXPECT_NEAR_OR_NAN(f(x), fi->funcp(&args)) << name << " at " << x;
    real deriv = 0;
    args.derivs = &deriv;
    fi->funcp(&args);
    EXPECT_NEAR_OR_NAN(dx(x), deriv) << name << " at " << x;
    real hes = 0;
    args.hes = &hes;
    fi->funcp(&args);
    EXPECT_NEAR_OR_NAN(dx2(x), hes) << name << " at " << x;
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
      EXPECT_NEAR_OR_NAN(f(x, y), fi->funcp(&args));
      real deriv[2] = {};
      args.derivs = deriv;
      fi->funcp(&args);
      EXPECT_NEAR_OR_NAN(dx(x, y), deriv[0]);
      EXPECT_NEAR_OR_NAN(dy(x, y), deriv[1]);
      real hes[3] = {};
      args.hes = hes;
      fi->funcp(&args);
      EXPECT_NEAR_OR_NAN(dx2(x, y),  hes[0]);
      EXPECT_NEAR_OR_NAN(dxdy(x, y), hes[1]);
      EXPECT_NEAR_OR_NAN(dy2(x, y),  hes[2]);
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
        arglist args = {0};
        real xy[] = {x, y, z};
        args.ra = xy;
        EXPECT_NEAR_OR_NAN(f(x, y, z), fi->funcp(&args));
        real deriv[3] = {};
        args.derivs = deriv;
        fi->funcp(&args);
        EXPECT_NEAR_OR_NAN(dx(x, y, z), deriv[0]);
        EXPECT_NEAR_OR_NAN(dy(x, y, z), deriv[1]);
        EXPECT_NEAR_OR_NAN(dz(x, y, z), deriv[2]);
        real hes[6] = {};
        args.hes = hes;
        fi->funcp(&args);
        EXPECT_NEAR_OR_NAN(dx2(x, y, z),  hes[0]);
        EXPECT_NEAR_OR_NAN(dxdy(x, y, z), hes[1]);
        EXPECT_NEAR_OR_NAN(dxdz(x, y, z), hes[2]);
        EXPECT_NEAR_OR_NAN(dy2(x, y, z),  hes[3]);
        EXPECT_NEAR_OR_NAN(dydz(x, y, z), hes[4]);
        EXPECT_NEAR_OR_NAN(dz2(x, y, z),  hes[5]);
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
  return x / sqrt(x * x + y * y + z * z);
}
double hypot3_dy(double x, double y, double z) {
  return y / sqrt(x * x + y * y + z * z);
}
double hypot3_dz(double x, double y, double z) {
  return z / sqrt(x * x + y * y + z * z);
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

#define TEST_FUNC(name) \
  TestFunc("gsl_" #name, gsl_##name, name##_dx, name##_dx2);

#define TEST_FUNC2(name) \
    TestFunc("gsl_" #name, gsl_##name, name##_dx, name##_dy, \
        name##_dx2, name##_dxdy, name##_dy2);

#define TEST_FUNC3(name) \
    TestFunc("gsl_" #name, gsl_##name, name##_dx, name##_dy, name##_dz, \
        name##_dx2, name##_dxdy, name##_dxdz, \
        name##_dy2, name##_dydz, name##_dz2);

TEST_F(GSLTest, Functions) {
  TEST_FUNC(log1p);
  TEST_FUNC(expm1);
  TEST_FUNC2(hypot);
  TEST_FUNC3(hypot3);
  TEST_FUNC(sf_airy_Ai);
  TEST_FUNC(sf_airy_Ai_scaled);
  ASSERT_NEAR(-0.00888609, sf_airy_Ai_scaled_dx(5), 1e-5);
  ASSERT_NEAR(0.00207512, sf_airy_Ai_scaled_dx2(5), 1e-5);
  TEST_FUNC(sf_airy_Bi);
  TEST_FUNC(sf_airy_Bi_scaled);
  ASSERT_NEAR(-0.0203063, sf_airy_Bi_scaled_dx(5), 1e-5);
  ASSERT_NEAR(0.00559418, sf_airy_Bi_scaled_dx2(5), 1e-5);
  TEST_FUNC(sf_bessel_J0);
}
}
