// GSL wrapper test.

#include <stdexcept>

#include "gtest/gtest.h"
#include "solvers/asl.h"
#include "tests/config.h"

namespace {

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

  static real Call(func_info *f, real x) {
    arglist args = {0};
    args.ra = &x;
    return f->funcp(&args);
  }

  static real Deriv(func_info *f, real x) {
    arglist args = {0};
    args.ra = &x;
    real deriv = 0;
    args.derivs = &deriv;
    f->funcp(&args);
    return deriv;
  }

  static real Deriv2(func_info *f, real x) {
    arglist args = {0};
    args.ra = &x;
    real deriv = 0;
    args.derivs = &deriv;
    real hes = 0;
    args.hes = &hes;
    f->funcp(&args);
    return hes;
  }

  static bool IsNaN(double x) {
    return x != x;
  }
};

TEST_F(GSLTest, gsl_log1p) {
  func_info *f = GetFunction("gsl_log1p");
  EXPECT_NEAR(0, Call(f, 0), 1e-5);
  EXPECT_NEAR(1.7917594, Call(f, 5), 1e-5);
  EXPECT_TRUE(IsNaN(Call(f, -5)));

  EXPECT_NEAR(1, Deriv(f, 0), 1e-5);
  EXPECT_NEAR(0.1666666, Deriv(f, 5), 1e-5);
  EXPECT_NEAR(-0.25, Deriv(f, -5), 1e-5);

  EXPECT_NEAR(-1, Deriv2(f, 0), 1e-5);
  EXPECT_NEAR(-0.0277777, Deriv2(f,  5), 1e-5);
  EXPECT_NEAR(-0.0625,    Deriv2(f, -5), 1e-5);
}

TEST_F(GSLTest, gsl_sf_bessel_J0) {
  func_info *f = GetFunction("gsl_sf_bessel_J0");
  EXPECT_NEAR(1, Call(f, 0), 1e-5);
  EXPECT_NEAR(-0.1775967, Call(f,  5), 1e-5);
  EXPECT_NEAR(-0.1775967, Call(f, -5), 1e-5);

  EXPECT_NEAR(0, Deriv(f, 0), 1e-5);
  EXPECT_NEAR( 0.3275791, Deriv(f,  5), 1e-5);
  EXPECT_NEAR(-0.3275791, Deriv(f, -5), 1e-5);

  EXPECT_NEAR(-0.5, Deriv2(f, 0), 1e-5);
  EXPECT_NEAR(0.1120809, Deriv2(f,  5), 1e-5);
  EXPECT_NEAR(0.1120809, Deriv2(f, -5), 1e-5);
}
}
