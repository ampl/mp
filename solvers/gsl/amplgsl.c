// AMPL bindings for the GNU Scientific Library.

#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_airy.h>
#include <gsl/gsl_sf_bessel.h>
#include "solvers/funcadd.h"

static real amplgsl_log1p(arglist *al) {
  real x = al->ra[0];
  if (al->derivs) {
    real deriv = *al->derivs = 1 / (x + 1);
    if (al->hes)
      *al->hes = -deriv * deriv;
  }
  return gsl_log1p(x);
}

static real amplgsl_expm1(arglist *al) {
  real x = al->ra[0];
  if (al->derivs) {
    real deriv = *al->derivs = exp(x);
    if (al->hes)
      *al->hes = deriv;
  }
  return gsl_expm1(x);
}

static real amplgsl_hypot(arglist *al) {
  real x = al->ra[0];
  real y = al->ra[1];
  real hypot = gsl_hypot(x, y);
  if (al->derivs) {
    real *derivs = al->derivs;
    derivs[0] = x / hypot;
    derivs[1] = y / hypot;
    if (al->hes) {
      real *hes = al->hes;
      hes[0] =  derivs[1] * derivs[1] / hypot;
      hes[1] = -derivs[0] * derivs[1] / hypot;
      hes[2] =  derivs[0] * derivs[0] / hypot;
    }
  }
  return hypot;
}

static real amplgsl_hypot3(arglist *al) {
  real x = al->ra[0];
  real y = al->ra[1];
  real z = al->ra[2];
  real hypot = gsl_hypot3(x, y, z);
  if (al->derivs) {
    real *derivs = al->derivs;
    derivs[0] = x / hypot;
    derivs[1] = y / hypot;
    derivs[2] = z / hypot;
    if (al->hes) {
      real *hes = al->hes;
      real dx2 = derivs[0] * derivs[0];
      real dy2 = derivs[1] * derivs[1];
      real dz2 = derivs[2] * derivs[2];
      hes[0] =  (dy2 + dz2) / hypot;
      hes[1] = -derivs[0] * derivs[1] / hypot;
      hes[2] = -derivs[0] * derivs[2] / hypot;
      hes[3] =  (dx2 + dz2) / hypot;
      hes[4] = -derivs[1] * derivs[2] / hypot;
      hes[5] =  (dx2 + dy2) / hypot;
    }
  }
  return hypot;
}

static real amplgsl_sf_airy_Ai(arglist *al) {
  real x = al->ra[0];
  real value = gsl_sf_airy_Ai(x, GSL_PREC_DOUBLE);
  if (al->derivs) {
    *al->derivs = gsl_sf_airy_Ai_deriv(x, GSL_PREC_DOUBLE);
    if (al->hes)
      *al->hes = x * value;
  }
  return value;
}

static real amplgsl_sf_airy_Bi(arglist *al) {
  real x = al->ra[0];
  real value = gsl_sf_airy_Bi(x, GSL_PREC_DOUBLE);
  if (al->derivs) {
    *al->derivs = gsl_sf_airy_Bi_deriv(x, GSL_PREC_DOUBLE);
    if (al->hes)
      *al->hes = x * value;
  }
  return value;
}

static real amplgsl_sf_airy_Ai_scaled(arglist *al) {
  real x = al->ra[0];
  real value = gsl_sf_airy_Ai_scaled(x, GSL_PREC_DOUBLE);
  if (al->derivs) {
    if (x > 0) {
      real sqrtx = sqrt(x);
      *al->derivs = gsl_sf_airy_Ai_deriv_scaled(x, GSL_PREC_DOUBLE) +
          sqrtx * gsl_sf_airy_Ai_scaled(x, GSL_PREC_DOUBLE);
      if (al->hes)
        *al->hes = (value + 4 * x * *al->derivs) / (2 * sqrtx);
    } else {
      *al->derivs = gsl_sf_airy_Ai_deriv(x, GSL_PREC_DOUBLE);
      if (al->hes)
        *al->hes = x * value;
    }
  }
  return value;
}

static real amplgsl_sf_airy_Bi_scaled(arglist *al) {
  real x = al->ra[0];
  real value = gsl_sf_airy_Bi_scaled(x, GSL_PREC_DOUBLE);
  if (al->derivs) {
    if (x > 0) {
      real sqrtx = sqrt(x);
      *al->derivs = gsl_sf_airy_Bi_deriv_scaled(x, GSL_PREC_DOUBLE) -
          sqrtx * gsl_sf_airy_Bi_scaled(x, GSL_PREC_DOUBLE);
      if (al->hes)
        *al->hes = -(value + 4 * x * *al->derivs) / (2 * sqrtx);
    } else {
      *al->derivs = gsl_sf_airy_Bi_deriv(x, GSL_PREC_DOUBLE);
      if (al->hes)
        *al->hes = x * value;
    }
  }
  return value;
}

static real amplgsl_sf_bessel_J0(arglist *al) {
  real x = al->ra[0];
  real j0 = gsl_sf_bessel_J0(x);
  if (al->derivs) {
    *al->derivs = -gsl_sf_bessel_J1(x);
    if (al->hes)
      *al->hes = 0.5 * (gsl_sf_bessel_Jn(2, x) - j0);
  }
  return j0;
}

static real amplgsl_sf_bessel_J1(arglist *al) {
  real x = al->ra[0];
  real j1 = gsl_sf_bessel_J1(x);
  if (al->derivs) {
    *al->derivs = 0.5 * (gsl_sf_bessel_J0(x) - gsl_sf_bessel_Jn(2, x));
    if (al->hes)
      *al->hes = 0.25 * (gsl_sf_bessel_Jn(3, x) - 3 * j1);
  }
  return j1;
}

static real amplgsl_sf_bessel_Y0(arglist *al) {
  real x = al->ra[0];
  real y0 = gsl_sf_bessel_Y0(x);
  if (al->derivs) {
    *al->derivs = -gsl_sf_bessel_Y1(x);
    if (al->hes)
      *al->hes = 0.5 * (gsl_sf_bessel_Yn(2, x) - y0);
  }
  return y0;
}

static real amplgsl_sf_bessel_Y1(arglist *al) {
  real x = al->ra[0];
  real y1 = gsl_sf_bessel_Y1(x);
  if (al->derivs) {
    *al->derivs = 0.5 * (gsl_sf_bessel_Y0(x) - gsl_sf_bessel_Yn(2, x));
    if (al->hes)
      *al->hes = 0.25 * (gsl_sf_bessel_Yn(3, x) - 3 * y1);
  }
  return y1;
}

static real amplgsl_sf_bessel_I0(arglist *al) {
  real x = al->ra[0];
  real i0 = gsl_sf_bessel_I0(x);
  if (al->derivs) {
    *al->derivs = gsl_sf_bessel_I1(x);
    if (al->hes)
      *al->hes = 0.5 * (gsl_sf_bessel_In(2, x) + i0);
  }
  return i0;
}

static real amplgsl_sf_bessel_I0_scaled(arglist *al) {
  real x = al->ra[0];
  real i0 = gsl_sf_bessel_I0_scaled(x);
  if (al->derivs) {
    real x_div_absx = x / abs(x);
    real i1 = gsl_sf_bessel_I1_scaled(x);
    *al->derivs = i1 - x_div_absx * i0;
    if (al->hes) {
      *al->hes = 1.5 * i0 - 2 * x_div_absx * i1 +
          0.5 * gsl_sf_bessel_In_scaled(2, x);
    }
  }
  return i0;
}

static real amplgsl_sf_bessel_I1_scaled(arglist *al) {
  real x = al->ra[0];
  real i1 = gsl_sf_bessel_I1_scaled(x);
  if (al->derivs) {
    real x_div_absx = x / abs(x);
    real i0 = gsl_sf_bessel_I0_scaled(x), i2 = gsl_sf_bessel_In_scaled(2, x);
    *al->derivs = 0.5 * i0 - x_div_absx * i1 + 0.5 * i2;
    if (al->hes) {
      *al->hes = -x_div_absx * i0 + 1.75 * i1 - x_div_absx * i2 +
          0.25 * gsl_sf_bessel_In_scaled(3, x);
    }
  }
  return i1;
}

static real amplgsl_sf_bessel_I1(arglist *al) {
  real x = al->ra[0];
  real i1 = gsl_sf_bessel_I1(x);
  if (al->derivs) {
    *al->derivs = 0.5 * (gsl_sf_bessel_I0(x) + gsl_sf_bessel_In(2, x));
    if (al->hes)
      *al->hes = 0.25 * (gsl_sf_bessel_In(3, x) + 3 * i1);
  }
  return i1;
}

static real amplgsl_sf_bessel_K0(arglist *al) {
  real x = al->ra[0];
  real k0 = gsl_sf_bessel_K0(x);
  if (al->derivs) {
    *al->derivs = -gsl_sf_bessel_K1(x);
    if (al->hes)
      *al->hes = 0.5 * (gsl_sf_bessel_Kn(2, x) + k0);
  }
  return k0;
}

static real amplgsl_sf_bessel_K0_scaled(arglist *al) {
  real x = al->ra[0];
  real k0 = gsl_sf_bessel_K0_scaled(x);
  if (al->derivs) {
    real k1 = gsl_sf_bessel_K1_scaled(x);
    *al->derivs = k0 - k1;
    if (al->hes)
      *al->hes = 1.5 * k0 - 2 * k1 + 0.5 * gsl_sf_bessel_Kn_scaled(2, x);
  }
  return k0;
}

static real amplgsl_sf_bessel_K1(arglist *al) {
  real x = al->ra[0];
  real k1 = gsl_sf_bessel_K1(x);
  if (al->derivs) {
    *al->derivs = -0.5 * (gsl_sf_bessel_K0(x) + gsl_sf_bessel_Kn(2, x));
    if (al->hes)
      *al->hes = 0.25 * (gsl_sf_bessel_Kn(3, x) + 3 * k1);
  }
  return k1;
}

static real amplgsl_sf_bessel_K1_scaled(arglist *al) {
  real x = al->ra[0];
  real k1 = gsl_sf_bessel_K1_scaled(x);
  if (al->derivs) {
    real k0 = gsl_sf_bessel_K0_scaled(x), k2 = gsl_sf_bessel_Kn_scaled(2, x);
    *al->derivs = -0.5 * k0 + k1 - 0.5 * k2;
    if (al->hes)
      *al->hes = -k0 + 1.75 * k1 - k2 + 0.25 * gsl_sf_bessel_Kn_scaled(3, x);
  }
  return k1;
}

static real amplgsl_sf_bessel_j0(arglist *al) {
  real x = al->ra[0];
  real j0 = gsl_sf_bessel_j0(x);
  if (al->derivs) {
    real x_squared = x * x;
    *al->derivs = (x * cos(x) - sin(x)) / x_squared;
    if (al->hes)
      *al->hes = ((2 - x_squared) * sin(x) - 2 * x * cos(x)) / (x_squared * x);
  }
  return j0;
}

void funcadd_ASL(AmplExports *ae) {
  // Don't call abort on error.
  gsl_set_error_handler_off();

  // Elementary Functions
  addfunc("gsl_log1p", amplgsl_log1p, 0, 1, 0);
  addfunc("gsl_expm1", amplgsl_expm1, 0, 1, 0);
  addfunc("gsl_hypot", amplgsl_hypot, 0, 2, 0);
  addfunc("gsl_hypot3", amplgsl_hypot3, 0, 3, 0);
  // AMPL has built-in functions acosh, asinh and atanh so wrappers
  // are not provided for their GSL equivalents.

  // Airy Functions
  addfunc("gsl_sf_airy_Ai", amplgsl_sf_airy_Ai, 0, 1, 0);
  addfunc("gsl_sf_airy_Bi", amplgsl_sf_airy_Bi, 0, 1, 0);
  addfunc("gsl_sf_airy_Ai_scaled", amplgsl_sf_airy_Ai_scaled, 0, 1, 0);
  addfunc("gsl_sf_airy_Bi_scaled", amplgsl_sf_airy_Bi_scaled, 0, 1, 0);

  // Bessel Functions
  addfunc("gsl_sf_bessel_J0", amplgsl_sf_bessel_J0, 0, 1, 0);
  addfunc("gsl_sf_bessel_J1", amplgsl_sf_bessel_J1, 0, 1, 0);
  // TODO: gsl_sf_bessel_Jn

  // Irregular Cylindrical Bessel Functions
  addfunc("gsl_sf_bessel_Y0", amplgsl_sf_bessel_Y0, 0, 1, 0);
  addfunc("gsl_sf_bessel_Y1", amplgsl_sf_bessel_Y1, 0, 1, 0);
  // TODO: gsl_sf_bessel_Yn

  // Regular Modified Cylindrical Bessel Functions
  addfunc("gsl_sf_bessel_I0", amplgsl_sf_bessel_I0, 0, 1, 0);
  addfunc("gsl_sf_bessel_I1", amplgsl_sf_bessel_I1, 0, 1, 0);
  addfunc("gsl_sf_bessel_I0_scaled", amplgsl_sf_bessel_I0_scaled, 0, 1, 0);
  addfunc("gsl_sf_bessel_I1_scaled", amplgsl_sf_bessel_I1_scaled, 0, 1, 0);
  // TODO: gsl_sf_bessel_In

  // Irregular Modified Cylindrical Bessel Functions
  addfunc("gsl_sf_bessel_K0", amplgsl_sf_bessel_K0, 0, 1, 0);
  addfunc("gsl_sf_bessel_K1", amplgsl_sf_bessel_K1, 0, 1, 0);
  addfunc("gsl_sf_bessel_K0_scaled", amplgsl_sf_bessel_K0_scaled, 0, 1, 0);
  addfunc("gsl_sf_bessel_K1_scaled", amplgsl_sf_bessel_K1_scaled, 0, 1, 0);
  // TODO: gsl_sf_bessel_Kn

  // Regular Spherical Bessel Functions
  addfunc("gsl_sf_bessel_j0", amplgsl_sf_bessel_j0, 0, 1, 0);
  // TODO: j1, j2, jl
  // TODO: y0, y1, y2, yl
  // TODO: i0_scaled, i1_scaled, i2_scaled, il_scaled
  // TODO: k0_scaled, k1_scaled, k2_scaled, kl_scaled
  // TODO: Jnu, Ynu, Inu, Inu_scaled, Knu, Knu_scaled
}
