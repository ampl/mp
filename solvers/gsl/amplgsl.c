// AMPL bindings for the GNU Scientific Library.

#include <math.h>

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

void funcadd_ASL(AmplExports *ae) {
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

  addfunc("gsl_sf_bessel_J0", amplgsl_sf_bessel_J0, 0, 1, 0);
}
