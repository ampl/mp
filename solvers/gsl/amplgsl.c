// AMPL bindings for the GNU Scientific Library.

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include "solvers/funcadd.h"

static real amplgsl_log1p(arglist *al) {
  real x = al->ra[0];
  real result = gsl_log1p(x);
  if (al->derivs) {
    real deriv = *al->derivs = 1 / (x + 1);
    if (al->hes)
      *al->hes = -deriv * deriv;
  }
  return result;
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
  addfunc("gsl_log1p", amplgsl_log1p, 0, 1, 0);
  addfunc("gsl_sf_bessel_J0", amplgsl_sf_bessel_J0, 0, 1, 0);
}
