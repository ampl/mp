/* AMPL bindings for the GNU Scientific Library. */

#include <math.h>
#include <stdarg.h>

#include <gsl/gsl_errno.h>
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

#include "solvers/funcadd.h"

enum { MAX_ERROR_MESSAGE_SIZE = 100 };

/* Formats the error message and stores it in al->Errmsg. */
static void error(arglist *al, const char *format, ...) {
  va_list args;
  al->Errmsg = al->AE->Tempmem(al->TMI, MAX_ERROR_MESSAGE_SIZE);
  va_start(args, format);
  al->AE->VsnprintF(al->Errmsg, MAX_ERROR_MESSAGE_SIZE, format, args);
  va_end(args);
}

/*
 * Checks the arguments of a zero function such as gsl_sf_airy_Ai_scaled:
 * - argument with the specified index should be representable as unsigned int
 * - al->derivs should be null
 */
static int check_zero_func_args(arglist *al, unsigned s_index) {
  real arg = al->ra[s_index];
  if ((unsigned)arg != arg) {
    error(al, "argument 's' can't be represented as unsigned int, s = %g", arg);
    return 0;
  }
  if (al->derivs) {
    /* Derivative information is requested, so the argument is not constant. */
    error(al, "argument 's' is not constant");
    return 0;
  }
  return 1;
}

/* Checks if the argument is within the bounds for derivative computation. */
static int check_deriv_arg(arglist *al, int arg, int min, int max) {
  if (arg < min) {
    error(al, "can't compute derivative: argument 'n' too small, n = %d", arg);
    return 0;
  }
  if (arg > max) {
    error(al, "can't compute derivative: argument 'n' too large, n = %d", arg);
    return 0;
  }
  return 1;
}

/* Flags for check_bessel_args */
enum {
  DERIV_INT_MIN = 1 /* Derivative can be computed for n = INT_MIN */
};

/*
 * Checks whether the first argument is constant and reports error if not.
 * Returns 1 iff the first argument is constant.
 */
static int check_const_arg(arglist *al, const char *name) {
  if (al->dig && al->dig[0])
    return 1;
  /* Derivative information is requested, so the argument is not constant. */
  error(al, "argument '%s' is not constant", name);
  return 0;
}

/* Checks if the argument with the specified index is representable as int. */
static int check_int_arg(arglist *al, unsigned index, const char *name) {
  real arg = al->ra[index];
  if ((int)arg != arg) {
    error(al, "argument '%s' can't be represented as int, %s = %g",
        name, name, arg);
    return 0;
  }
  return 1;
}

/* Checks the arguments of a Bessel function. */
static int check_bessel_args(arglist *al, int flags) {
  int n = al->ra[0];
  if (!check_int_arg(al, 0, "n"))
    return 0;
  if (al->derivs) {
    int deriv_min = INT_MIN + ((flags & DERIV_INT_MIN) != 0 ? 0 : 1);
    if (!al->dig || !al->dig[0]) {
      /* Can't compute derivative with respect to an integer argument. */
      error(al, "argument 'n' is not constant");
      return 0;
    }
    if ((al->hes && !check_deriv_arg(al, n, INT_MIN + 2, INT_MAX - 2)) ||
        !check_deriv_arg(al, n, deriv_min, INT_MAX - 1))
      return 0;
  }
  return 1;
}

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

static real amplgsl_sf_airy_zero_Ai(arglist *al) {
  return check_zero_func_args(al, 0) ? gsl_sf_airy_zero_Ai(al->ra[0]) : 0;
}

static real amplgsl_sf_airy_zero_Bi(arglist *al) {
  return check_zero_func_args(al, 0) ? gsl_sf_airy_zero_Bi(al->ra[0]) : 0;
}

static real amplgsl_sf_airy_zero_Ai_deriv(arglist *al) {
  return check_zero_func_args(al, 0) ? gsl_sf_airy_zero_Ai_deriv(al->ra[0]) : 0;
}

static real amplgsl_sf_airy_zero_Bi_deriv(arglist *al) {
  return check_zero_func_args(al, 0) ? gsl_sf_airy_zero_Bi_deriv(al->ra[0]) : 0;
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

static real amplgsl_sf_bessel_Jn(arglist *al) {
  int n = al->ra[0];
  real x = al->ra[1];
  real jn = 0;
  if (!check_bessel_args(al, 0))
    return 0;
  jn = gsl_sf_bessel_Jn(n, x);
  if (al->derivs) {
    al->derivs[1] = 0.5 *
        (gsl_sf_bessel_Jn(n - 1, x) - gsl_sf_bessel_Jn(n + 1, x));
    if (al->hes) {
      al->hes[2] = 0.25 *
          (gsl_sf_bessel_Jn(n - 2, x) - 2 * jn + gsl_sf_bessel_Jn(n + 2, x));
    }
  }
  return jn;
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

static real amplgsl_sf_bessel_Yn(arglist *al) {
  int n = al->ra[0];
  real x = al->ra[1];
  real yn = 0;
  if (!check_bessel_args(al, 0))
    return 0;
  yn = gsl_sf_bessel_Yn(n, x);
  if (al->derivs) {
    al->derivs[1] = 0.5 *
        (gsl_sf_bessel_Yn(n - 1, x) - gsl_sf_bessel_Yn(n + 1, x));
    if (al->hes) {
      al->hes[2] = 0.25 *
          (gsl_sf_bessel_Yn(n - 2, x) - 2 * yn + gsl_sf_bessel_Yn(n + 2, x));
    }
  }
  return yn;
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

static real amplgsl_sf_bessel_In(arglist *al) {
  int n = al->ra[0];
  real x = al->ra[1];
  real in = 0;
  if (!check_bessel_args(al, 0))
    return 0;
  in = gsl_sf_bessel_In(n, x);
  if (al->derivs) {
    al->derivs[1] = 0.5 *
        (gsl_sf_bessel_In(n - 1, x) + gsl_sf_bessel_In(n + 1, x));
    if (al->hes) {
      al->hes[2] = 0.25 *
          (gsl_sf_bessel_In(n - 2, x) + 2 * in + gsl_sf_bessel_In(n + 2, x));
    }
  }
  return in;
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

static real amplgsl_sf_bessel_In_scaled(arglist *al) {
  int n = al->ra[0];
  real x = al->ra[1];
  real in = 0;
  if (!check_bessel_args(al, 0))
    return 0;
  in = gsl_sf_bessel_In_scaled(n, x);
  if (al->derivs) {
    real in_minus1 = gsl_sf_bessel_In_scaled(n - 1, x);
    real in_plus1 = gsl_sf_bessel_In_scaled(n + 1, x);
    al->derivs[1] = 0.5 * (in_minus1 - (2 * x * in) / abs(x) + in_plus1);
    if (al->hes) {
      al->hes[2] =
          (abs(x) * (gsl_sf_bessel_In_scaled(n - 2, x) + 6 * in +
                     gsl_sf_bessel_In_scaled(n + 2, x)) -
           4 * x * (in_minus1 + in_plus1)) / (4 * abs(x));
    }
  }
  return in;
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

static real amplgsl_sf_bessel_Kn(arglist *al) {
  int n = al->ra[0];
  real x = al->ra[1];
  real kn = 0;
  if (!check_bessel_args(al, 0))
    return 0;
  kn = gsl_sf_bessel_Kn(n, x);
  if (al->derivs) {
    al->derivs[1] = -0.5 *
        (gsl_sf_bessel_Kn(n - 1, x) + gsl_sf_bessel_Kn(n + 1, x));
    if (al->hes) {
      al->hes[2] = 0.25 *
          (gsl_sf_bessel_Kn(n - 2, x) + 2 * kn + gsl_sf_bessel_Kn(n + 2, x));
    }
  }
  return kn;
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

static real amplgsl_sf_bessel_Kn_scaled(arglist *al) {
  int n = al->ra[0];
  real x = al->ra[1];
  real kn = 0;
  if (!check_bessel_args(al, 0))
    return 0;
  kn = gsl_sf_bessel_Kn_scaled(n, x);
  if (al->derivs) {
    real kn_minus1 = gsl_sf_bessel_Kn_scaled(n - 1, x);
    real kn_plus1 = gsl_sf_bessel_Kn_scaled(n + 1, x);
    al->derivs[1] = -0.5 * (kn_minus1 - 2 * kn + kn_plus1);
    if (al->hes) {
      al->hes[2] = 0.25 *
          (gsl_sf_bessel_Kn_scaled(n - 2, x) - 4 * kn_minus1 + 6 * kn -
              4 * kn_plus1 + gsl_sf_bessel_Kn_scaled(n + 2, x));
    }
  }
  return kn;
}

static real amplgsl_sf_bessel_j0(arglist *al) {
  real x = al->ra[0];
  if (al->derivs) {
    real x_squared = x * x;
    *al->derivs = (x * cos(x) - sin(x)) / x_squared;
    if (al->hes)
      *al->hes = ((2 - x_squared) * sin(x) - 2 * x * cos(x)) / (x_squared * x);
  }
  return gsl_sf_bessel_j0(x);
}

static real amplgsl_sf_bessel_j1(arglist *al) {
  real x = al->ra[0];
  real j1 = gsl_sf_bessel_j1(x);
  if (al->derivs) {
    real sinx = sin(x);
    *al->derivs = (sinx - 2 * j1) / x;
    if (al->hes) {
      real x_squared = x * x;
      *al->hes = (x * (x_squared - 6) * cos(x) - 3 * (x_squared - 2) * sinx) /
        (x_squared * x_squared);
    }
  }
  return j1;
}

static real amplgsl_sf_bessel_j2(arglist *al) {
  real x = al->ra[0];
  real j2 = gsl_sf_bessel_j2(x);
  if (al->derivs) {
    *al->derivs = gsl_sf_bessel_j1(x) - 3 * j2 / x;
    if (al->hes) {
      real x_pow2 = x * x, x_pow4 = x_pow2 * x_pow2;
      *al->hes = (x * (5 * x_pow2 - 36) * cos(x) +
          (x_pow4 - 17 * x_pow2 + 36) * sin(x)) / (x_pow4 * x);
    }
  }
  return j2;
}

static real amplgsl_sf_bessel_jl(arglist *al) {
  int n = al->ra[0];
  real x = al->ra[1];
  real jn = 0;
  if (!check_bessel_args(al, DERIV_INT_MIN))
    return 0;
  jn = gsl_sf_bessel_jl(n, x);
  if (al->derivs) {
    real jn_plus1 = gsl_sf_bessel_jl(n + 1, x);
    al->derivs[1] = n * jn / x - jn_plus1;
    if (al->hes) {
      real x_squared = x * x;
      al->hes[2] = (
          x_squared * gsl_sf_bessel_jl(n - 2, x) -
          2 * x * gsl_sf_bessel_jl(n - 1, x) -
          (2 * x_squared - 3) * jn + 2 * x * jn_plus1 +
          x_squared * gsl_sf_bessel_jl(n + 2, x)) / (4 * x_squared);
    }
  }
  return jn;
}

static real amplgsl_sf_bessel_y0(arglist *al) {
  real x = al->ra[0];
  if (al->derivs) {
    real x_squared = x * x;
    *al->derivs = (x * sin(x) + cos(x)) / x_squared;
    if (al->hes)
      *al->hes = ((x_squared - 2) * cos(x) - 2 * x * sin(x)) / (x_squared * x);
  }
  return gsl_sf_bessel_y0(x);
}

static real amplgsl_sf_bessel_y1(arglist *al) {
  real x = al->ra[0];
  real y1 = gsl_sf_bessel_y1(x);
  if (al->derivs) {
    *al->derivs = -(2 * y1 + cos(x)) / x;
    if (al->hes) {
      real x_squared = x * x;
      *al->hes = (x * (x_squared - 6) * sin(x) +
          3 * (x_squared - 2) * cos(x)) / (x_squared * x_squared);
    }
  }
  return y1;
}

static real amplgsl_sf_bessel_y2(arglist *al) {
  real x = al->ra[0];
  real y2 = gsl_sf_bessel_y2(x);
  if (al->derivs) {
    real y1 = gsl_sf_bessel_y1(x);
    *al->derivs = y1 - (3 * y2) / x;
    if (al->hes) {
      real x_squared = x * x;
      *al->hes = ((36 - 5 * x_squared) * y1 -
          (x_squared - 12) * cos(x)) / (x_squared * x);
    }
  }
  return y2;
}

static real amplgsl_sf_bessel_yl(arglist *al) {
  int n = al->ra[0];
  real x = al->ra[1];
  real yn = 0;
  if (!check_bessel_args(al, 0))
    return 0;
  yn = gsl_sf_bessel_yl(n, x);
  if (al->derivs) {
    real yn_minus1 = gsl_sf_bessel_yl(n - 1, x);
    real yn_plus1 = gsl_sf_bessel_yl(n + 1, x);
    al->derivs[1] = 0.5 * (yn_minus1 - yn / x - yn_plus1);
    if (al->hes) {
      real x_squared = x * x;
      al->hes[2] = (
          x_squared * gsl_sf_bessel_yl(n - 2, x) - 2 * x * yn_minus1 -
          (2 * x_squared - 3) * yn + 2 * x * yn_plus1 +
          x_squared * gsl_sf_bessel_yl(n + 2, x)) / (4 * x_squared);
    }
  }
  return yn;
}

static real amplgsl_sf_bessel_i0_scaled(arglist *al) {
  real x = al->ra[0];
  real i0 = gsl_sf_bessel_i0_scaled(x);
  if (al->derivs) {
    real hyp_coef = exp(-abs(x)) * sqrt(1 / x) / sqrt(x);
    real i_minus1 = hyp_coef * cosh(x);
    real i1 = gsl_sf_bessel_i1_scaled(x);
    real coef = -(1 + 2 * abs(x)) / x;
    *al->derivs = 0.5 * (i_minus1 + coef * i0 + i1);
    if (al->hes) {
      coef *= 2;
      *al->hes = 0.25 * (
          hyp_coef * sinh(x) - i_minus1 / x +
          coef * i_minus1 +
          (3 + 6 * x * x + 4 * abs(x)) * i0 / (x * x) +
          coef * i1 +
          gsl_sf_bessel_il_scaled(2, x));
    }
  }
  return i0;
}

static real amplgsl_sf_bessel_i1_scaled(arglist *al) {
  real x = al->ra[0];
  real i1 = gsl_sf_bessel_i1_scaled(x);
  if (al->derivs) {
    real i0 = gsl_sf_bessel_i0_scaled(x);
    real i2 = gsl_sf_bessel_i2_scaled(x);
    real coef = -(1 + 2 * abs(x)) / x;
    *al->derivs = 0.5 * (i0 + coef * i1 + i2);
    if (al->hes) {
      coef *= 2;
      *al->hes = 0.25 * (
          exp(-abs(x)) * sqrt(1 / x) * cosh(x) / sqrt(x) +
          coef * i0 +
          (3 + 6 * x * x + 4 * abs(x)) * i1 / (x * x) +
          coef * i2 +
          gsl_sf_bessel_il_scaled(3, x));
    }
  }
  return i1;
}

static real amplgsl_sf_bessel_i2_scaled(arglist *al) {
  real x = al->ra[0];
  real i2 = gsl_sf_bessel_i2_scaled(x);
  if (al->derivs) {
    real i1 = gsl_sf_bessel_i1_scaled(x);
    real i3 = gsl_sf_bessel_il_scaled(3, x);
    real coef = -(1 + 2 * abs(x)) / x;
    *al->derivs = 0.5 * (i1 + coef * i2 + i3);
    if (al->hes) {
      coef *= 2;
      *al->hes = 0.25 * (
          gsl_sf_bessel_i0_scaled(x) +
          coef * i1 +
          (3 + 6 * x * x + 4 * abs(x)) * i2 / (x * x) +
          coef * i3 +
          gsl_sf_bessel_il_scaled(4, x));
    }
  }
  return i2;
}

static real amplgsl_sf_bessel_il_scaled(arglist *al) {
  int n = al->ra[0];
  real x = al->ra[1];
  real in = 0;
  if (!check_bessel_args(al, 0))
    return 0;
  in = gsl_sf_bessel_il_scaled(n, x);
  if (al->derivs) {
    real in_minus1 = gsl_sf_bessel_il_scaled(n - 1, x);
    real in_plus1 = gsl_sf_bessel_il_scaled(n + 1, x);
    real coef = -(1 + 2 * abs(x)) / x;
    al->derivs[1] = 0.5 * (in_minus1 + coef * in + in_plus1);
    if (al->hes) {
      coef *= 2;
      al->hes[2] = 0.25 * (
          gsl_sf_bessel_il_scaled(n - 2, x) +
          coef * in_minus1 +
          (3 + 4 * abs(x) + 6 * x * x) * in / (x * x) +
          coef * in_plus1 +
          gsl_sf_bessel_il_scaled(n + 2, x));
    }
  }
  return in;
}

static real amplgsl_sf_bessel_k0_scaled(arglist *al) {
  real x = al->ra[0];
  if (al->derivs) {
    real pi_sqrt_inv_x = M_PI * sqrt(1 / x);
    *al->derivs = -pi_sqrt_inv_x / (2 * pow(x, 1.5));
    if (al->hes)
      *al->hes = pi_sqrt_inv_x / pow(x, 2.5);
  }
  return gsl_sf_bessel_k0_scaled(x);
}

static real amplgsl_sf_bessel_k1_scaled(arglist *al) {
  real x = al->ra[0];
  if (al->derivs) {
    real pi_sqrt_inv_x = M_PI * sqrt(1 / x);
    *al->derivs = -(pi_sqrt_inv_x * (x + 2)) / (2 * pow(x, 2.5));
    if (al->hes)
      *al->hes = (pi_sqrt_inv_x * (x + 3)) / pow(x, 3.5);
  }
  return gsl_sf_bessel_k1_scaled(x);
}

static real amplgsl_sf_bessel_k2_scaled(arglist *al) {
  real x = al->ra[0];
  if (al->derivs) {
    real pi_sqrt_inv_x = M_PI * sqrt(1 / x);
    *al->derivs = -pi_sqrt_inv_x * (x + 3) * (x + 3) / (2 * pow(x, 3.5));
    if (al->hes)
      *al->hes = pi_sqrt_inv_x * (x * x + 9 * x + 18) / pow(x, 4.5);
  }
  return gsl_sf_bessel_k2_scaled(x);
}

static real amplgsl_sf_bessel_kl_scaled(arglist *al) {
  int n = al->ra[0];
  real x = al->ra[1];
  real kn = 0;
  if (!check_bessel_args(al, 0))
    return 0;
  kn = gsl_sf_bessel_kl_scaled(n, x);
  if (al->derivs) {
    real kn_minus1 = gsl_sf_bessel_kl_scaled(n - 1, x);
    real kn_plus1 = gsl_sf_bessel_kl_scaled(n + 1, x);
    real coef = (1 - 2 * x) / x;
    al->derivs[1] = -0.5 * (kn_minus1 + coef * kn + kn_plus1);
    if (al->hes) {
      coef *= 2;
      al->hes[2] = 0.25 * (
          gsl_sf_bessel_kl_scaled(n - 2, x) +
          coef * kn_minus1 +
          (3 - 4 * x + 6 * x * x) * kn / (x * x) +
          coef * kn_plus1 +
          gsl_sf_bessel_kl_scaled(n + 2, x));
    }
  }
  return kn;
}

static real amplgsl_sf_bessel_Jnu(arglist *al) {
  real n = al->ra[0];
  real x = al->ra[1];
  real jn = gsl_sf_bessel_Jnu(n, x);
  if (al->derivs) {
    if (!check_const_arg(al, "nu"))
      return 0;
    al->derivs[1] = 0.5 *
        (gsl_sf_bessel_Jnu(n - 1, x) - gsl_sf_bessel_Jnu(n + 1, x));
    if (al->hes) {
      al->hes[2] = 0.25 *
          (gsl_sf_bessel_Jnu(n - 2, x) - 2 * jn + gsl_sf_bessel_Jnu(n + 2, x));
    }
  }
  return jn;
}

static real amplgsl_sf_bessel_Ynu(arglist *al) {
  real n = al->ra[0];
  real x = al->ra[1];
  real yn = gsl_sf_bessel_Ynu(n, x);
  if (al->derivs) {
    if (!check_const_arg(al, "nu"))
      return 0;
    al->derivs[1] = 0.5 *
        (gsl_sf_bessel_Ynu(n - 1, x) - gsl_sf_bessel_Ynu(n + 1, x));
    if (al->hes) {
      al->hes[2] = 0.25 *
          (gsl_sf_bessel_Ynu(n - 2, x) - 2 * yn + gsl_sf_bessel_Ynu(n + 2, x));
    }
  }
  return yn;
}

static real amplgsl_sf_bessel_Inu(arglist *al) {
  real n = al->ra[0];
  real x = al->ra[1];
  real in = gsl_sf_bessel_Inu(n, x);
  if (al->derivs) {
    if (!check_const_arg(al, "nu"))
      return 0;
    al->derivs[1] = 0.5 *
        (gsl_sf_bessel_Inu(n - 1, x) + gsl_sf_bessel_Inu(n + 1, x));
    if (al->hes) {
      al->hes[2] = 0.25 *
          (gsl_sf_bessel_Inu(n - 2, x) + 2 * in + gsl_sf_bessel_Inu(n + 2, x));
    }
  }
  return in;
}

static real amplgsl_sf_bessel_Inu_scaled(arglist *al) {
  real n = al->ra[0];
  real x = al->ra[1];
  real in = gsl_sf_bessel_Inu_scaled(n, x);
  if (al->derivs) {
    real in_minus1 = 0, in_plus1 = 0;
    if (!check_const_arg(al, "nu"))
      return 0;
    in_minus1 = gsl_sf_bessel_Inu_scaled(n - 1, x);
    in_plus1 = gsl_sf_bessel_Inu_scaled(n + 1, x);
    al->derivs[1] = 0.5 * (in_minus1 - (2 * x * in) / abs(x) + in_plus1);
    if (al->hes) {
      al->hes[2] =
          (abs(x) * (gsl_sf_bessel_Inu_scaled(n - 2, x) + 6 * in +
                     gsl_sf_bessel_Inu_scaled(n + 2, x)) -
           4 * x * (in_minus1 + in_plus1)) / (4 * abs(x));
    }
  }
  return in;
}

static real amplgsl_sf_bessel_Knu(arglist *al) {
  real n = al->ra[0];
  real x = al->ra[1];
  real kn = gsl_sf_bessel_Knu(n, x);
  if (al->derivs) {
    if (!check_const_arg(al, "nu"))
      return 0;
    al->derivs[1] = -0.5 *
        (gsl_sf_bessel_Knu(n - 1, x) + gsl_sf_bessel_Knu(n + 1, x));
    if (al->hes) {
      al->hes[2] = 0.25 *
          (gsl_sf_bessel_Knu(n - 2, x) + 2 * kn + gsl_sf_bessel_Knu(n + 2, x));
    }
  }
  return kn;
}

static real amplgsl_sf_bessel_lnKnu(arglist *al) {
  real n = al->ra[0];
  real x = al->ra[1];
  if (al->derivs) {
    real kn = 0, kn_minus1_plus1 = 0;
    if (!check_const_arg(al, "nu"))
      return 0;
    kn = gsl_sf_bessel_Knu(n, x);
    kn_minus1_plus1 =
        gsl_sf_bessel_Knu(n - 1, x) + gsl_sf_bessel_Knu(n + 1, x);
    al->derivs[1] = -0.5 * kn_minus1_plus1 / kn;
    if (al->hes) {
      al->hes[2] = 0.25 *
          (kn * (gsl_sf_bessel_Knu(n - 2, x) + 2 * kn +
          gsl_sf_bessel_Knu(n + 2, x)) - kn_minus1_plus1 * kn_minus1_plus1) /
          (kn * kn);
    }
  }
  return gsl_sf_bessel_lnKnu(n, x);
}

static real amplgsl_sf_bessel_Knu_scaled(arglist *al) {
  real n = al->ra[0];
  real x = al->ra[1];
  real kn = gsl_sf_bessel_Knu_scaled(n, x);
  if (al->derivs) {
    real kn_minus1 = 0, kn_plus1 = 0;
    if (!check_const_arg(al, "nu"))
      return 0;
    kn_minus1 = gsl_sf_bessel_Knu_scaled(n - 1, x);
    kn_plus1 = gsl_sf_bessel_Knu_scaled(n + 1, x);
    al->derivs[1] = -0.5 * (kn_minus1 - 2 * kn + kn_plus1);
    if (al->hes) {
      al->hes[2] = 0.25 *
          (gsl_sf_bessel_Knu_scaled(n - 2, x) - 4 * kn_minus1 + 6 * kn -
              4 * kn_plus1 + gsl_sf_bessel_Knu_scaled(n + 2, x));
    }
  }
  return kn;
}

static real amplgsl_sf_bessel_zero_J0(arglist *al) {
  return check_zero_func_args(al, 0) ? gsl_sf_bessel_zero_J0(al->ra[0]) : 0;
}

static real amplgsl_sf_bessel_zero_J1(arglist *al) {
  return check_zero_func_args(al, 0) ? gsl_sf_bessel_zero_J1(al->ra[0]) : 0;
}

static real amplgsl_sf_bessel_zero_Jnu(arglist *al) {
  return check_zero_func_args(al, 1) ?
      gsl_sf_bessel_zero_Jnu(al->ra[0], al->ra[1]) : 0;
}

static real amplgsl_sf_clausen(arglist *al) {
  real x = al->ra[0];
  if (al->derivs) {
    *al->derivs = -log(2 * sin(0.5 * x));
    if (al->hes)
      *al->hes = -0.5 * tan(0.5 * M_PI - x);
  }
  return gsl_sf_clausen(x);
}

static real amplgsl_sf_hydrogenicR_1(arglist *al) {
  real Z = al->ra[0], r = al->ra[1];
  if (al->derivs) {
    real *derivs = al->derivs;
    real exp_minusZr = exp(-Z * r);
    derivs[0] = sqrt(Z) * exp_minusZr * (3 - 2 * r * Z);
    derivs[1] = -2 * pow(Z, 2.5) * exp_minusZr;
    if (al->hes) {
      real *hes = al->hes;
      hes[0] = (exp_minusZr * (4 * r * r * Z * Z - 12 * r * Z + 3)) /
          (2 *sqrt(Z));
      hes[1] = pow(Z, 1.5) * exp_minusZr * (2 * r * Z - 5);
      hes[2] = 2 * pow(Z, 3.5) * exp_minusZr;
    }
  }
  return gsl_sf_hydrogenicR_1(Z, r);
}

static real amplgsl_sf_hydrogenicR(arglist *al) {
  if (!check_int_arg(al, 0, "n") || !check_int_arg(al, 1, "l"))
    return 0;
  if (al->derivs) {
    error(al, "derivative is not provided");
    return 0;
  }
  return gsl_sf_hydrogenicR(al->ra[0], al->ra[1], al->ra[2], al->ra[3]);
}

static real amplgsl_sf_coulomb_CL(arglist *al) {
  gsl_sf_result result = {0};
  if (al->derivs) {
    error(al, "derivative is not provided");
    return 0;
  }
  return gsl_sf_coulomb_CL_e(al->ra[0], al->ra[1], &result) ?
      GSL_NAN : result.val;
}

static int check_coupling_args(arglist *al, const char *const* arg_names) {
  unsigned i = 0, n_args = al->n;
  for (; i < n_args; ++i) {
    if (!check_int_arg(al, i, arg_names[i]))
      return 0;
  }
  if (al->derivs) {
    error(al, "arguments are not constant");
    return 0;
  }
  return 1;
}

static real amplgsl_sf_coupling_3j(arglist *al) {
  static const char *ARG_NAMES[] = {
      "two_ja", "two_jb", "two_jc",
      "two_ma", "two_mb", "two_mc"
  };
  if (!check_coupling_args(al, ARG_NAMES))
    return 0;
  return gsl_sf_coupling_3j(
      al->ra[0], al->ra[1], al->ra[2], al->ra[3], al->ra[4], al->ra[5]);
}

static real amplgsl_sf_coupling_6j(arglist *al) {
  static const char *ARG_NAMES[] = {
      "two_ja", "two_jb", "two_jc",
      "two_jd", "two_je", "two_jf"
  };
  if (!check_coupling_args(al, ARG_NAMES))
    return 0;
  return gsl_sf_coupling_6j(
      al->ra[0], al->ra[1], al->ra[2], al->ra[3], al->ra[4], al->ra[5]);
}

static real amplgsl_sf_coupling_9j(arglist *al) {
  static const char *ARG_NAMES[] = {
      "two_ja", "two_jb", "two_jc",
      "two_jd", "two_je", "two_jf",
      "two_jg", "two_jh", "two_ji"
  };
  if (!check_coupling_args(al, ARG_NAMES))
    return 0;
  return gsl_sf_coupling_9j(
      al->ra[0], al->ra[1], al->ra[2], al->ra[3], al->ra[4], al->ra[5],
      al->ra[6], al->ra[7], al->ra[8]);
}

static real amplgsl_sf_dawson(arglist *al) {
  real x = al->ra[0];
  real f = gsl_sf_dawson(x);
  if (al->derivs) {
    real deriv = *al->derivs = 1 - 2 * x * f;
    if (al->hes)
      *al->hes = - 2 * (f + x * deriv);
  }
  return f;
}

static real debye(arglist *al, int n, double (*func)(double)) {
  real x = al->ra[0];
  real f = func(x);
  if (al->derivs) {
    real exp_x = exp(x);
    real deriv = *al->derivs = n * (1 / (exp_x - 1) - f / x);
    if (al->hes) {
      *al->hes = n * (-exp_x / ((exp_x - 1) * (exp_x - 1)) +
          f / (x * x) - deriv / x);
    }
  }
  return f;
}

static real amplgsl_sf_debye_1(arglist *al) {
  return debye(al, 1, gsl_sf_debye_1);
}

static real amplgsl_sf_debye_2(arglist *al) {
  return debye(al, 2, gsl_sf_debye_2);
}

static real amplgsl_sf_debye_3(arglist *al) {
  return debye(al, 3, gsl_sf_debye_3);
}

static real amplgsl_sf_debye_4(arglist *al) {
  return debye(al, 4, gsl_sf_debye_4);
}

static real amplgsl_sf_debye_5(arglist *al) {
  return debye(al, 5, gsl_sf_debye_5);
}

static real amplgsl_sf_debye_6(arglist *al) {
  return debye(al, 6, gsl_sf_debye_6);
}

static real amplgsl_sf_dilog(arglist *al) {
  real x = al->ra[0];
  if (al->derivs) {
    real deriv = *al->derivs = x != 0 ?
        -GSL_REAL(gsl_complex_log(gsl_complex_rect(1 - x, 0))) / x : 1;
    if (al->hes)
      *al->hes = x != 0 ? (1 / (1 - x) - deriv) / x : 0.5;
  }
  return gsl_sf_dilog(x);
}

void funcadd_ASL(AmplExports *ae) {
  /* Don't call abort on error. */
  gsl_set_error_handler_off();

  /* Elementary Functions */
  addfunc("gsl_log1p", amplgsl_log1p, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_expm1", amplgsl_expm1, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_hypot", amplgsl_hypot, FUNCADD_REAL_VALUED, 2, 0);
  addfunc("gsl_hypot3", amplgsl_hypot3, FUNCADD_REAL_VALUED, 3, 0);

  /* AMPL has built-in functions acosh, asinh and atanh so wrappers
     are not provided for their GSL equivalents. */

  /* Wrappers for functions operating on complex numbers are not provided
     since this requires support for structures/tuples as function arguments. */

  /* Airy Functions */
  addfunc("gsl_sf_airy_Ai", amplgsl_sf_airy_Ai, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_airy_Bi", amplgsl_sf_airy_Bi, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_airy_Ai_scaled", amplgsl_sf_airy_Ai_scaled,
      FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_airy_Bi_scaled", amplgsl_sf_airy_Bi_scaled,
      FUNCADD_REAL_VALUED, 1, 0);

  /* Zeros of Airy Functions */
  addfunc("gsl_sf_airy_zero_Ai", amplgsl_sf_airy_zero_Ai,
      FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_airy_zero_Bi", amplgsl_sf_airy_zero_Bi,
      FUNCADD_REAL_VALUED, 1, 0);

  /* Zeros of Derivatives of Airy Functions */
  addfunc("gsl_sf_airy_zero_Ai_deriv", amplgsl_sf_airy_zero_Ai_deriv,
      FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_airy_zero_Bi_deriv", amplgsl_sf_airy_zero_Bi_deriv,
      FUNCADD_REAL_VALUED, 1, 0);

  /* Bessel Functions */
  addfunc("gsl_sf_bessel_J0", amplgsl_sf_bessel_J0, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_bessel_J1", amplgsl_sf_bessel_J1, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_bessel_Jn", amplgsl_sf_bessel_Jn, FUNCADD_REAL_VALUED, 2, 0);

  /* Irregular Cylindrical Bessel Functions */
  addfunc("gsl_sf_bessel_Y0", amplgsl_sf_bessel_Y0, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_bessel_Y1", amplgsl_sf_bessel_Y1, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_bessel_Yn", amplgsl_sf_bessel_Yn, FUNCADD_REAL_VALUED, 2, 0);

  /* Regular Modified Cylindrical Bessel Functions */
  addfunc("gsl_sf_bessel_I0", amplgsl_sf_bessel_I0, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_bessel_I1", amplgsl_sf_bessel_I1, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_bessel_In", amplgsl_sf_bessel_In, FUNCADD_REAL_VALUED, 2, 0);
  addfunc("gsl_sf_bessel_I0_scaled", amplgsl_sf_bessel_I0_scaled,
      FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_bessel_I1_scaled", amplgsl_sf_bessel_I1_scaled,
      FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_bessel_In_scaled", amplgsl_sf_bessel_In_scaled,
      FUNCADD_REAL_VALUED, 2, 0);

  /* Irregular Modified Cylindrical Bessel Functions */
  addfunc("gsl_sf_bessel_K0", amplgsl_sf_bessel_K0, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_bessel_K1", amplgsl_sf_bessel_K1, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_bessel_Kn", amplgsl_sf_bessel_Kn, FUNCADD_REAL_VALUED, 2, 0);
  addfunc("gsl_sf_bessel_K0_scaled", amplgsl_sf_bessel_K0_scaled,
      FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_bessel_K1_scaled", amplgsl_sf_bessel_K1_scaled,
      FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_bessel_Kn_scaled", amplgsl_sf_bessel_Kn_scaled,
      FUNCADD_REAL_VALUED, 2, 0);

  /* Regular Spherical Bessel Functions */
  addfunc("gsl_sf_bessel_j0", amplgsl_sf_bessel_j0, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_bessel_j1", amplgsl_sf_bessel_j1, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_bessel_j2", amplgsl_sf_bessel_j2, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_bessel_jl", amplgsl_sf_bessel_jl, FUNCADD_REAL_VALUED, 2, 0);

  /* Irregular Spherical Bessel Functions */
  addfunc("gsl_sf_bessel_y0", amplgsl_sf_bessel_y0, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_bessel_y1", amplgsl_sf_bessel_y1, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_bessel_y2", amplgsl_sf_bessel_y2, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_bessel_yl", amplgsl_sf_bessel_yl, FUNCADD_REAL_VALUED, 2, 0);

  /* Regular Modified Spherical Bessel Functions */
  addfunc("gsl_sf_bessel_i0_scaled", amplgsl_sf_bessel_i0_scaled,
      FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_bessel_i1_scaled", amplgsl_sf_bessel_i1_scaled,
      FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_bessel_i2_scaled", amplgsl_sf_bessel_i2_scaled,
      FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_bessel_il_scaled", amplgsl_sf_bessel_il_scaled,
      FUNCADD_REAL_VALUED, 2, 0);

  /* Irregular Modified Spherical Bessel Functions */
  addfunc("gsl_sf_bessel_k0_scaled", amplgsl_sf_bessel_k0_scaled,
      FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_bessel_k1_scaled", amplgsl_sf_bessel_k1_scaled,
      FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_bessel_k2_scaled", amplgsl_sf_bessel_k2_scaled,
      FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_bessel_kl_scaled", amplgsl_sf_bessel_kl_scaled,
      FUNCADD_REAL_VALUED, 2, 0);

  /* Regular Bessel Function - Fractional Order */
  addfunc("gsl_sf_bessel_Jnu", amplgsl_sf_bessel_Jnu,
      FUNCADD_REAL_VALUED, 2, 0);

  /* Irregular Bessel Functions - Fractional Order */
  addfunc("gsl_sf_bessel_Ynu", amplgsl_sf_bessel_Ynu,
      FUNCADD_REAL_VALUED, 2, 0);

  /* Regular Modified Bessel Functions - Fractional Order */
  addfunc("gsl_sf_bessel_Inu", amplgsl_sf_bessel_Inu,
      FUNCADD_REAL_VALUED, 2, 0);
  addfunc("gsl_sf_bessel_Inu_scaled", amplgsl_sf_bessel_Inu_scaled,
      FUNCADD_REAL_VALUED, 2, 0);

  /* Irregular Modified Bessel Functions - Fractional Order */
  addfunc("gsl_sf_bessel_Knu", amplgsl_sf_bessel_Knu,
      FUNCADD_REAL_VALUED, 2, 0);
  addfunc("gsl_sf_bessel_lnKnu", amplgsl_sf_bessel_lnKnu,
      FUNCADD_REAL_VALUED, 2, 0);
  addfunc("gsl_sf_bessel_Knu_scaled", amplgsl_sf_bessel_Knu_scaled,
      FUNCADD_REAL_VALUED, 2, 0);

  /* Zeros of Regular Bessel Functions */
  addfunc("gsl_sf_bessel_zero_J0", amplgsl_sf_bessel_zero_J0,
      FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_bessel_zero_J1", amplgsl_sf_bessel_zero_J1,
      FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_bessel_zero_Jnu", amplgsl_sf_bessel_zero_Jnu,
      FUNCADD_REAL_VALUED, 2, 0);

  /* Clausen Function */
  addfunc("gsl_sf_clausen", amplgsl_sf_clausen, FUNCADD_REAL_VALUED, 1, 0);

  /* Normalized Hydrogenic Bound States */
  addfunc("gsl_sf_hydrogenicR_1", amplgsl_sf_hydrogenicR_1,
      FUNCADD_REAL_VALUED, 2, 0);
  addfunc("gsl_sf_hydrogenicR", amplgsl_sf_hydrogenicR,
      FUNCADD_REAL_VALUED, 4, 0);

  /* Coulomb Wave Function Normalization Constant */
  addfunc("gsl_sf_coulomb_CL", amplgsl_sf_coulomb_CL,
      FUNCADD_REAL_VALUED, 2, 0);

  /* Coupling Coefficients */
  addfunc("gsl_sf_coupling_3j", amplgsl_sf_coupling_3j,
      FUNCADD_REAL_VALUED, 6, 0);
  addfunc("gsl_sf_coupling_6j", amplgsl_sf_coupling_6j,
      FUNCADD_REAL_VALUED, 6, 0);
  addfunc("gsl_sf_coupling_9j", amplgsl_sf_coupling_9j,
      FUNCADD_REAL_VALUED, 9, 0);

  /* Dawson Function */
  addfunc("gsl_sf_dawson", amplgsl_sf_dawson, FUNCADD_REAL_VALUED, 1, 0);

  /* Debye Functions */
  addfunc("gsl_sf_debye_1", amplgsl_sf_debye_1, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_debye_2", amplgsl_sf_debye_2, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_debye_3", amplgsl_sf_debye_3, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_debye_4", amplgsl_sf_debye_4, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_debye_5", amplgsl_sf_debye_5, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_debye_6", amplgsl_sf_debye_6, FUNCADD_REAL_VALUED, 1, 0);

  /* Dilogarithm */
  addfunc("gsl_sf_dilog", amplgsl_sf_dilog, FUNCADD_REAL_VALUED, 1, 0);

  /* Elliptic Integrals */
  // TODO

  /* Elliptic Functions (Jacobi) */
  // TODO

  /* Error Functions */
  // TODO

  /* Complementary Error Function */
  // TODO

  /* Log Complementary Error Function */
  // TODO
}
