/*
 AMPL bindings for GNU Scientific Library.

 Copyright (C) 2012 AMPL Optimization LLC

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization LLC disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Victor Zverovich
 */

#include <math.h>
#include <stdarg.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf.h>

#include "solvers/funcadd.h"

enum { MAX_ERROR_MESSAGE_SIZE = 100 };

static const char *const DERIVS_NOT_PROVIDED = "derivatives are not provided";

/* Computes (x / fabs(x)) * y. Returns 0 if y is 0. */
static double mul_by_sign(double x, double y) {
  return y != 0 ? (x / fabs(x)) * y : 0;
}

/* Formats the error message and stores it in al->Errmsg. */
static void error(arglist *al, const char *format, ...) {
  va_list args;
  al->Errmsg = al->AE->Tempmem(al->TMI, MAX_ERROR_MESSAGE_SIZE);
  va_start(args, format);
  al->AE->VsnprintF(al->Errmsg, MAX_ERROR_MESSAGE_SIZE, format, args);
  va_end(args);
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

/* Reports a function evaluation error printing the specified suffix
 * after the function name. */
static void eval_error_with_suffix(arglist *al, const char *suffix) {
  int n = 0, i = 0;
  al->Errmsg = al->AE->Tempmem(al->TMI, MAX_ERROR_MESSAGE_SIZE);
  n += al->AE->SnprintF(al->Errmsg, MAX_ERROR_MESSAGE_SIZE,
      "can't evaluate %s%s(", al->funcinfo, suffix);
  for (i = 0; i < al->n - 1; ++i) {
    n += al->AE->SnprintF(al->Errmsg + n, MAX_ERROR_MESSAGE_SIZE - n,
        "%g, ", al->ra[i]);
  }
  al->AE->SnprintF(al->Errmsg + n, MAX_ERROR_MESSAGE_SIZE - n,
      "%g)", al->ra[al->n - 1]);
}

/* Reports a function evaluation error. */
static void eval_error(arglist *al) {
  eval_error_with_suffix(al, "");
}

static int check_args(arglist *al) {
  int i = 0;
  for (; i < al->n; ++i) {
    if (gsl_isnan(al->ra[i])) {
      eval_error(al);
      return 0;
    }
  }
  return 1;
}

static double check_result(arglist *al, double result) {
  int i = 0, n = 0;
  if (gsl_isnan(result)) {
    eval_error(al);
    return 0;
  }
  for (i = 0; i < al->n; ++i) {
    if (gsl_isnan(al->ra[i])) {
      eval_error(al);
      return 0;
    }
  }
  if (al->derivs) {
    for (i = 0; i < al->n; ++i) {
      if (gsl_isnan(al->derivs[i])) {
        eval_error_with_suffix(al, "'");
        return 0;
      }
    }
    if (al->hes) {
      for (i = 0, n = al->n * (al->n + 1) / 2; i < n; ++i) {
        if (gsl_isnan(al->hes[i])) {
          eval_error_with_suffix(al, "''");
          return 0;
        }
      }
    }
  }
  return result;
}

/* Flags for check_bessel_args */
enum {
  DERIV_INT_MIN = 1 /* Derivative can be computed for n = INT_MIN */
};

/*
 * Checks whether the first argument is constant and reports error if not.
 * Returns 1 iff the first argument is constant.
 */
static int check_const_arg(arglist *al, unsigned index, const char *name) {
  if (al->dig && al->dig[index])
    return 1;
  /* Derivative information is requested, so the argument is not constant. */
  error(al, "argument '%s' is not constant", name);
  return 0;
}

/* Checks if the argument with the specified index is representable as int. */
static int check_int_arg(arglist *al, unsigned index, const char *name) {
  double arg = al->ra[index];
  if ((int)arg != arg) {
    error(al, "argument '%s' can't be represented as int, %s = %g",
        name, name, arg);
    return 0;
  }
  return al->derivs ? check_const_arg(al, index, name) : 1;
}

/*
 * Checks the arguments of a zero function such as gsl_sf_airy_Ai_scaled:
 * - argument with the specified index should be representable as unsigned int
 * - al->derivs should be null
 */
static int check_zero_func_args(arglist *al, unsigned s_index) {
  double arg = al->ra[s_index];
  if ((unsigned)arg != arg) {
    error(al, "argument 's' can't be represented as unsigned int, s = %g", arg);
    return 0;
  }
  if (al->derivs) {
    if (check_const_arg(al, s_index, "s"))
      error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return 1;
}

/* Checks the arguments of a Bessel function. */
static int check_bessel_args(arglist *al, int flags, const char *arg_name) {
  int n = al->ra[0];
  if (!check_int_arg(al, 0, arg_name))
    return 0;
  if (al->derivs) {
    int deriv_min = INT_MIN + ((flags & DERIV_INT_MIN) != 0 ? 0 : 1);
    if (!al->dig || !al->dig[0]) {
      /* Can't compute derivative with respect to an integer argument. */
      error(al, "argument '%s' is not constant", arg_name);
      return 0;
    }
    if ((al->hes && !check_deriv_arg(al, n, INT_MIN + 2, INT_MAX - 2)) ||
        !check_deriv_arg(al, n, deriv_min, INT_MAX - 1))
      return 0;
  }
  return 1;
}

static double amplgsl_log1p(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    double deriv = *al->derivs = 1 / (x + 1);
    if (al->hes)
      *al->hes = -deriv * deriv;
  }
  return check_result(al, gsl_log1p(x));
}

static double amplgsl_expm1(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    double deriv = *al->derivs = exp(x);
    if (al->hes)
      *al->hes = deriv;
  }
  return check_result(al, gsl_expm1(x));
}

static double amplgsl_hypot(arglist *al) {
  double x = al->ra[0];
  double y = al->ra[1];
  double hypot = gsl_hypot(x, y);
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
  return check_result(al, hypot);
}

static double amplgsl_hypot3(arglist *al) {
  double x = al->ra[0];
  double y = al->ra[1];
  double z = al->ra[2];
  double hypot = gsl_hypot3(x, y, z);
  if (al->derivs) {
    real *derivs = al->derivs;
    derivs[0] = x / hypot;
    derivs[1] = y / hypot;
    derivs[2] = z / hypot;
    if (al->hes) {
      real *hes = al->hes;
      double dx2 = derivs[0] * derivs[0];
      double dy2 = derivs[1] * derivs[1];
      double dz2 = derivs[2] * derivs[2];
      hes[0] =  (dy2 + dz2) / hypot;
      hes[1] = -derivs[0] * derivs[1] / hypot;
      hes[2] = -derivs[0] * derivs[2] / hypot;
      hes[3] =  (dx2 + dz2) / hypot;
      hes[4] = -derivs[1] * derivs[2] / hypot;
      hes[5] =  (dx2 + dy2) / hypot;
    }
  }
  return check_result(al, hypot);
}

static double amplgsl_sf_airy_Ai(arglist *al) {
  double x = al->ra[0];
  double value = gsl_sf_airy_Ai(x, GSL_PREC_DOUBLE);
  if (al->derivs) {
    *al->derivs = gsl_sf_airy_Ai_deriv(x, GSL_PREC_DOUBLE);
    if (al->hes)
      *al->hes = x * value;
  }
  return check_result(al, value);
}

static double amplgsl_sf_airy_Bi(arglist *al) {
  double x = al->ra[0];
  double value = gsl_sf_airy_Bi(x, GSL_PREC_DOUBLE);
  if (al->derivs) {
    *al->derivs = gsl_sf_airy_Bi_deriv(x, GSL_PREC_DOUBLE);
    if (al->hes)
      *al->hes = x * value;
  }
  return check_result(al, value);
}

static double amplgsl_sf_airy_Ai_scaled(arglist *al) {
  double x = al->ra[0];
  double value = gsl_sf_airy_Ai_scaled(x, GSL_PREC_DOUBLE);
  if (al->derivs) {
    if (x > 0) {
      double sqrtx = sqrt(x);
      *al->derivs = gsl_sf_airy_Ai_deriv_scaled(x, GSL_PREC_DOUBLE) +
          sqrtx * gsl_sf_airy_Ai_scaled(x, GSL_PREC_DOUBLE);
      if (al->hes)
        *al->hes = (value + 4 * x * *al->derivs) / (2 * sqrtx);
    } else {
      *al->derivs = gsl_sf_airy_Ai_deriv(x, GSL_PREC_DOUBLE);
      if (al->hes) {
        /* Return NaN for x = 0 since the right derivative is infinity
           and the left derivative is 0 at this point. */
        *al->hes = x != 0 ? x * value : GSL_NAN;
      }
    }
  }
  return check_result(al, value);
}

static double amplgsl_sf_airy_Bi_scaled(arglist *al) {
  double x = al->ra[0];
  double value = gsl_sf_airy_Bi_scaled(x, GSL_PREC_DOUBLE);
  if (al->derivs) {
    if (x > 0) {
      double sqrtx = sqrt(x);
      *al->derivs = gsl_sf_airy_Bi_deriv_scaled(x, GSL_PREC_DOUBLE) -
          sqrtx * gsl_sf_airy_Bi_scaled(x, GSL_PREC_DOUBLE);
      if (al->hes)
        *al->hes = -(value + 4 * x * *al->derivs) / (2 * sqrtx);
    } else {
      *al->derivs = gsl_sf_airy_Bi_deriv(x, GSL_PREC_DOUBLE);
      if (al->hes) {
        /* Return NaN for x = 0 since the right derivative is -infinity
           and left derivative is 0 at this point. */
        *al->hes = x != 0 ? x * value : GSL_NAN;
      }
    }
  }
  return check_result(al, value);
}

static double amplgsl_sf_airy_zero_Ai(arglist *al) {
  if (!check_zero_func_args(al, 0))
    return 0;
  return check_result(al, gsl_sf_airy_zero_Ai(al->ra[0]));
}

static double amplgsl_sf_airy_zero_Bi(arglist *al) {
  if (!check_zero_func_args(al, 0))
    return 0;
  return check_result(al, gsl_sf_airy_zero_Bi(al->ra[0]));
}

static double amplgsl_sf_airy_zero_Ai_deriv(arglist *al) {
  if (!check_zero_func_args(al, 0))
    return 0;
  return check_result(al, gsl_sf_airy_zero_Ai_deriv(al->ra[0]));
}

static double amplgsl_sf_airy_zero_Bi_deriv(arglist *al) {
  if (!check_zero_func_args(al, 0))
    return 0;
  return check_result(al, gsl_sf_airy_zero_Bi_deriv(al->ra[0]));
}

static double amplgsl_sf_bessel_J0(arglist *al) {
  double x = al->ra[0];
  double j0 = gsl_sf_bessel_J0(x);
  if (al->derivs) {
    *al->derivs = -gsl_sf_bessel_J1(x);
    if (al->hes)
      *al->hes = 0.5 * (gsl_sf_bessel_Jn(2, x) - j0);
  }
  return check_result(al, j0);
}

static double amplgsl_sf_bessel_J1(arglist *al) {
  double x = al->ra[0];
  double j1 = gsl_sf_bessel_J1(x);
  if (al->derivs) {
    *al->derivs = 0.5 * (gsl_sf_bessel_J0(x) - gsl_sf_bessel_Jn(2, x));
    if (al->hes)
      *al->hes = 0.25 * (gsl_sf_bessel_Jn(3, x) - 3 * j1);
  }
  return check_result(al, j1);
}

static double amplgsl_sf_bessel_Jn(arglist *al) {
  int n = al->ra[0];
  double x = al->ra[1];
  double jn = 0;
  if (!check_bessel_args(al, 0, "n"))
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
  return check_result(al, jn);
}

static double amplgsl_sf_bessel_Y0(arglist *al) {
  double x = al->ra[0];
  double y0 = gsl_sf_bessel_Y0(x);
  if (al->derivs) {
    *al->derivs = -gsl_sf_bessel_Y1(x);
    if (al->hes)
      *al->hes = 0.5 * (gsl_sf_bessel_Yn(2, x) - y0);
  }
  return check_result(al, y0);
}

static double amplgsl_sf_bessel_Y1(arglist *al) {
  double x = al->ra[0];
  double y1 = gsl_sf_bessel_Y1(x);
  if (al->derivs) {
    *al->derivs = 0.5 * (gsl_sf_bessel_Y0(x) - gsl_sf_bessel_Yn(2, x));
    if (al->hes)
      *al->hes = 0.25 * (gsl_sf_bessel_Yn(3, x) - 3 * y1);
  }
  return check_result(al, y1);
}

static double amplgsl_sf_bessel_Yn(arglist *al) {
  int n = al->ra[0];
  double x = al->ra[1];
  double yn = 0;
  if (!check_bessel_args(al, 0, "n"))
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
  return check_result(al, yn);
}

static double amplgsl_sf_bessel_I0(arglist *al) {
  double x = al->ra[0];
  double i0 = gsl_sf_bessel_I0(x);
  if (al->derivs) {
    *al->derivs = gsl_sf_bessel_I1(x);
    if (al->hes)
      *al->hes = 0.5 * (gsl_sf_bessel_In(2, x) + i0);
  }
  return check_result(al, i0);
}

static double amplgsl_sf_bessel_I1(arglist *al) {
  double x = al->ra[0];
  double i1 = gsl_sf_bessel_I1(x);
  if (al->derivs) {
    *al->derivs = 0.5 * (gsl_sf_bessel_I0(x) + gsl_sf_bessel_In(2, x));
    if (al->hes)
      *al->hes = 0.25 * (gsl_sf_bessel_In(3, x) + 3 * i1);
  }
  return check_result(al, i1);
}

static double amplgsl_sf_bessel_In(arglist *al) {
  int n = al->ra[0];
  double x = al->ra[1];
  double in = 0;
  if (!check_bessel_args(al, 0, "n"))
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
  return check_result(al, in);
}

static double amplgsl_sf_bessel_I0_scaled(arglist *al) {
  double x = al->ra[0];
  double i0 = gsl_sf_bessel_I0_scaled(x);
  if (al->derivs) {
    double i1 = gsl_sf_bessel_I1_scaled(x);
    *al->derivs = i1 - mul_by_sign(x, i0);
    if (al->hes) {
      *al->hes = 1.5 * i0 - 2 * fabs(x) * i1 / x +
          0.5 * gsl_sf_bessel_In_scaled(2, x);
    }
  }
  return check_result(al, i0);
}

static double amplgsl_sf_bessel_I1_scaled(arglist *al) {
  double x = al->ra[0];
  double i1 = gsl_sf_bessel_I1_scaled(x);
  if (al->derivs) {
    double i0 = gsl_sf_bessel_I0_scaled(x), i2 = gsl_sf_bessel_In_scaled(2, x);
    *al->derivs = 0.5 * i0 - mul_by_sign(x, i1) + 0.5 * i2;
    if (al->hes) {
      *al->hes = -fabs(x) * i0 / x + 1.75 * i1 - fabs(x) * i2 / x +
          0.25 * gsl_sf_bessel_In_scaled(3, x);
    }
  }
  return check_result(al, i1);
}

static double amplgsl_sf_bessel_In_scaled(arglist *al) {
  int n = al->ra[0];
  double x = al->ra[1];
  double in = 0;
  if (!check_bessel_args(al, 0, "n"))
    return 0;
  in = gsl_sf_bessel_In_scaled(n, x);
  if (al->derivs) {
    double in_minus_1 = gsl_sf_bessel_In_scaled(n - 1, x);
    double in_plus_1 = gsl_sf_bessel_In_scaled(n + 1, x);
    al->derivs[1] = 0.5 * in_minus_1 - mul_by_sign(x, in) + 0.5 * in_plus_1;
    if (al->hes) {
      al->hes[2] = 0.25 *
          (gsl_sf_bessel_In_scaled(n - 2, x) + 6 * in +
           gsl_sf_bessel_In_scaled(n + 2, x)) -
           mul_by_sign(x, in_minus_1 + in_plus_1);
    }
  }
  return check_result(al, in);
}

static double amplgsl_sf_bessel_K0(arglist *al) {
  double x = al->ra[0];
  double k0 = gsl_sf_bessel_K0(x);
  if (al->derivs) {
    *al->derivs = -gsl_sf_bessel_K1(x);
    if (al->hes)
      *al->hes = 0.5 * (gsl_sf_bessel_Kn(2, x) + k0);
  }
  return check_result(al, k0);
}

static double amplgsl_sf_bessel_K1(arglist *al) {
  double x = al->ra[0];
  double k1 = gsl_sf_bessel_K1(x);
  if (al->derivs) {
    *al->derivs = -0.5 * (gsl_sf_bessel_K0(x) + gsl_sf_bessel_Kn(2, x));
    if (al->hes)
      *al->hes = 0.25 * (gsl_sf_bessel_Kn(3, x) + 3 * k1);
  }
  return check_result(al, k1);
}

static double amplgsl_sf_bessel_Kn(arglist *al) {
  int n = al->ra[0];
  double x = al->ra[1];
  double kn = 0;
  if (!check_bessel_args(al, 0, "n"))
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
  return check_result(al, kn);
}

static double amplgsl_sf_bessel_K0_scaled(arglist *al) {
  double x = al->ra[0];
  double k0 = gsl_sf_bessel_K0_scaled(x);
  if (al->derivs) {
    double k1 = gsl_sf_bessel_K1_scaled(x);
    *al->derivs = k0 - k1;
    if (al->hes)
      *al->hes = 1.5 * k0 - 2 * k1 + 0.5 * gsl_sf_bessel_Kn_scaled(2, x);
  }
  return check_result(al, k0);
}

static double amplgsl_sf_bessel_K1_scaled(arglist *al) {
  double x = al->ra[0];
  double k1 = gsl_sf_bessel_K1_scaled(x);
  if (al->derivs) {
    double k0 = gsl_sf_bessel_K0_scaled(x), k2 = gsl_sf_bessel_Kn_scaled(2, x);
    *al->derivs = -0.5 * k0 + k1 - 0.5 * k2;
    if (al->hes)
      *al->hes = -k0 + 1.75 * k1 - k2 + 0.25 * gsl_sf_bessel_Kn_scaled(3, x);
  }
  return check_result(al, k1);
}

static double amplgsl_sf_bessel_Kn_scaled(arglist *al) {
  int n = al->ra[0];
  double x = al->ra[1];
  double kn = 0;
  if (!check_bessel_args(al, 0, "n"))
    return 0;
  kn = gsl_sf_bessel_Kn_scaled(n, x);
  if (al->derivs) {
    double kn_minus_1 = gsl_sf_bessel_Kn_scaled(n - 1, x);
    double kn_plus_1 = gsl_sf_bessel_Kn_scaled(n + 1, x);
    al->derivs[1] = -0.5 * (kn_minus_1 - 2 * kn + kn_plus_1);
    if (al->hes) {
      al->hes[2] = 0.25 *
          (gsl_sf_bessel_Kn_scaled(n - 2, x) - 4 * kn_minus_1 + 6 * kn -
              4 * kn_plus_1 + gsl_sf_bessel_Kn_scaled(n + 2, x));
    }
  }
  return check_result(al, kn);
}

static double amplgsl_sf_bessel_j0(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = x != 0 ? (x * cos(x) - sin(x)) / gsl_pow_2(x) : 0;
    if (al->hes)
      *al->hes = x != 0 ?
          ((2 - gsl_pow_2(x)) * sin(x) - 2 * x * cos(x)) / gsl_pow_3(x) :
          -1.0 / 3.0;
  }
  return check_result(al, gsl_sf_bessel_j0(x));
}

static double amplgsl_sf_bessel_j1(arglist *al) {
  double x = al->ra[0];
  double j1 = gsl_sf_bessel_j1(x);
  if (al->derivs) {
    *al->derivs = x != 0 ? (sin(x) - 2 * j1) / x : 1.0 / 3.0;
    if (al->hes) {
      *al->hes = x != 0 ? (x * (gsl_pow_2(x) - 6) * cos(x) -
          3 * (gsl_pow_2(x) - 2) * sin(x)) / gsl_pow_4(x) : 0;
    }
  }
  return check_result(al, j1);
}

static double amplgsl_sf_bessel_j2(arglist *al) {
  double x = al->ra[0];
  double j2 = gsl_sf_bessel_j2(x);
  if (al->derivs) {
    *al->derivs = x != 0 ? gsl_sf_bessel_j1(x) - 3 * j2 / x : 0;
    if (al->hes) {
      *al->hes = x != 0 ? (x * (5 * gsl_pow_2(x) - 36) * cos(x) +
          (gsl_pow_4(x) - 17 * gsl_pow_2(x) + 36) * sin(x)) / gsl_pow_5(x) :
              2.0 / 15.0;
    }
  }
  return check_result(al, j2);
}

static double amplgsl_sf_bessel_jl(arglist *al) {
  int el = al->ra[0];
  double x = al->ra[1];
  double jl = 0;
  if (!check_bessel_args(al, DERIV_INT_MIN, "l"))
    return 0;
  jl = gsl_sf_bessel_jl(el, x);
  if (al->derivs) {
    double jl_plus_1 = gsl_sf_bessel_jl(el + 1, x);
    if (x == 0)
      al->derivs[1] = el == 1 ? 1.0 / 3.0 : 0;
    else
      al->derivs[1] = el * jl / x - jl_plus_1;
    if (al->hes) {
      double hes = 0;
      if (x == 0) {
        if (el == 0)
          hes = -1.0 / 3.0;
        else if (el == 2)
          hes = 2.0 / 15.0;
      } else {
        double jl_minus_1 = 0, jl_minus_2 = 0;
        if (el == 0) {
          jl_minus_1 = cos(x) / x;
          jl_minus_2 = -(cos(x) / x + sin(x)) / x;
        } else if (el == 1) {
          jl_minus_1 = gsl_sf_bessel_jl(el - 1, x);
          jl_minus_2 = cos(x) / x;
        } else {
          jl_minus_1 = gsl_sf_bessel_jl(el - 1, x);
          jl_minus_2 = gsl_sf_bessel_jl(el - 2, x);
        }
        hes = (
          gsl_pow_2(x) * jl_minus_2 -
          2 * x * jl_minus_1 -
          (2 * gsl_pow_2(x) - 3) * jl + 2 * x * jl_plus_1 +
          gsl_pow_2(x) * gsl_sf_bessel_jl(el + 2, x)) / (4 * gsl_pow_2(x));
      }
      al->hes[2] = hes;
    }
  }
  return check_result(al, jl);
}

static double amplgsl_sf_bessel_y0(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = (x * sin(x) + cos(x)) / gsl_pow_2(x);
    if (al->hes) {
      *al->hes = ((gsl_pow_2(x) - 2) * cos(x) - 2 * x * sin(x)) /
          gsl_pow_3(x);
    }
  }
  return check_result(al, gsl_sf_bessel_y0(x));
}

static double amplgsl_sf_bessel_y1(arglist *al) {
  double x = al->ra[0];
  double y1 = gsl_sf_bessel_y1(x);
  if (al->derivs) {
    *al->derivs = -(2 * y1 + cos(x)) / x;
    if (al->hes) {
      *al->hes = (x * (gsl_pow_2(x) - 6) * sin(x) +
          3 * (gsl_pow_2(x) - 2) * cos(x)) / gsl_pow_4(x);
    }
  }
  return check_result(al, y1);
}

static double amplgsl_sf_bessel_y2(arglist *al) {
  double x = al->ra[0];
  double y2 = gsl_sf_bessel_y2(x);
  if (al->derivs) {
    double y1 = gsl_sf_bessel_y1(x);
    *al->derivs = y1 - (3 * y2) / x;
    if (al->hes) {
      *al->hes = ((36 - 5 * gsl_pow_2(x)) * y1 -
          (gsl_pow_2(x) - 12) * cos(x)) / gsl_pow_3(x);
    }
  }
  return check_result(al, y2);
}

static double amplgsl_sf_bessel_yl(arglist *al) {
  int el = al->ra[0];
  double x = al->ra[1];
  double yl = 0;
  if (!check_bessel_args(al, 0, "l"))
    return 0;
  yl = gsl_sf_bessel_yl(el, x);
  if (al->derivs) {
    double yl_minus_1 = el != 0 ? gsl_sf_bessel_yl(el - 1, x) : sin(x) / x;
    double yl_plus_1 = gsl_sf_bessel_yl(el + 1, x);
    al->derivs[1] = 0.5 * (yl_minus_1 - yl / x - yl_plus_1);
    if (al->hes) {
      double yl_minus_2 = 0;
      if (el == 0)
        yl_minus_2 = (cos(x) - sin(x) / x) / x;
      else if (el == 1)
        yl_minus_2 = sin(x) / x;
      else
        yl_minus_2 = gsl_sf_bessel_yl(el - 2, x);
      al->hes[2] = (
          gsl_pow_2(x) * yl_minus_2 - 2 * x * yl_minus_1 -
          (2 * gsl_pow_2(x) - 3) * yl + 2 * x * yl_plus_1 +
          gsl_pow_2(x) * gsl_sf_bessel_yl(el + 2, x)) / (4 * gsl_pow_2(x));
    }
  }
  return check_result(al, yl);
}

static double amplgsl_sf_bessel_i0_scaled(arglist *al) {
  double x = al->ra[0];
  double i0 = gsl_sf_bessel_i0_scaled(x);
  if (al->derivs) {
    /* Contrary to the documentation, gsl_sf_bessel_i0_scaled
       implements \exp(-|x|) \sqrt{\pi}/\sqrt{2x} I_{1/2}(x)
       and not \exp(-|x|) \sqrt{\pi/(2x)} I_{1/2}(x).
       These are different since \sqrt(1/x) != \sqrt(x) for negative x. */
    double hyp_coef = exp(-fabs(x)) / x;
    double i_minus_1 = hyp_coef * cosh(x);
    double i1 = gsl_sf_bessel_i1_scaled(x);
    double coef = -(1 + 2 * fabs(x)) / x;
    *al->derivs = 0.5 * (i_minus_1 + coef * i0 + i1);
    if (al->hes) {
      coef *= 2;
      *al->hes = 0.25 * (
          hyp_coef * sinh(x) - i_minus_1 / x +
          coef * i_minus_1 +
          (3 + 6 * gsl_pow_2(x) + 4 * fabs(x)) * i0 / gsl_pow_2(x) +
          coef * i1 +
          gsl_sf_bessel_il_scaled(2, x));
    }
  }
  return check_result(al, i0);
}

static double amplgsl_sf_bessel_i1_scaled(arglist *al) {
  double x = al->ra[0];
  double i1 = gsl_sf_bessel_i1_scaled(x);
  if (al->derivs) {
    /* Contrary to the documentation, gsl_sf_bessel_i1_scaled
       implements \exp(-|x|) \sqrt{\pi}/\sqrt{2x} I_{1+1/2}(x)
       and not \exp(-|x|) \sqrt{\pi/(2x)} I_{1+1/2}(x).
       These are different since \sqrt(1/x) != \sqrt(x) for negative x. */
    double i0 = gsl_sf_bessel_i0_scaled(x);
    double i2 = gsl_sf_bessel_i2_scaled(x);
    double coef = -(1 + 2 * fabs(x)) / x;
    *al->derivs = x != 0 ? 0.5 * (i0 + coef * i1 + i2) : 1.0 / 3.0;
    if (al->hes) {
      coef *= 2;
      *al->hes = 0.25 * (
          exp(-fabs(x)) * cosh(x) / x +
          coef * i0 +
          (3 + 6 * gsl_pow_2(x) + 4 * fabs(x)) * i1 / gsl_pow_2(x) +
          coef * i2 +
          gsl_sf_bessel_il_scaled(3, x));
    }
  }
  return check_result(al, i1);
}

static double amplgsl_sf_bessel_i2_scaled(arglist *al) {
  double x = al->ra[0];
  double i2 = gsl_sf_bessel_i2_scaled(x);
  if (al->derivs) {
    /* Contrary to the documentation, gsl_sf_bessel_i2_scaled
       implements \exp(-|x|) \sqrt{\pi}/\sqrt{2x} I_{2+1/2}(x)
       and not \exp(-|x|) \sqrt{\pi/(2x)} I_{2+1/2}(x).
       These are different since \sqrt(1/x) != \sqrt(x) for negative x. */
    double i1 = gsl_sf_bessel_i1_scaled(x);
    double i3 = gsl_sf_bessel_il_scaled(3, x);
    double coef = -(1 + 2 * fabs(x)) / x;
    *al->derivs = x != 0 ? 0.5 * (i1 + coef * i2 + i3) : 0;
    if (al->hes) {
      coef *= 2;
      *al->hes = x != 0 ? 0.25 * (
          gsl_sf_bessel_i0_scaled(x) +
          coef * i1 +
          (3 + 6 * gsl_pow_2(x) + 4 * fabs(x)) * i2 / gsl_pow_2(x) +
          coef * i3 +
          gsl_sf_bessel_il_scaled(4, x)) : 2.0 / 15.0;
    }
  }
  return check_result(al, i2);
}

static double amplgsl_sf_bessel_il_scaled(arglist *al) {
  int el = al->ra[0];
  double x = al->ra[1];
  double il = 0;
  if (!check_bessel_args(al, 0, "l"))
    return 0;
  il = gsl_sf_bessel_il_scaled(el, x);
  if (al->derivs) {
    double il_minus_1 = el != 0 ?
        gsl_sf_bessel_il_scaled(el - 1, x) : exp(-fabs(x)) * cosh(x) / x;
    double il_plus_1 = gsl_sf_bessel_il_scaled(el + 1, x);
    double coef = -(1 + 2 * fabs(x)) / x;
    double deriv = GSL_NAN;
    if (x == 0) {
      /* If el <= 0, keep deriv equal to GSL_NAN. */
      if (el == 1)
        deriv = 1.0 / 3.0;
      else if (el > 1)
        deriv = 0;
    } else {
      deriv = 0.5 * (il_minus_1 + coef * il + il_plus_1);
    }
    al->derivs[1] = deriv;
    if (al->hes) {
      double hes = GSL_NAN;
      double il_minus_2 = 0;
      if (el == 0)
        il_minus_2 = (exp(-fabs(x)) * (sinh(x) - cosh(x) / x)) / x;
      else if (el == 1)
        il_minus_2 = exp(-fabs(x)) * cosh(x) / x;
      else
        il_minus_2 = gsl_sf_bessel_il_scaled(el - 2, x);
      coef *= 2;
      if (x == 0) {
        /* If el == 1 or el < 0, keep hes equal to GSL_NAN. */
        if (el == 0)
          hes = 4.0 / 3.0;
        else if (el == 2)
          hes = 2.0 / 15.0;
        else if (el > 2)
          hes = 0;
      } else {
        hes = 0.25 * (
          il_minus_2 +
          coef * il_minus_1 +
          (3 + 4 * fabs(x) + 6 * gsl_pow_2(x)) * il / gsl_pow_2(x) +
          coef * il_plus_1 +
          gsl_sf_bessel_il_scaled(el + 2, x));
      }
      al->hes[2] = hes;
    }
  }
  return check_result(al, il);
}

static double amplgsl_sf_bessel_k0_scaled(arglist *al) {
  double x = al->ra[0];
  double k0 = gsl_sf_bessel_k0_scaled(x);
  if (al->derivs) {
    double pi_sqrt_inv_x = M_PI * sqrt(1 / x);
    *al->derivs = -pi_sqrt_inv_x / (2 * pow(x, 1.5));
    if (al->hes)
      *al->hes = pi_sqrt_inv_x / pow(x, 2.5);
  }
  return check_result(al, k0);
}

static double amplgsl_sf_bessel_k1_scaled(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    double pi_sqrt_inv_x = M_PI * sqrt(1 / x);
    *al->derivs = -(pi_sqrt_inv_x * (x + 2)) / (2 * pow(x, 2.5));
    if (al->hes)
      *al->hes = (pi_sqrt_inv_x * (x + 3)) / pow(x, 3.5);
  }
  return check_result(al, gsl_sf_bessel_k1_scaled(x));
}

static double amplgsl_sf_bessel_k2_scaled(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    double pi_sqrt_inv_x = M_PI * sqrt(1 / x);
    *al->derivs = -pi_sqrt_inv_x * (x + 3) * (x + 3) / (2 * pow(x, 3.5));
    if (al->hes)
      *al->hes = pi_sqrt_inv_x * (x * x + 9 * x + 18) / pow(x, 4.5);
  }
  return check_result(al, gsl_sf_bessel_k2_scaled(x));
}

static double amplgsl_sf_bessel_kl_scaled(arglist *al) {
  int el = al->ra[0];
  double x = al->ra[1];
  double kl = 0;
  if (!check_bessel_args(al, 0, "l"))
    return 0;
  kl = gsl_sf_bessel_kl_scaled(el, x);
  if (al->derivs) {
    double kl_minus_1 = el != 0 ?
        gsl_sf_bessel_kl_scaled(el - 1, x) : M_PI_2 / x;
    double kl_plus_1 = gsl_sf_bessel_kl_scaled(el + 1, x);
    double coef = (1 - 2 * x) / x;
    al->derivs[1] = -0.5 * (kl_minus_1 + coef * kl + kl_plus_1);
    if (al->hes) {
      double kl_minus_2 = 0;
      if (el == 0)
        kl_minus_2 = M_PI_2 * (1 / x + 1) / x;
      else if (el == 1)
        kl_minus_2 = M_PI_2 / x;
      else
        kl_minus_2 = gsl_sf_bessel_kl_scaled(el - 2, x);
      coef *= 2;
      al->hes[2] = 0.25 * (
          kl_minus_2 +
          coef * kl_minus_1 +
          (3 - 4 * x + 6 * gsl_pow_2(x)) * kl / gsl_pow_2(x) +
          coef * kl_plus_1 +
          gsl_sf_bessel_kl_scaled(el + 2, x));
    }
  }
  return check_result(al, kl);
}

static double amplgsl_sf_bessel_Jnu(arglist *al) {
  double n = al->ra[0];
  double x = al->ra[1];
  double jn = gsl_sf_bessel_Jnu(n, x);
  if (al->derivs) {
    if (!check_const_arg(al, 0, "nu"))
      return 0;
    al->derivs[1] = 0.5 *
        (gsl_sf_bessel_Jnu(n - 1, x) - gsl_sf_bessel_Jnu(n + 1, x));
    if (al->hes) {
      al->hes[2] = 0.25 *
          (gsl_sf_bessel_Jnu(n - 2, x) - 2 * jn + gsl_sf_bessel_Jnu(n + 2, x));
    }
  }
  return check_result(al, jn);
}

static double amplgsl_sf_bessel_Ynu(arglist *al) {
  double n = al->ra[0];
  double x = al->ra[1];
  double yn = gsl_sf_bessel_Ynu(n, x);
  if (al->derivs) {
    if (!check_const_arg(al, 0, "nu"))
      return 0;
    al->derivs[1] = 0.5 *
        (gsl_sf_bessel_Ynu(n - 1, x) - gsl_sf_bessel_Ynu(n + 1, x));
    if (al->hes) {
      al->hes[2] = 0.25 *
          (gsl_sf_bessel_Ynu(n - 2, x) - 2 * yn + gsl_sf_bessel_Ynu(n + 2, x));
    }
  }
  return check_result(al, yn);
}

static double amplgsl_sf_bessel_Inu(arglist *al) {
  double n = al->ra[0];
  double x = al->ra[1];
  double in = gsl_sf_bessel_Inu(n, x);
  if (al->derivs) {
    if (!check_const_arg(al, 0, "nu"))
      return 0;
    al->derivs[1] = 0.5 *
        (gsl_sf_bessel_Inu(n - 1, x) + gsl_sf_bessel_Inu(n + 1, x));
    if (al->hes) {
      al->hes[2] = 0.25 *
          (gsl_sf_bessel_Inu(n - 2, x) + 2 * in + gsl_sf_bessel_Inu(n + 2, x));
    }
  }
  return check_result(al, in);
}

static double amplgsl_sf_bessel_Inu_scaled(arglist *al) {
  double n = al->ra[0];
  double x = al->ra[1];
  double in = gsl_sf_bessel_Inu_scaled(n, x);
  if (al->derivs) {
    double in_minus_1 = 0, in_plus_1 = 0;
    if (!check_const_arg(al, 0, "nu"))
      return 0;
    in_minus_1 = gsl_sf_bessel_Inu_scaled(n - 1, x);
    in_plus_1 = gsl_sf_bessel_Inu_scaled(n + 1, x);
    al->derivs[1] = 0.5 * in_minus_1 - mul_by_sign(x, in) + 0.5 * in_plus_1;
    if (al->hes) {
      al->hes[2] = 0.25 *
          (gsl_sf_bessel_Inu_scaled(n - 2, x) + 6 * in +
              gsl_sf_bessel_Inu_scaled(n + 2, x)) -
              mul_by_sign(x, in_minus_1 + in_plus_1);
    }
  }
  return check_result(al, in);
}

static double amplgsl_sf_bessel_Knu(arglist *al) {
  double n = al->ra[0];
  double x = al->ra[1];
  double kn = gsl_sf_bessel_Knu(n, x);
  if (al->derivs) {
    if (!check_const_arg(al, 0, "nu"))
      return 0;
    al->derivs[1] = -0.5 *
        (gsl_sf_bessel_Knu(n - 1, x) + gsl_sf_bessel_Knu(n + 1, x));
    if (al->hes) {
      al->hes[2] = 0.25 *
          (gsl_sf_bessel_Knu(n - 2, x) + 2 * kn + gsl_sf_bessel_Knu(n + 2, x));
    }
  }
  return check_result(al, kn);
}

static double amplgsl_sf_bessel_lnKnu(arglist *al) {
  double n = al->ra[0];
  double x = al->ra[1];
  if (al->derivs) {
    double kn = 0, kn_minus_1_plus_1 = 0;
    if (!check_const_arg(al, 0, "nu"))
      return 0;
    kn = gsl_sf_bessel_Knu(n, x);
    kn_minus_1_plus_1 =
        gsl_sf_bessel_Knu(n - 1, x) + gsl_sf_bessel_Knu(n + 1, x);
    al->derivs[1] = -0.5 * kn_minus_1_plus_1 / kn;
    if (al->hes) {
      al->hes[2] = 0.25 *
          (kn * (gsl_sf_bessel_Knu(n - 2, x) + 2 * kn +
          gsl_sf_bessel_Knu(n + 2, x)) -
              kn_minus_1_plus_1 * kn_minus_1_plus_1) / (kn * kn);
    }
  }
  return check_result(al, gsl_sf_bessel_lnKnu(n, x));
}

static double amplgsl_sf_bessel_Knu_scaled(arglist *al) {
  double n = al->ra[0];
  double x = al->ra[1];
  double kn = gsl_sf_bessel_Knu_scaled(n, x);
  if (al->derivs) {
    double kn_minus_1 = 0, kn_plus_1 = 0;
    if (!check_const_arg(al, 0, "nu"))
      return 0;
    kn_minus_1 = gsl_sf_bessel_Knu_scaled(n - 1, x);
    kn_plus_1 = gsl_sf_bessel_Knu_scaled(n + 1, x);
    al->derivs[1] = -0.5 * (kn_minus_1 - 2 * kn + kn_plus_1);
    if (al->hes) {
      al->hes[2] = 0.25 *
          (gsl_sf_bessel_Knu_scaled(n - 2, x) - 4 * kn_minus_1 + 6 * kn -
              4 * kn_plus_1 + gsl_sf_bessel_Knu_scaled(n + 2, x));
    }
  }
  return check_result(al, kn);
}

static double amplgsl_sf_bessel_zero_J0(arglist *al) {
  if (!check_zero_func_args(al, 0))
    return 0;
  return check_result(al, gsl_sf_bessel_zero_J0(al->ra[0]));
}

static double amplgsl_sf_bessel_zero_J1(arglist *al) {
  if (!check_zero_func_args(al, 0))
    return 0;
  return check_result(al, gsl_sf_bessel_zero_J1(al->ra[0]));
}

static double amplgsl_sf_bessel_zero_Jnu(arglist *al) {
  double nu = al->ra[0];
  if (!check_zero_func_args(al, 1))
    return 0;
  return check_result(al, gsl_sf_bessel_zero_Jnu(nu, al->ra[1]));
}

static double amplgsl_sf_clausen(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = -log(2 * sin(0.5 * fabs(x)));
    if (al->hes)
      *al->hes = fmod(x, M_PI) != 0 ? -0.5 * tan(M_PI_2 - 0.5 * x) : GSL_NAN;
  }
  return check_result(al, gsl_sf_clausen(x));
}

static double amplgsl_sf_hydrogenicR_1(arglist *al) {
  double Z = al->ra[0], r = al->ra[1];
  if (al->derivs) {
    real *derivs = al->derivs;
    double exp_minusZr = exp(-Z * r);
    derivs[0] = sqrt(Z) * exp_minusZr * (3 - 2 * r * Z);
    derivs[1] = -2 * pow(Z, 2.5) * exp_minusZr;
    if (al->hes) {
      real *hes = al->hes;
      hes[0] = (exp_minusZr * (4 * gsl_pow_2(r * Z) - 12 * r * Z + 3)) /
          (2 *sqrt(Z));
      hes[1] = pow(Z, 1.5) * exp_minusZr * (2 * r * Z - 5);
      hes[2] = 2 * pow(Z, 3.5) * exp_minusZr;
    }
  }
  return check_result(al, gsl_sf_hydrogenicR_1(Z, r));
}

static double amplgsl_sf_hydrogenicR(arglist *al) {
  if (!check_int_arg(al, 0, "n") || !check_int_arg(al, 1, "l"))
    return 0;
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al, gsl_sf_hydrogenicR(
      al->ra[0], al->ra[1], al->ra[2], al->ra[3]));
}

static double amplgsl_sf_coulomb_CL(arglist *al) {
  gsl_sf_result result = {0};
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  if (gsl_sf_coulomb_CL_e(al->ra[0], al->ra[1], &result)) {
    eval_error(al);
    return 0;
  }
  return check_result(al, result.val);
}

static int check_coupling_args(arglist *al, const char *const* arg_names) {
  unsigned i = 0, n_args = al->n;
  for (; i < n_args; ++i) {
    if (!check_int_arg(al, i, arg_names[i]))
      return 0;
  }
  return 1;
}

static double amplgsl_sf_coupling_3j(arglist *al) {
  static const char *ARG_NAMES[] = {
      "two_ja", "two_jb", "two_jc",
      "two_ma", "two_mb", "two_mc"
  };
  if (!check_coupling_args(al, ARG_NAMES))
    return 0;
  return check_result(al, gsl_sf_coupling_3j(
      al->ra[0], al->ra[1], al->ra[2], al->ra[3], al->ra[4], al->ra[5]));
}

static double amplgsl_sf_coupling_6j(arglist *al) {
  static const char *ARG_NAMES[] = {
      "two_ja", "two_jb", "two_jc",
      "two_jd", "two_je", "two_jf"
  };
  if (!check_coupling_args(al, ARG_NAMES))
    return 0;
  return check_result(al, gsl_sf_coupling_6j(
      al->ra[0], al->ra[1], al->ra[2], al->ra[3], al->ra[4], al->ra[5]));
}

static double amplgsl_sf_coupling_9j(arglist *al) {
  static const char *ARG_NAMES[] = {
      "two_ja", "two_jb", "two_jc",
      "two_jd", "two_je", "two_jf",
      "two_jg", "two_jh", "two_ji"
  };
  if (!check_coupling_args(al, ARG_NAMES))
    return 0;
  return check_result(al, gsl_sf_coupling_9j(
      al->ra[0], al->ra[1], al->ra[2], al->ra[3], al->ra[4], al->ra[5],
      al->ra[6], al->ra[7], al->ra[8]));
}

static double amplgsl_sf_dawson(arglist *al) {
  double x = al->ra[0];
  double f = gsl_sf_dawson(x);
  if (al->derivs) {
    double deriv = *al->derivs = 1 - 2 * x * f;
    if (al->hes)
      *al->hes = - 2 * (f + x * deriv);
  }
  return check_result(al, f);
}

/* Values of the right derivatives of the Debye functions at 0. */
static const double DEBYE_DERIV_AT_0[] = {
    -1.0 / 4.0, -1.0 / 3.0,  -3.0 / 8.0,
    -2.0 / 5.0, -5.0 / 12.0, -3.0 / 7.0
};

/* Values of the second right derivatives of the Debye functions at 0. */
static const double DEBYE_DERIV2_AT_0[] = {
    1.0 / 18.0, 1.0 / 12.0, 1.0 / 10.0,
    1.0 / 9.0,  5.0 / 42.0, 1.0 / 8.0
};

static double debye(arglist *al, int n, double (*func)(double)) {
  double x = al->ra[0];
  double f = func(x);
  if (al->derivs) {
    double exp_x = exp(x);
    double deriv = *al->derivs = x != 0 ?
        n * (1 / (exp_x - 1) - f / x) : DEBYE_DERIV_AT_0[n - 1];
    if (al->hes) {
      *al->hes = x != 0 ? n * (-exp_x / gsl_pow_2(exp_x - 1) +
          f / gsl_pow_2(x) - deriv / x) : DEBYE_DERIV2_AT_0[n - 1];
    }
  }
  return check_result(al, f);
}

#define DEBYE(n) \
  static double amplgsl_sf_debye_##n(arglist *al) { \
    return debye(al, n, gsl_sf_debye_##n); \
  }

DEBYE(1)
DEBYE(2)
DEBYE(3)
DEBYE(4)
DEBYE(5)
DEBYE(6)

static double amplgsl_sf_dilog(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    double deriv = 0;
    if (x == 0) {
      deriv = 1;
    } else if (x == 1) {
      deriv = GSL_POSINF;
    } else {
      gsl_complex log = gsl_complex_log(gsl_complex_rect(1 - x, 0));
      deriv = -GSL_REAL(log) / x;
    }
    *al->derivs = deriv;
    if (al->hes)
      *al->hes = x != 0 ? (1 / (1 - x) - deriv) / x : 0.5;
  }
  return check_result(al, gsl_sf_dilog(x));
}

static double amplgsl_sf_ellint_Kcomp(arglist *al) {
  double k = al->ra[0];
  double kcomp = gsl_sf_ellint_Kcomp(k, GSL_PREC_DOUBLE);
  if (al->derivs) {
    double ecomp = gsl_sf_ellint_Ecomp(k, GSL_PREC_DOUBLE);
    double divisor = k * (1 - k * k);
    *al->derivs = k != 0 ? ecomp / divisor - kcomp / k : 0;
    if (al->hes) {
      *al->hes = k != 0 ?
          ((2 * gsl_pow_4(k) - 3 * gsl_pow_2(k) + 1) * kcomp +
          (3 * gsl_pow_2(k) - 1) * ecomp) / gsl_pow_2(divisor) : M_PI_4;
    }
  }
  return check_result(al, kcomp);
}

static double amplgsl_sf_ellint_Ecomp(arglist *al) {
  double k = al->ra[0];
  double ecomp = gsl_sf_ellint_Ecomp(k, GSL_PREC_DOUBLE);
  if (al->derivs) {
    double kcomp = gsl_sf_ellint_Kcomp(k, GSL_PREC_DOUBLE);
    *al->derivs = k != 0 ? (ecomp - kcomp) / k : 0;
    if (al->hes) {
      *al->hes = k != 0 ?
          ((k * k - 1) * kcomp + ecomp) / (k * k * (k * k - 1)) : -M_PI_4;
    }
  }
  return check_result(al, ecomp);
}

static double amplgsl_sf_ellint_Pcomp(arglist *al) {
  double k = al->ra[0], n = al->ra[1];
  double pcomp = gsl_sf_ellint_Pcomp(k, n, GSL_PREC_DOUBLE);
  if (al->derivs) {
    double ecomp = gsl_sf_ellint_Ecomp(k, GSL_PREC_DOUBLE);
    double kcomp = gsl_sf_ellint_Kcomp(k, GSL_PREC_DOUBLE);
    double divisor = (k * k - 1) * (k * k + n);
    if (k != 0 || n != 0) {
      al->derivs[0] = -k * ((k * k - 1) * pcomp + ecomp) / divisor;
      if (n != 0) {
        al->derivs[1] = (-kcomp * (k * k + n) +
          (k * k - n * n) * pcomp + n * ecomp) /
          (2 * n * (n + 1) * (k * k + n));
      } else {
        al->derivs[1] =
            -(4 * kcomp + M_PI * k * k * gsl_sf_hyperg_2F1(0.5, 1.5, 2, k * k) -
                4 * ecomp) / (8 * k * k);
      }
    } else {
      al->derivs[0] = 0;
      al->derivs[1] = -M_PI_4;
    }
    if (al->hes) {
      if (k != 0 || n != 0) {
        al->hes[0] = ((k * k - 1) * (kcomp * (k * k + n) +
            (k * k - 1) * (2 * k * k - n) * pcomp) +
            (3 * gsl_pow_4(k) - k * k + 2 * n) * ecomp) / gsl_pow_2(divisor);
        al->hes[1] = (k * ((k * k - 1) * (kcomp * (k * k + n) +
            (n * (3 * n + 2) - k * k) * pcomp) +
            n * (-k * k + 2 * n + 3) * ecomp)) /
            (2 * divisor * n * (n + 1) * (k * k + n));
        al->hes[2] = (kcomp * (gsl_pow_4(k) * (4 * n + 1) +
            3 * k * k * n * (3 * n + 1) + n * n * (5 * n + 2)) +
            n * (k * k * (1 - 2 * n) - n * (5 * n + 2)) * ecomp -
            (gsl_pow_4(k) * (4 * n + 1) + 2 * k * k * n * (5 * n + 2) -
                3 * gsl_pow_4(n)) * pcomp) /
                (4 * gsl_pow_2(n * (n + 1) * (k * k + n)));
      } else {
        al->hes[0] = M_PI_4;
        al->hes[1] = 0;
        al->hes[2] = 3 * M_PI / 8;
      }
    }
  }
  return check_result(al, pcomp);
}

static double amplgsl_sf_ellint_F(arglist *al) {
  double phi = al->ra[0], k = al->ra[1];
  double f = gsl_sf_ellint_F(phi, k, GSL_PREC_DOUBLE);
  if (al->derivs) {
    double e = gsl_sf_ellint_E(phi, k, GSL_PREC_DOUBLE);
    al->derivs[0] = 1 / sqrt(1 - gsl_pow_2(k * sin(phi)));
    if (!al->dig || !al->dig[1]) {
      if (k == 0) {
        al->derivs[1] = 0;
      } else if (fabs(k) == 1) {
        double sec_phi = 1 / cos(phi);
        al->derivs[1] = 0.5 * k * (atanh(sin(phi)) -
           sec_phi * sec_phi * ((1 + cos(2 * phi)) * log(sec_phi + tan(phi)) -
               sin(phi)));
      } else {
        al->derivs[1] = (e + (k * k - 1) * f -
            (k * k * cos(phi) * sin(phi)) / sqrt(1 - gsl_pow_2(k * sin(phi)))) /
            (k - gsl_pow_3(k));
      }
    }
    if (al->hes) {
      double k2 = k * k;
      al->hes[0] = (k2 * sin(phi) * cos(phi)) /
          pow(1 - gsl_pow_2(k * sin(phi)), 1.5);
      al->hes[1] = k * gsl_pow_2(sin(phi)) /
          pow(1 - gsl_pow_2(k * sin(phi)), 1.5);
      if (k == 0) {
        al->hes[2] = 0.5 * (phi - cos(phi) * sin(phi));
      } else if (fabs(k) == 1) {
        double sec_phi = 1 / cos(phi);
        al->hes[2] = sec_phi * sec_phi *
            (atanh(sin(phi)) * (2 - 46 * cos(2 * phi)) +
                8 * (1 + 7 * cos(2 * phi)) * log(sec_phi + tan(phi)) +
                sec_phi * (-11 * sec_phi * sin(3 * phi) + 13 * tan(phi))) / 32;
      } else {
        /* sub1 and sub2 are just common subexpressions */
        double sub1 = 1 - 3 * k2;
        double sub2 = M_SQRT2 * pow(2 - k2 + k2 * cos(2 * phi), 1.5);
        al->hes[2] = -(sub1 * sub2 * e - (sub1 + 2 * gsl_pow_4(k)) * sub2 * f +
          4 * gsl_pow_4(k) * (sub1 * cos(phi) * gsl_pow_3(sin(phi)) +
              sin(2 * phi))) / (gsl_pow_2(k * (k2 - 1)) * sub2);
      }
    }
  }
  return check_result(al, f);
}

static double amplgsl_sf_ellint_E(arglist *al) {
  double phi = al->ra[0], k = al->ra[1];
  double e = gsl_sf_ellint_E(phi, k, GSL_PREC_DOUBLE);
  if (al->derivs) {
    double f = gsl_sf_ellint_F(phi, k, GSL_PREC_DOUBLE);
    double d_phi = al->derivs[0] = sqrt(1 - gsl_pow_2(k * sin(phi)));
    al->derivs[1] = k != 0 ? (e - f) / k : 0;
    if (al->hes) {
      double k2 = k * k;
      al->hes[0] = -k2 * cos(phi) * sin(phi) / d_phi;
      al->hes[1] = -k * gsl_pow_2(sin(phi)) / d_phi;
      if (k == 0) {
        al->hes[2] = -0.5 * phi + 0.25 * sin(2 * phi);
      } else if (fabs(k) == 1) {
        double sec_phi = 1 / cos(phi), tan_phi = tan(phi);
        al->hes[2] = -0.5 * atanh(sin(phi)) + log(sec_phi + tan_phi) -
            0.5 * sec_phi * tan_phi;
      } else {
        al->hes[2] = ((k2 - 1) * sqrt(4 - 2 * k2 + 2 * k2 * cos(2 * phi)) * f +
          2 * e * d_phi - k2 * sin(2 * phi)) / (2 * k2 * (k2 - 1) * d_phi);
      }
    }
  }
  return check_result(al, e);
}

static double amplgsl_sf_ellint_P(arglist *al) {
  double phi = al->ra[0], k = al->ra[1], n = al->ra[2];
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al, gsl_sf_ellint_P(phi, k, n, GSL_PREC_DOUBLE));
}

static double amplgsl_sf_ellint_D(arglist *al) {
  double phi = al->ra[0], k = al->ra[1], n = al->ra[2];
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al, gsl_sf_ellint_D(phi, k, n, GSL_PREC_DOUBLE));
}

static double amplgsl_sf_ellint_RC(arglist *al) {
  double x = al->ra[0], y = al->ra[1];
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al, gsl_sf_ellint_RC(x, y, GSL_PREC_DOUBLE));
}

static double amplgsl_sf_ellint_RD(arglist *al) {
  double x = al->ra[0], y = al->ra[1], z = al->ra[2];
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al, gsl_sf_ellint_RD(x, y, z, GSL_PREC_DOUBLE));
}

static double amplgsl_sf_ellint_RF(arglist *al) {
  double x = al->ra[0], y = al->ra[1], z = al->ra[2];
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al, gsl_sf_ellint_RF(x, y, z, GSL_PREC_DOUBLE));
}

static double amplgsl_sf_ellint_RJ(arglist *al) {
  double x = al->ra[0], y = al->ra[1], z = al->ra[2], p = al->ra[3];
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al, gsl_sf_ellint_RJ(x, y, z, p, GSL_PREC_DOUBLE));
}

static double amplgsl_sf_erf(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = 2 * exp(-x * x) / sqrt(M_PI);
    if (al->hes)
      *al->hes = -2 * x * *al->derivs;
  }
  return check_result(al, gsl_sf_erf(x));
}

static double amplgsl_sf_erfc(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = -2 * exp(-x * x) / sqrt(M_PI);
    if (al->hes)
      *al->hes = -2 * x * *al->derivs;
  }
  return check_result(al, gsl_sf_erfc(x));
}

static double amplgsl_sf_log_erfc(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    double erfc = gsl_sf_erfc(x);
    *al->derivs = -2 * exp(-x * x) / (sqrt(M_PI) * erfc);
    if (al->hes) {
      *al->hes = -2 * x * *al->derivs -
          ((4 * exp(-2 * x * x)) / (M_PI * erfc * erfc));
    }
  }
  return check_result(al, gsl_sf_log_erfc(x));
}

static double amplgsl_sf_erf_Z(arglist *al) {
  double x = al->ra[0];
  double z = gsl_sf_erf_Z(x);
  if (al->derivs) {
    *al->derivs = -x * z;
    if (al->hes)
      *al->hes = -(z + x * *al->derivs);
  }
  return check_result(al, z);
}

static double amplgsl_sf_erf_Q(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    double deriv = *al->derivs = -gsl_sf_erf_Z(x);
    if (al->hes)
      *al->hes = -x * deriv;
  }
  return check_result(al, gsl_sf_erf_Q(x));
}

static double amplgsl_sf_hazard(arglist *al) {
  double x = al->ra[0];
  double hazard = gsl_sf_hazard(x);
  if (al->derivs) {
    *al->derivs = (hazard - x) * hazard;
    if (al->hes)
      *al->hes = hazard * (hazard * (2 * hazard - 3 * x) + (x * x - 1));
  }
  return check_result(al, hazard);
}

static double amplgsl_sf_expint_E1(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = -exp(-x) / x;
    if (al->hes)
      *al->hes = -*al->derivs * (1 / x + 1);
  }
  return check_result(al, gsl_sf_expint_E1(x));
}

static double amplgsl_sf_expint_E2(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = -gsl_sf_expint_E1(x);
    if (al->hes)
      *al->hes = exp(-x) / x;
  }
  return check_result(al, gsl_sf_expint_E2(x));
}

static double amplgsl_sf_expint_En(arglist *al) {
  int n = (int)al->ra[0];
  double x = al->ra[1];
  if (!check_int_arg(al, 0, "n"))
    return 0;
  if (al->derivs) {
    al->derivs[1] = n != 0 ?
        -gsl_sf_expint_En(n - 1, x) : -exp(-x) * (1 / x + 1) / x;
    if (al->hes) {
      if (n == 0)
        al->hes[2] = exp(-x) * (1 + 2 * (1 + 1 / x) / x) / x;
      else if (n == 1)
        al->hes[2] = exp(-x) * (1 / x + 1) / x;
      else
        al->hes[2] = gsl_sf_expint_En(n - 2, x);
    }
  }
  return check_result(al, gsl_sf_expint_En(n, x));
}

static double amplgsl_sf_expint_Ei(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = exp(x) / x;
    if (al->hes)
      *al->hes = *al->derivs * (1 - 1 / x);
  }
  return check_result(al, gsl_sf_expint_Ei(x));
}

static double amplgsl_sf_Shi(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = x != 0 ? sinh(x) / x : 1;
    if (al->hes)
      *al->hes = x != 0 ? (cosh(x) - *al->derivs) / x : 0;
  }
  return check_result(al, gsl_sf_Shi(x));
}

static double amplgsl_sf_Chi(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = cosh(x) / x;
    if (al->hes)
      *al->hes = (sinh(x) - *al->derivs) / x;
  }
  return check_result(al, gsl_sf_Chi(x));
}

static double amplgsl_sf_expint_3(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = exp(-gsl_pow_3(x));
    if (al->hes)
      *al->hes = -3 * x * x * *al->derivs;
  }
  return check_result(al, gsl_sf_expint_3(x));
}

static double amplgsl_sf_Si(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = x != 0 ? sin(x) / x : 1;
    if (al->hes)
      *al->hes = x != 0 ? (cos(x) - *al->derivs) / x : 0;
  }
  return check_result(al, gsl_sf_Si(x));
}

static double amplgsl_sf_Ci(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = cos(x) / x;
    if (al->hes)
      *al->hes = -(sin(x) + *al->derivs) / x;
  }
  return check_result(al, gsl_sf_Ci(x));
}

static double amplgsl_sf_atanint(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = x != 0 ? atan(x) / x : 1;
    if (al->hes)
      *al->hes = x != 0 ? (1 / (x * x + 1) - *al->derivs) / x : 0;
  }
  return check_result(al, gsl_sf_atanint(x));
}

static double amplgsl_sf_fermi_dirac_m1(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = exp(x) / gsl_pow_2(exp(x) + 1);
    if (al->hes)
      *al->hes = -(exp(x) * (exp(x) - 1)) / gsl_pow_3(exp(x) + 1);
  }
  return check_result(al, gsl_sf_fermi_dirac_m1(x));
}

static double amplgsl_sf_fermi_dirac_0(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = gsl_sf_fermi_dirac_m1(x);
    if (al->hes)
      *al->hes = exp(x) / gsl_pow_2(exp(x) + 1);
  }
  return check_result(al, gsl_sf_fermi_dirac_0(x));
}

static double amplgsl_sf_fermi_dirac_1(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = gsl_sf_fermi_dirac_0(x);
    if (al->hes)
      *al->hes = gsl_sf_fermi_dirac_m1(x);
  }
  return check_result(al, gsl_sf_fermi_dirac_1(x));
}

static double amplgsl_sf_fermi_dirac_2(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = gsl_sf_fermi_dirac_1(x);
    if (al->hes)
      *al->hes = gsl_sf_fermi_dirac_0(x);
  }
  return check_result(al, gsl_sf_fermi_dirac_2(x));
}

static double amplgsl_sf_fermi_dirac_int(arglist *al) {
  int j = al->ra[0];
  double x = al->ra[1];
  if (!check_int_arg(al, 0, "j"))
    return 0;
  if (al->derivs) {
    al->derivs[1] = gsl_sf_fermi_dirac_int(j - 1, x);
    if (al->hes)
      al->hes[2] = gsl_sf_fermi_dirac_int(j - 2, x);
  }
  return check_result(al, gsl_sf_fermi_dirac_int(j, x));
}

static double amplgsl_sf_fermi_dirac_mhalf(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al, gsl_sf_fermi_dirac_mhalf(x));
}

static double amplgsl_sf_fermi_dirac_half(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al, gsl_sf_fermi_dirac_half(x));
}

static double amplgsl_sf_fermi_dirac_3half(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = gsl_sf_fermi_dirac_half(x);
    if (al->hes)
      *al->hes = gsl_sf_fermi_dirac_mhalf(x);
  }
  return check_result(al, gsl_sf_fermi_dirac_3half(x));
}

static double amplgsl_sf_fermi_dirac_inc_0(arglist *al) {
  double x = al->ra[0], b = al->ra[1];
  if (al->derivs) {
    double exp_x = exp(x), exp_b = exp(b);
    al->derivs[0] = exp_x / (exp_b + exp_x);
    al->derivs[1] = -al->derivs[0];
    if (al->hes) {
      al->hes[0] = al->hes[2] = al->derivs[0] * exp_b / (exp_b + exp_x);
      al->hes[1] = -al->hes[0];
    }
  }
  return check_result(al, gsl_sf_fermi_dirac_inc_0(x, b));
}

static double amplgsl_sf_gamma(arglist *al) {
  double x = al->ra[0];
  double gamma = gsl_sf_gamma(x);
  if (al->derivs) {
    double psi0 = gsl_sf_psi(x);
    *al->derivs = gamma * psi0;
    if (al->hes)
      *al->hes = *al->derivs * psi0 + gamma * gsl_sf_psi_1(x);
  }
  return check_result(al, gamma);
}

static double amplgsl_sf_lngamma(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = x >= 0 || round(x) != x ? gsl_sf_psi(x) : GSL_NAN;
    if (al->hes)
      *al->hes = gsl_sf_psi_1(x);
  }
  return check_result(al, gsl_sf_lngamma(x));
}

static double amplgsl_sf_gammastar(arglist *al) {
  double x = al->ra[0];
  double gammastar = gsl_sf_gammastar(x);
  if (al->derivs) {
    double coef = (0.5 / x - log(x) + gsl_sf_psi(x));
    *al->derivs = coef * gammastar;
    if (al->hes) {
      *al->hes = coef * *al->derivs +
          (gsl_sf_psi_1(x) - (1 + 0.5 / x) / x) * gammastar;
    }
  }
  return check_result(al, gammastar);
}

static double amplgsl_sf_gammainv(arglist *al) {
  double x = al->ra[0];
  double gammainv = gsl_sf_gammainv(x);
  if (al->derivs) {
    if (x > 0 || round(x) != x) {
      double psi0 = gsl_sf_psi(x);
      *al->derivs = -gammainv * psi0;
      if (al->hes)
        *al->hes = -*al->derivs * psi0 - gammainv * gsl_sf_psi_1(x);
    } else {
      *al->derivs = pow(-1, -x) * gsl_sf_gamma(1 - x);
      if (al->hes)
        *al->hes = -2 * *al->derivs * gsl_sf_psi(1 - x);
    }
  }
  return check_result(al, gammainv);
}

static double amplgsl_sf_poch(arglist *al) {
  double a = al->ra[0], x = al->ra[1];
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al, gsl_sf_poch(a, x));
}

static double amplgsl_sf_lnpoch(arglist *al) {
  double a = al->ra[0], x = al->ra[1];
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al, gsl_sf_lnpoch(a, x));
}

static double amplgsl_sf_pochrel(arglist *al) {
  double a = al->ra[0], x = al->ra[1];
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al, gsl_sf_pochrel(a, x));
}

static double amplgsl_sf_gamma_inc(arglist *al) {
  double a = al->ra[0], x = al->ra[1];
  if (al->derivs) {
    if (!check_const_arg(al, 0, "a"))
      return 0;
    al->derivs[1] = x != 0 ? -exp(-x) * pow(x, a - 1) : GSL_NAN;
    if (al->hes)
      al->hes[2] = al->derivs[1] * (a - x - 1) / x;
  }
  return check_result(al, gsl_sf_gamma_inc(a, x));
}

static double amplgsl_sf_gamma_inc_Q(arglist *al) {
  double a = al->ra[0], x = al->ra[1];
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al, gsl_sf_gamma_inc_Q(a, x));
}

static double amplgsl_sf_gamma_inc_P(arglist *al) {
  double a = al->ra[0], x = al->ra[1];
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al, gsl_sf_gamma_inc_P(a, x));
}

static double amplgsl_sf_beta(arglist *al) {
  double a = al->ra[0], b = al->ra[1];
  double beta = gsl_sf_beta(a, b);
  if (al->derivs) {
    double psi_a_plus_b = gsl_sf_psi(a + b);
    double da_coef = 0, db_coef = 0;
    int need_da = 1, need_db = 1;
    if (al->dig) {
      need_da = !al->dig[0];
      need_db = !al->dig[1];
    }
    if (need_da) {
      da_coef = gsl_sf_psi(a) - psi_a_plus_b;
      al->derivs[0] = beta * da_coef;
    }
    if (need_db) {
      db_coef = gsl_sf_psi(b) - psi_a_plus_b;
      al->derivs[1] = beta * db_coef;
    }
    if (al->hes) {
      double psi1_a_plus_b = gsl_sf_psi_1(a + b);
      if (need_da) {
        al->hes[0] = al->derivs[0] * da_coef +
            beta * (gsl_sf_psi_1(a) - psi1_a_plus_b);
        if (need_db)
          al->hes[1] = al->derivs[0] * db_coef - beta * psi1_a_plus_b;
      }
      if (need_db) {
       al->hes[2] = al->derivs[1] * db_coef +
           beta * (gsl_sf_psi_1(b) - psi1_a_plus_b);
      }
    }
  }
  return check_result(al, beta);
}

static double amplgsl_sf_lnbeta(arglist *al) {
  double a = al->ra[0], b = al->ra[1];
  if (al->derivs) {
    double psi_a_plus_b = gsl_sf_psi(a + b);
    int need_da = 1, need_db = 1;
    if (al->dig) {
      need_da = !al->dig[0];
      need_db = !al->dig[1];
    }
    if (need_da)
      al->derivs[0] = gsl_sf_psi(a) - psi_a_plus_b;
    if (need_db)
      al->derivs[1] = gsl_sf_psi(b) - psi_a_plus_b;
    if (al->hes) {
      double psi1_a_plus_b = gsl_sf_psi_1(a + b);
      if (need_da) {
        al->hes[0] = gsl_sf_psi_1(a) - psi1_a_plus_b;
        if (need_db)
          al->hes[1] = -psi1_a_plus_b;
      }
      if (need_db)
       al->hes[2] = gsl_sf_psi_1(b) - psi1_a_plus_b;
    }
  }
  return check_result(al, gsl_sf_lnbeta(a, b));
}

static double amplgsl_sf_beta_inc(arglist *al) {
  double a = al->ra[0], b = al->ra[1], x = al->ra[2];
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al, gsl_sf_beta_inc(a, b, x));
}

static double amplgsl_sf_gegenpoly_1(arglist *al) {
  double lambda = al->ra[0], x = al->ra[1];
  if (al->derivs) {
    al->derivs[0] = 2 * x;
    /* For unclear reason gsl_sf_gegenpoly_1(0, x) returns 2 * x. */
    al->derivs[1] = lambda != 0 ? 2 * lambda : 2;
    if (al->hes) {
      al->hes[0] = al->hes[2] = 0;
      al->hes[1] = 2;
    }
  }
  return check_result(al, gsl_sf_gegenpoly_1(lambda, x));
}

static double amplgsl_sf_gegenpoly_2(arglist *al) {
  double lambda = al->ra[0], x = al->ra[1];
  if (al->derivs) {
    double coef1 = (1 + 2 * lambda) * x;
    double coef2 = 4 * lambda * (lambda + 1);
    al->derivs[0] = 2 * coef1 * x - 1;
    /* For unclear reason gsl_sf_gegenpoly_2(0, x) returns 2 * x^2 - 1. */
    al->derivs[1] = lambda != 0 ? coef2 * x : 4 * x;
    if (al->hes) {
      al->hes[0] = 4 * x * x;
      al->hes[1] = 4 * coef1;
      al->hes[2] = lambda != 0 ? coef2 : 4;
    }
  }
  return check_result(al, gsl_sf_gegenpoly_2(lambda, x));
}

static double amplgsl_sf_gegenpoly_3(arglist *al) {
  double lambda = al->ra[0], x = al->ra[1];
  if (al->derivs) {
    double x2 = x * x;
    al->derivs[0] =
        x * (4 * (2.0 / 3.0 + lambda * (lambda + 2)) * x2 -
        2 * (2 * lambda + 1));
    /* For unclear reason gsl_sf_gegenpoly_3(0, x) returns
       x * (-2.0 + 4.0 / 3.0 * x * x). */
    al->derivs[1] = lambda != 0 ?
        2 * lambda * (lambda + 1) * (2 * (2 + lambda) * x2 - 1) :
        4 * x2 - 2;
    if (al->hes) {
      al->hes[0] = 4 * x * (2 * (lambda + 1) * x2 - 1);
      al->hes[1] = 2 * (x2 * (6 * lambda * lambda + 4) +
          2 * lambda * (6 * x2 - 1) - 1);
      al->hes[2] = lambda != 0 ?
          8 * lambda * (lambda + 1) * (lambda + 2) * x : 8 * x;
    }
  }
  return check_result(al, gsl_sf_gegenpoly_3(lambda, x));
}

static double amplgsl_sf_gegenpoly_n(arglist *al) {
  int n = (int)al->ra[0];
  double lambda = al->ra[1], x = al->ra[2];
  if (!check_int_arg(al, 0, "n"))
    return 0;
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al, gsl_sf_gegenpoly_n(n, lambda, x));
}

static double amplgsl_sf_hyperg_0F1(arglist *al) {
  double c = 0, x = 0;
  if (!check_args(al))
    return 0;
  c = al->ra[0];
  x = al->ra[1];
  if (al->derivs) {
    if (!check_const_arg(al, 0, "c"))
      return 0;
    al->derivs[1] = gsl_sf_hyperg_0F1(c + 1, x) / c;
    if (al->hes)
      al->hes[2] = gsl_sf_hyperg_0F1(c + 2, x) / (c * (c + 1));
  }
  return check_result(al, gsl_sf_hyperg_0F1(c, x));
}

static double amplgsl_sf_hyperg_1F1_int(arglist *al) {
  int m = 0, n = 0;
  double x = 0;
  if (!check_args(al) || !check_int_arg(al, 0, "m") ||
      !check_int_arg(al, 1, "n")) {
    return 0;
  }
  m = (int)al->ra[0];
  n = (int)al->ra[1];
  x = al->ra[2];
  if (al->derivs) {
    /* If n is an integer <= 0, then 1F1(m; n; x) is undefined.
       See http://mathworld.wolfram.com/
       ConfluentHypergeometricFunctionoftheFirstKind.html */
    al->derivs[2] = n > 0 ?
        m * gsl_sf_hyperg_1F1_int(m + 1, n + 1, x) / n : GSL_NAN;
    if (al->hes) {
      al->hes[5] =
          m * (m + 1) * gsl_sf_hyperg_1F1_int(m + 2, n + 2, x) / (n * (n + 1));
    }
  }
  return check_result(al, gsl_sf_hyperg_1F1_int(m, n, x));
}

static double amplgsl_sf_hyperg_1F1(arglist *al) {
  if (!check_args(al))
    return 0;
  /* gsl_sf_hyperg_1F1 seems to be broken, in particular,
     gsl_sf_hyperg_1F1(-2, -0.23, -5) returns 183.641 instead of -183.641. */
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al, gsl_sf_hyperg_1F1(al->ra[0], al->ra[1], al->ra[2]));
}

static double amplgsl_sf_hyperg_U_int(arglist *al) {
  if (!check_args(al) || !check_int_arg(al, 0, "m") ||
      !check_int_arg(al, 1, "n")) {
    return 0;
  }
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al,
      gsl_sf_hyperg_U_int((int)al->ra[0], (int)al->ra[1], al->ra[2]));
}

static double amplgsl_sf_hyperg_U(arglist *al) {
  if (!check_args(al))
    return 0;
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al, gsl_sf_hyperg_U(al->ra[0], al->ra[1], al->ra[2]));
}

static double amplgsl_sf_hyperg_2F1(arglist *al) {
  if (!check_args(al))
    return 0;
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al,
      gsl_sf_hyperg_2F1(al->ra[0], al->ra[1], al->ra[2], al->ra[3]));
}

static double amplgsl_sf_hyperg_2F1_conj(arglist *al) {
  if (!check_args(al))
    return 0;
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al,
      gsl_sf_hyperg_2F1_conj(al->ra[0], al->ra[1], al->ra[2], al->ra[3]));
}

static double amplgsl_sf_hyperg_2F1_renorm(arglist *al) {
  if (!check_args(al))
    return 0;
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al,
      gsl_sf_hyperg_2F1_renorm(al->ra[0], al->ra[1], al->ra[2], al->ra[3]));
}

static double amplgsl_sf_hyperg_2F0(arglist *al) {
  if (!check_args(al))
    return 0;
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al, gsl_sf_hyperg_2F0(al->ra[0], al->ra[1], al->ra[2]));
}

static double amplgsl_sf_laguerre_1(arglist *al) {
  double a = 0, x = 0;
  if (!check_args(al))
    return 0;
  a = al->ra[0];
  x = al->ra[1];
  if (al->derivs) {
    al->derivs[0] =  1;
    al->derivs[1] = -1;
    if (al->hes)
      al->hes[0] = al->hes[1] = al->hes[2] = 0;
  }
  return check_result(al, gsl_sf_laguerre_1(a, x));
}

static double amplgsl_sf_laguerre_2(arglist *al) {
  double a = 0, x = 0;
  if (!check_args(al))
    return 0;
  a = al->ra[0];
  x = al->ra[1];
  if (al->derivs) {
    al->derivs[0] =  a - x + 1.5;
    al->derivs[1] = -a + x - 2;
    if (al->hes) {
      al->hes[0] = al->hes[2] = 1;
      al->hes[1] = -1;
    }
  }
  return check_result(al, gsl_sf_laguerre_2(a, x));
}

static double amplgsl_sf_laguerre_3(arglist *al) {
  double a = 0, x = 0;
  if (!check_args(al))
    return 0;
  a = al->ra[0];
  x = al->ra[1];
  if (al->derivs) {
    al->derivs[0] = (11 + 3 * a * (a - 2 * (x - 2)) + 3 * x * (x - 5)) / 6;
    al->derivs[1] = (-0.5 * (a + 2) + x) * (a + 3) - 0.5 * x * x;
    if (al->hes) {
      al->hes[0] =  a - x + 2;
      al->hes[1] = -a + x - 2.5;
      al->hes[2] =  a - x + 3;
    }
  }
  return check_result(al, gsl_sf_laguerre_3(a, x));
}

static double amplgsl_sf_laguerre_n(arglist *al) {
  int n = 0;
  double a = 0, x = 0;
  if (!check_args(al) || !check_int_arg(al, 0, "n"))
    return 0;
  n = (int)al->ra[0];
  a = al->ra[1];
  x = al->ra[2];
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al, gsl_sf_laguerre_n(n, a, x));
}

static double amplgsl_sf_lambert_W0(arglist *al) {
  double x = al->ra[0];
  double value = gsl_sf_lambert_W0(x);
  if (al->derivs) {
    if (x < -1 / M_E)
      *al->derivs = GSL_NAN;
    else
      *al->derivs = x != 0 ? value / (x * (value + 1)) : 1;
    if (al->hes)
      *al->hes = -*al->derivs * *al->derivs * (value + 2) / (value + 1);
  }
  return check_result(al, value);
}

static double amplgsl_sf_lambert_Wm1(arglist *al) {
  double x = al->ra[0];
  double value = gsl_sf_lambert_Wm1(x);
  if (al->derivs) {
    if (x < -1 / M_E)
      *al->derivs = GSL_NAN;
    else
      *al->derivs = value / (x * (value + 1));
    if (al->hes)
      *al->hes = -*al->derivs * *al->derivs * (value + 2) / (value + 1);
  }
  return check_result(al, value);
}

static double amplgsl_sf_legendre_P1(arglist *al) {
  if (al->derivs) {
    *al->derivs = 1;
    if (al->hes)
      *al->hes = 0;
  }
  return check_result(al, gsl_sf_legendre_P1(al->ra[0]));
}

static double amplgsl_sf_legendre_P2(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = 3 * x;
    if (al->hes)
      *al->hes = 3;
  }
  return check_result(al, gsl_sf_legendre_P2(x));
}

static double amplgsl_sf_legendre_P3(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = 7.5 * x * x - 1.5;
    if (al->hes)
      *al->hes = 15 * x;
  }
  return check_result(al, gsl_sf_legendre_P3(x));
}

static double amplgsl_sf_legendre_Pl(arglist *al) {
  int el = 0;
  double x = 0, pl = 0;
  if (!check_int_arg(al, 0, "l"))
    return 0;
  el = al->ra[0];
  x = al->ra[1];
  pl = gsl_sf_legendre_Pl(el, x);
  if (al->derivs) {
    if (fabs(x) != 1) {
      double pl_plus_1 = gsl_sf_legendre_Pl(el + 1, x);
      double coef = (el + 1) / (x * x - 1);
      al->derivs[1] = -coef * (x * pl - pl_plus_1);
      if (al->hes) {
        al->hes[2] =
            coef * ((x * x * (el + 2) + 1) * pl - (2 * el + 5) * x * pl_plus_1 +
            (el + 2) * gsl_sf_legendre_Pl(el + 2, x)) / (x * x - 1);
      }
    } else {
      double coef = 0.5 * el * (el + 1);
      al->derivs[1] = pow(x, el + 1) * coef;
      if (al->hes)
        al->hes[2] = pow(x, el) * 0.25 * coef * (el * el + el - 2);
    }
  }
  return check_result(al, pl);
}

static double amplgsl_sf_legendre_Q0(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = 1 / (1 - x * x);
    if (al->hes)
      *al->hes = 2 * x * *al->derivs * *al->derivs;
  }
  return check_result(al, gsl_sf_legendre_Q0(x));
}

static double amplgsl_sf_legendre_Q1(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    double coef = 1 / (1 - x * x);
    *al->derivs = coef * x + 0.5 * (log(1 + x) - log(fabs(1 - x)));
    if (al->hes)
      *al->hes = 2 * coef * coef;
  }
  return check_result(al, gsl_sf_legendre_Q1(x));
}

static double amplgsl_sf_legendre_Ql(arglist *al) {
  int el = 0;
  double x = 0, ql = 0;
  if (!check_int_arg(al, 0, "l"))
    return 0;
  el = al->ra[0];
  x = al->ra[1];
  ql = gsl_sf_legendre_Ql(el, x);
  if (al->derivs) {
    double coef = (el + 1) / (x * x - 1);
    double ql_plus_1 = gsl_sf_legendre_Ql(el + 1, x);
    al->derivs[1] = coef * (ql_plus_1 - x * ql);
    if (al->hes) {
      al->hes[2] =
          coef * ((x * x * (el + 2) + 1) * ql - (2 * el + 5) * x * ql_plus_1 +
          (el + 2) * gsl_sf_legendre_Ql(el + 2, x)) / (x * x - 1);
    }
  }
  return check_result(al, ql);
}

#define ADDFUNC(name, num_args) \
    addfunc(#name, ampl##name, FUNCADD_REAL_VALUED, num_args, #name);

void funcadd_ASL(AmplExports *ae) {
  /* Don't call abort on error. */
  gsl_set_error_handler_off();

  /* Elementary Functions */
  ADDFUNC(gsl_log1p, 1);
  ADDFUNC(gsl_expm1, 1);
  ADDFUNC(gsl_hypot, 2);
  ADDFUNC(gsl_hypot3, 3);

  /* AMPL has built-in functions acosh, asinh and atanh so wrappers
     are not provided for their GSL equivalents. */

  /* Wrappers for functions operating on complex numbers are not provided
     since this requires support for structures/tuples as function arguments. */

  /* Airy Functions */
  ADDFUNC(gsl_sf_airy_Ai, 1);;
  ADDFUNC(gsl_sf_airy_Bi, 1);
  ADDFUNC(gsl_sf_airy_Ai_scaled, 1);
  ADDFUNC(gsl_sf_airy_Bi_scaled, 1);

  /* Zeros of Airy Functions */
  ADDFUNC(gsl_sf_airy_zero_Ai, 1);
  ADDFUNC(gsl_sf_airy_zero_Bi, 1);

  /* Zeros of Derivatives of Airy Functions */
  ADDFUNC(gsl_sf_airy_zero_Ai_deriv, 1);
  ADDFUNC(gsl_sf_airy_zero_Bi_deriv, 1);

  /* Bessel Functions */
  ADDFUNC(gsl_sf_bessel_J0, 1);
  ADDFUNC(gsl_sf_bessel_J1, 1);
  ADDFUNC(gsl_sf_bessel_Jn, 2);

  /* Irregular Cylindrical Bessel Functions */
  ADDFUNC(gsl_sf_bessel_Y0, 1);
  ADDFUNC(gsl_sf_bessel_Y1, 1);
  ADDFUNC(gsl_sf_bessel_Yn, 2);

  /* Regular Modified Cylindrical Bessel Functions */
  ADDFUNC(gsl_sf_bessel_I0, 1);
  ADDFUNC(gsl_sf_bessel_I1, 1);
  ADDFUNC(gsl_sf_bessel_In, 2);
  ADDFUNC(gsl_sf_bessel_I0_scaled, 1);
  ADDFUNC(gsl_sf_bessel_I1_scaled, 1);
  ADDFUNC(gsl_sf_bessel_In_scaled, 2);

  /* Irregular Modified Cylindrical Bessel Functions */
  ADDFUNC(gsl_sf_bessel_K0, 1);
  ADDFUNC(gsl_sf_bessel_K1, 1);
  ADDFUNC(gsl_sf_bessel_Kn, 2);
  ADDFUNC(gsl_sf_bessel_K0_scaled, 1);
  ADDFUNC(gsl_sf_bessel_K1_scaled, 1);
  ADDFUNC(gsl_sf_bessel_Kn_scaled, 2);

  /* Regular Spherical Bessel Functions */
  ADDFUNC(gsl_sf_bessel_j0, 1);
  ADDFUNC(gsl_sf_bessel_j1, 1);
  ADDFUNC(gsl_sf_bessel_j2, 1);
  ADDFUNC(gsl_sf_bessel_jl, 2);

  /* Irregular Spherical Bessel Functions */
  ADDFUNC(gsl_sf_bessel_y0, 1);
  ADDFUNC(gsl_sf_bessel_y1, 1);
  ADDFUNC(gsl_sf_bessel_y2, 1);
  ADDFUNC(gsl_sf_bessel_yl, 2);

  /* Regular Modified Spherical Bessel Functions */
  ADDFUNC(gsl_sf_bessel_i0_scaled, 1);
  ADDFUNC(gsl_sf_bessel_i1_scaled, 1);
  ADDFUNC(gsl_sf_bessel_i2_scaled, 1);
  ADDFUNC(gsl_sf_bessel_il_scaled, 2);

  /* Irregular Modified Spherical Bessel Functions */
  ADDFUNC(gsl_sf_bessel_k0_scaled, 1);
  ADDFUNC(gsl_sf_bessel_k1_scaled, 1);
  ADDFUNC(gsl_sf_bessel_k2_scaled, 1);
  ADDFUNC(gsl_sf_bessel_kl_scaled, 2);

  /* Regular Bessel Function - Fractional Order */
  ADDFUNC(gsl_sf_bessel_Jnu, 2);

  /* Irregular Bessel Functions - Fractional Order */
  ADDFUNC(gsl_sf_bessel_Ynu, 2);

  /* Regular Modified Bessel Functions - Fractional Order */
  ADDFUNC(gsl_sf_bessel_Inu, 2);
  ADDFUNC(gsl_sf_bessel_Inu_scaled, 2);

  /* Irregular Modified Bessel Functions - Fractional Order */
  ADDFUNC(gsl_sf_bessel_Knu, 2);
  ADDFUNC(gsl_sf_bessel_lnKnu, 2);
  ADDFUNC(gsl_sf_bessel_Knu_scaled, 2);

  /* Zeros of Regular Bessel Functions */
  ADDFUNC(gsl_sf_bessel_zero_J0, 1);
  ADDFUNC(gsl_sf_bessel_zero_J1, 1);
  ADDFUNC(gsl_sf_bessel_zero_Jnu, 2);

  /* Clausen Function */
  ADDFUNC(gsl_sf_clausen, 1);

  /* Normalized Hydrogenic Bound States */
  ADDFUNC(gsl_sf_hydrogenicR_1, 2);
  ADDFUNC(gsl_sf_hydrogenicR, 4);

  /* Coulomb Wave Function Normalization Constant */
  ADDFUNC(gsl_sf_coulomb_CL, 2);

  /* Coupling Coefficients */
  ADDFUNC(gsl_sf_coupling_3j, 6);
  ADDFUNC(gsl_sf_coupling_6j, 6);
  ADDFUNC(gsl_sf_coupling_9j, 9);

  /* Dawson Function */
  ADDFUNC(gsl_sf_dawson, 1);

  /* Debye Functions */
  ADDFUNC(gsl_sf_debye_1, 1);
  ADDFUNC(gsl_sf_debye_2, 1);
  ADDFUNC(gsl_sf_debye_3, 1);
  ADDFUNC(gsl_sf_debye_4, 1);
  ADDFUNC(gsl_sf_debye_5, 1);
  ADDFUNC(gsl_sf_debye_6, 1);

  /* Dilogarithm */
  ADDFUNC(gsl_sf_dilog, 1);

  /* Legendre Form of Complete Elliptic Integrals */
  ADDFUNC(gsl_sf_ellint_Kcomp, 1);
  ADDFUNC(gsl_sf_ellint_Ecomp, 1);
  ADDFUNC(gsl_sf_ellint_Pcomp, 2);

  /* Legendre Form of Incomplete Elliptic Integrals */
  ADDFUNC(gsl_sf_ellint_F, 2);
  ADDFUNC(gsl_sf_ellint_E, 2);
  ADDFUNC(gsl_sf_ellint_P, 3);
  ADDFUNC(gsl_sf_ellint_D, 3);

  /* Carlson Forms */
  ADDFUNC(gsl_sf_ellint_RC, 2);
  ADDFUNC(gsl_sf_ellint_RD, 3);
  ADDFUNC(gsl_sf_ellint_RF, 3);
  ADDFUNC(gsl_sf_ellint_RJ, 4);

  /* Elliptic Functions (Jacobi) */
  /* Wrapper for gsl_sf_elljac_e is not provided since the latter produces
     multiple values (through output parameters). */

  /* Error Functions */
  ADDFUNC(gsl_sf_erf, 1);
  ADDFUNC(gsl_sf_erfc, 1);
  ADDFUNC(gsl_sf_log_erfc, 1);
  ADDFUNC(gsl_sf_erf_Z, 1);
  ADDFUNC(gsl_sf_erf_Q, 1);
  ADDFUNC(gsl_sf_hazard, 1);

  /* Exponential Integrals */
  ADDFUNC(gsl_sf_expint_E1, 1);
  ADDFUNC(gsl_sf_expint_E2, 1);
  ADDFUNC(gsl_sf_expint_En, 2);
  ADDFUNC(gsl_sf_expint_Ei, 1);
  ADDFUNC(gsl_sf_Shi, 1);
  ADDFUNC(gsl_sf_Chi, 1);
  ADDFUNC(gsl_sf_expint_3, 1);
  ADDFUNC(gsl_sf_Si, 1);
  ADDFUNC(gsl_sf_Ci, 1);
  ADDFUNC(gsl_sf_atanint, 1);

  /* Fermi-Dirac Function */
  ADDFUNC(gsl_sf_fermi_dirac_m1, 1);
  ADDFUNC(gsl_sf_fermi_dirac_0, 1);
  ADDFUNC(gsl_sf_fermi_dirac_1, 1);
  ADDFUNC(gsl_sf_fermi_dirac_2, 1);
  ADDFUNC(gsl_sf_fermi_dirac_int, 2);
  ADDFUNC(gsl_sf_fermi_dirac_mhalf, 1);
  ADDFUNC(gsl_sf_fermi_dirac_half, 1);
  ADDFUNC(gsl_sf_fermi_dirac_3half, 1);
  ADDFUNC(gsl_sf_fermi_dirac_inc_0, 2);

  /* Gamma Functions */
  ADDFUNC(gsl_sf_gamma, 1);
  ADDFUNC(gsl_sf_lngamma, 1);
  ADDFUNC(gsl_sf_gammastar, 1);
  ADDFUNC(gsl_sf_gammainv, 1);

  /* Wrapper for factorials are not provided since these are easily
     implemented using built-in AMPL features like the prod operator. */

  /* Pochhammer Symbol */
  ADDFUNC(gsl_sf_poch, 2);
  ADDFUNC(gsl_sf_lnpoch, 2);
  ADDFUNC(gsl_sf_pochrel, 2);

  /* Incomplete Gamma Functions */
  ADDFUNC(gsl_sf_gamma_inc, 2);
  ADDFUNC(gsl_sf_gamma_inc_Q, 2);
  ADDFUNC(gsl_sf_gamma_inc_P, 2);

  /* Beta Functions */
  ADDFUNC(gsl_sf_beta, 2);
  ADDFUNC(gsl_sf_lnbeta, 2);

  /* Incomplete Beta Function */
  ADDFUNC(gsl_sf_beta_inc, 3);

  /* Gegenbauer Functions */
  ADDFUNC(gsl_sf_gegenpoly_1, 2);
  ADDFUNC(gsl_sf_gegenpoly_2, 2);
  ADDFUNC(gsl_sf_gegenpoly_3, 2);
  ADDFUNC(gsl_sf_gegenpoly_n, 3);

  /* Hypergeometric Functions */
  ADDFUNC(gsl_sf_hyperg_0F1, 2);
  ADDFUNC(gsl_sf_hyperg_1F1_int, 3);
  ADDFUNC(gsl_sf_hyperg_1F1, 3);
  ADDFUNC(gsl_sf_hyperg_U_int, 3);
  ADDFUNC(gsl_sf_hyperg_U, 3);
  ADDFUNC(gsl_sf_hyperg_2F1, 4);
  ADDFUNC(gsl_sf_hyperg_2F1_conj, 4);
  ADDFUNC(gsl_sf_hyperg_2F1_renorm, 4);
  ADDFUNC(gsl_sf_hyperg_2F0, 3);

  /* Laguerre Functions */
  ADDFUNC(gsl_sf_laguerre_1, 2);
  ADDFUNC(gsl_sf_laguerre_2, 2);
  ADDFUNC(gsl_sf_laguerre_3, 2);
  ADDFUNC(gsl_sf_laguerre_n, 3);

  /* Lambert W Functions */
  ADDFUNC(gsl_sf_lambert_W0, 1);
  ADDFUNC(gsl_sf_lambert_Wm1, 1);

  /* Legendre Polynomials */
  ADDFUNC(gsl_sf_legendre_P1, 1);
  ADDFUNC(gsl_sf_legendre_P2, 1);
  ADDFUNC(gsl_sf_legendre_P3, 1);
  ADDFUNC(gsl_sf_legendre_Pl, 2);
  ADDFUNC(gsl_sf_legendre_Q0, 1);
  ADDFUNC(gsl_sf_legendre_Q1, 1);
  ADDFUNC(gsl_sf_legendre_Ql, 2);
  // TODO
}
