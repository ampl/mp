/* AMPL bindings for the GNU Scientific Library. */

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
  if (al->derivs)
    *al->derivs = GSL_NAN;
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

/* Reports a function evaluation error printing the specified suffix
 * after the function name. */
static void eval_error_with_suffix(arglist *al,
    const char *func_name, const char *suffix) {
  int n = 0, i = 0;
  al->Errmsg = al->AE->Tempmem(al->TMI, MAX_ERROR_MESSAGE_SIZE);
  n += al->AE->SnprintF(al->Errmsg, MAX_ERROR_MESSAGE_SIZE,
      "can't evaluate %s%s(", func_name, suffix);
  for (i = 0; i < al->n - 1; ++i) {
    n += al->AE->SnprintF(al->Errmsg + n, MAX_ERROR_MESSAGE_SIZE - n,
        "%g, ", al->ra[i]);
  }
  al->AE->SnprintF(al->Errmsg + n, MAX_ERROR_MESSAGE_SIZE - n,
      "%g)", al->ra[al->n - 1]);
}

/* Reports a function evaluation error. */
static void eval_error(arglist *al, const char *func_name) {
  eval_error_with_suffix(al, func_name, "");
}

static double check_result(arglist *al, double result, const char *func_name) {
  int i = 0, n = 0;
  if (gsl_isnan(result)) {
    eval_error(al, func_name);
    return 0;
  }
  if (al->derivs) {
    for (i = 0; i < al->n; ++i) {
      if (gsl_isnan(al->derivs[i])) {
        eval_error_with_suffix(al, func_name, "'");
        return 0;
      }
    }
    if (al->hes) {
      for (i = 0, n = al->n * (al->n + 1) / 2; i < n; ++i) {
        if (gsl_isnan(al->hes[i])) {
          eval_error_with_suffix(al, func_name, "''");
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
  return check_result(al, gsl_log1p(x), "gsl_log1p");
}

static double amplgsl_expm1(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    double deriv = *al->derivs = exp(x);
    if (al->hes)
      *al->hes = deriv;
  }
  return gsl_expm1(x);
}

static double amplgsl_hypot(arglist *al) {
  double x = al->ra[0];
  double y = al->ra[1];
  double hypot = gsl_hypot(x, y);
  if (al->derivs) {
    real *derivs = al->derivs;
    if (hypot == 0) {
      eval_error(al, "gsl_hypot'");
      return 0;
    }
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

static double amplgsl_hypot3(arglist *al) {
  double x = al->ra[0];
  double y = al->ra[1];
  double z = al->ra[2];
  double hypot = gsl_hypot3(x, y, z);
  if (al->derivs) {
    real *derivs = al->derivs;
    if (hypot == 0) {
      eval_error(al, "gsl_hypot3'");
      return 0;
    }
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
  return hypot;
}

static double amplgsl_sf_airy_Ai(arglist *al) {
  double x = al->ra[0];
  double value = gsl_sf_airy_Ai(x, GSL_PREC_DOUBLE);
  if (al->derivs) {
    *al->derivs = gsl_sf_airy_Ai_deriv(x, GSL_PREC_DOUBLE);
    if (al->hes)
      *al->hes = x * value;
  }
  return value;
}

static double amplgsl_sf_airy_Bi(arglist *al) {
  double x = al->ra[0];
  double value = gsl_sf_airy_Bi(x, GSL_PREC_DOUBLE);
  if (al->derivs) {
    *al->derivs = gsl_sf_airy_Bi_deriv(x, GSL_PREC_DOUBLE);
    if (al->hes)
      *al->hes = x * value;
  }
  return value;
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
  return check_result(al, value, "gsl_sf_airy_Ai_scaled");
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
  return check_result(al, value, "gsl_sf_airy_Bi_scaled");
}

static double amplgsl_sf_airy_zero_Ai(arglist *al) {
  if (!check_zero_func_args(al, 0))
    return 0;
  return check_result(al, gsl_sf_airy_zero_Ai(al->ra[0]),
      "gsl_sf_airy_zero_Ai");
}

static double amplgsl_sf_airy_zero_Bi(arglist *al) {
  if (!check_zero_func_args(al, 0))
    return 0;
  return check_result(al, gsl_sf_airy_zero_Bi(al->ra[0]),
      "gsl_sf_airy_zero_Bi");
}

static double amplgsl_sf_airy_zero_Ai_deriv(arglist *al) {
  if (!check_zero_func_args(al, 0))
    return 0;
  return check_result(al, gsl_sf_airy_zero_Ai_deriv(al->ra[0]),
      "gsl_sf_airy_zero_Ai_deriv");
}

static double amplgsl_sf_airy_zero_Bi_deriv(arglist *al) {
  if (!check_zero_func_args(al, 0))
    return 0;
  return check_result(al, gsl_sf_airy_zero_Bi_deriv(al->ra[0]),
      "gsl_sf_airy_zero_Bi_deriv");
}

static double amplgsl_sf_bessel_J0(arglist *al) {
  double x = al->ra[0];
  double j0 = gsl_sf_bessel_J0(x);
  if (al->derivs) {
    *al->derivs = -gsl_sf_bessel_J1(x);
    if (al->hes)
      *al->hes = 0.5 * (gsl_sf_bessel_Jn(2, x) - j0);
  }
  return j0;
}

static double amplgsl_sf_bessel_J1(arglist *al) {
  double x = al->ra[0];
  double j1 = gsl_sf_bessel_J1(x);
  if (al->derivs) {
    *al->derivs = 0.5 * (gsl_sf_bessel_J0(x) - gsl_sf_bessel_Jn(2, x));
    if (al->hes)
      *al->hes = 0.25 * (gsl_sf_bessel_Jn(3, x) - 3 * j1);
  }
  return j1;
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
  return jn;
}

static double amplgsl_sf_bessel_Y0(arglist *al) {
  double x = al->ra[0];
  double y0 = gsl_sf_bessel_Y0(x);
  if (al->derivs) {
    *al->derivs = -gsl_sf_bessel_Y1(x);
    if (al->hes)
      *al->hes = 0.5 * (gsl_sf_bessel_Yn(2, x) - y0);
  }
  return check_result(al, y0, "gsl_sf_bessel_Y0");
}

static double amplgsl_sf_bessel_Y1(arglist *al) {
  double x = al->ra[0];
  double y1 = gsl_sf_bessel_Y1(x);
  if (al->derivs) {
    *al->derivs = 0.5 * (gsl_sf_bessel_Y0(x) - gsl_sf_bessel_Yn(2, x));
    if (al->hes)
      *al->hes = 0.25 * (gsl_sf_bessel_Yn(3, x) - 3 * y1);
  }
  return check_result(al, y1, "gsl_sf_bessel_Y1");
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
  return check_result(al, yn, "gsl_sf_bessel_Yn");
}

static double amplgsl_sf_bessel_I0(arglist *al) {
  double x = al->ra[0];
  double i0 = gsl_sf_bessel_I0(x);
  if (al->derivs) {
    *al->derivs = gsl_sf_bessel_I1(x);
    if (al->hes)
      *al->hes = 0.5 * (gsl_sf_bessel_In(2, x) + i0);
  }
  return i0;
}

static double amplgsl_sf_bessel_I1(arglist *al) {
  double x = al->ra[0];
  double i1 = gsl_sf_bessel_I1(x);
  if (al->derivs) {
    *al->derivs = 0.5 * (gsl_sf_bessel_I0(x) + gsl_sf_bessel_In(2, x));
    if (al->hes)
      *al->hes = 0.25 * (gsl_sf_bessel_In(3, x) + 3 * i1);
  }
  return i1;
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
  return in;
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
  return check_result(al, i0, "gsl_sf_bessel_I0_scaled");
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
  return check_result(al, i1, "gsl_sf_bessel_I1_scaled");
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
  return check_result(al, in, "gsl_sf_bessel_In_scaled");
}

static double amplgsl_sf_bessel_K0(arglist *al) {
  double x = al->ra[0];
  double k0 = gsl_sf_bessel_K0(x);
  if (al->derivs) {
    *al->derivs = -gsl_sf_bessel_K1(x);
    if (al->hes)
      *al->hes = 0.5 * (gsl_sf_bessel_Kn(2, x) + k0);
  }
  return check_result(al, k0, "gsl_sf_bessel_K0");
}

static double amplgsl_sf_bessel_K1(arglist *al) {
  double x = al->ra[0];
  double k1 = gsl_sf_bessel_K1(x);
  if (al->derivs) {
    *al->derivs = -0.5 * (gsl_sf_bessel_K0(x) + gsl_sf_bessel_Kn(2, x));
    if (al->hes)
      *al->hes = 0.25 * (gsl_sf_bessel_Kn(3, x) + 3 * k1);
  }
  return check_result(al, k1, "gsl_sf_bessel_K1");
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
  return check_result(al, kn, "gsl_sf_bessel_Kn");
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
  return check_result(al, k0, "gsl_sf_bessel_K0_scaled");
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
  return check_result(al, k1, "gsl_sf_bessel_K1_scaled");
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
  return check_result(al, kn, "gsl_sf_bessel_Kn_scaled");
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
  return check_result(al, gsl_sf_bessel_j0(x), "gsl_sf_bessel_j0");
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
  return check_result(al, j1, "gsl_sf_bessel_j1");
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
  return check_result(al, j2, "gsl_sf_bessel_j2");
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
  return check_result(al, jl, "gsl_sf_bessel_jl");
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
  return check_result(al, gsl_sf_bessel_y0(x), "gsl_sf_bessel_y0");
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
  return check_result(al, y1, "gsl_sf_bessel_y1");
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
  return check_result(al, y2, "gsl_sf_bessel_y2");
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
  return check_result(al, yl, "gsl_sf_bessel_yl");
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
    if (!check_result(al, *al->derivs, "gsl_sf_bessel_i0_scaled'"))
      return 0;
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
  return check_result(al, i0, "gsl_sf_bessel_i0_scaled");
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
  return check_result(al, i1, "gsl_sf_bessel_i1_scaled");
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
  return check_result(al, i2, "gsl_sf_bessel_i2_scaled");
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
    } else deriv = 0.5 * (il_minus_1 + coef * il + il_plus_1);
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
  return check_result(al, il, "gsl_sf_bessel_il_scaled");
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
  return check_result(al, k0, "gsl_sf_bessel_k0_scaled");
}

static double amplgsl_sf_bessel_k1_scaled(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    double pi_sqrt_inv_x = M_PI * sqrt(1 / x);
    *al->derivs = -(pi_sqrt_inv_x * (x + 2)) / (2 * pow(x, 2.5));
    if (al->hes)
      *al->hes = (pi_sqrt_inv_x * (x + 3)) / pow(x, 3.5);
  }
  return check_result(al, gsl_sf_bessel_k1_scaled(x),
      "gsl_sf_bessel_k1_scaled");
}

static double amplgsl_sf_bessel_k2_scaled(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    double pi_sqrt_inv_x = M_PI * sqrt(1 / x);
    *al->derivs = -pi_sqrt_inv_x * (x + 3) * (x + 3) / (2 * pow(x, 3.5));
    if (al->hes)
      *al->hes = pi_sqrt_inv_x * (x * x + 9 * x + 18) / pow(x, 4.5);
  }
  return check_result(al, gsl_sf_bessel_k2_scaled(x),
      "gsl_sf_bessel_k2_scaled");
}

static double amplgsl_sf_bessel_kl_scaled(arglist *al) {
  int el = al->ra[0];
  double x = al->ra[1];
  double kl = 0;
  if (!check_bessel_args(al, 0, "l"))
    return 0;
  kl = gsl_sf_bessel_kl_scaled(el, x);
  if (al->derivs) {
    double kl_minus_1 = el != 0 ? gsl_sf_bessel_kl_scaled(el - 1, x) : M_PI_2 / x;
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
  return check_result(al, kl, "gsl_sf_bessel_kl_scaled");
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
  return check_result(al, jn, "gsl_sf_bessel_Jnu");
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
  return check_result(al, yn, "gsl_sf_bessel_Ynu");
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
  return check_result(al, in, "gsl_sf_bessel_Inu");
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
  return check_result(al, in, "gsl_sf_bessel_Inu_scaled");
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
  return check_result(al, kn, "gsl_sf_bessel_Knu");
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
  return check_result(al, gsl_sf_bessel_lnKnu(n, x), "gsl_sf_bessel_lnKnu");
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
  return check_result(al, kn, "gsl_sf_bessel_Knu_scaled");
}

static double amplgsl_sf_bessel_zero_J0(arglist *al) {
  if (!check_zero_func_args(al, 0))
    return 0;
  return check_result(al, gsl_sf_bessel_zero_J0(al->ra[0]),
      "gsl_sf_bessel_zero_J0");
}

static double amplgsl_sf_bessel_zero_J1(arglist *al) {
  if (!check_zero_func_args(al, 0))
    return 0;
  return check_result(al, gsl_sf_bessel_zero_J1(al->ra[0]),
      "gsl_sf_bessel_zero_J1");
}

static double amplgsl_sf_bessel_zero_Jnu(arglist *al) {
  if (!check_zero_func_args(al, 1))
    return 0;
  return check_result(al, gsl_sf_bessel_zero_Jnu(al->ra[0], al->ra[1]),
      "gsl_sf_bessel_zero_Jnu");
}

static double amplgsl_sf_clausen(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = -log(2 * sin(0.5 * fabs(x)));
    if (al->hes)
      *al->hes = fmod(x, M_PI) != 0 ? -0.5 * tan(M_PI_2 - 0.5 * x) : GSL_NAN;
  }
  return check_result(al, gsl_sf_clausen(x), "gsl_sf_clausen");
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
  return check_result(al, gsl_sf_hydrogenicR_1(Z, r), "gsl_sf_hydrogenicR_1");
}

static double amplgsl_sf_hydrogenicR(arglist *al) {
  if (!check_int_arg(al, 0, "n") || !check_int_arg(al, 1, "l"))
    return 0;
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al, gsl_sf_hydrogenicR(
      al->ra[0], al->ra[1], al->ra[2], al->ra[3]), "gsl_sf_hydrogenicR");
}

static double amplgsl_sf_coulomb_CL(arglist *al) {
  gsl_sf_result result = {0};
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  if (gsl_sf_coulomb_CL_e(al->ra[0], al->ra[1], &result)) {
    eval_error(al, "gsl_sf_coulomb_CL");
    return 0;
  }
  return result.val;
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
      al->ra[0], al->ra[1], al->ra[2], al->ra[3], al->ra[4], al->ra[5]),
      "gsl_sf_coupling_3j");
}

static double amplgsl_sf_coupling_6j(arglist *al) {
  static const char *ARG_NAMES[] = {
      "two_ja", "two_jb", "two_jc",
      "two_jd", "two_je", "two_jf"
  };
  if (!check_coupling_args(al, ARG_NAMES))
    return 0;
  return check_result(al, gsl_sf_coupling_6j(
      al->ra[0], al->ra[1], al->ra[2], al->ra[3], al->ra[4], al->ra[5]),
      "gsl_sf_coupling_6j");
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
      al->ra[6], al->ra[7], al->ra[8]), "gsl_sf_coupling_9j");
}

static double amplgsl_sf_dawson(arglist *al) {
  double x = al->ra[0];
  double f = gsl_sf_dawson(x);
  if (al->derivs) {
    double deriv = *al->derivs = 1 - 2 * x * f;
    if (al->hes)
      *al->hes = - 2 * (f + x * deriv);
  }
  return f;
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

static double debye(arglist *al, int n,
    double (*func)(double), const char *name) {
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
  return check_result(al, f, name);
}

#define DEBYE(n) \
  static double amplgsl_sf_debye_##n(arglist *al) { \
    return debye(al, n, gsl_sf_debye_##n, "gsl_sf_debye_" #n); \
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
  return check_result(al, gsl_sf_dilog(x), "gsl_sf_dilog");
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
  return check_result(al, kcomp, "gsl_sf_ellint_Kcomp");
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
  return check_result(al, ecomp, "gsl_sf_ellint_Ecomp");
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
  return check_result(al, pcomp, "gsl_sf_ellint_Pcomp");
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
  return check_result(al, f, "gsl_sf_ellint_F");
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
  return check_result(al, e, "gsl_sf_ellint_E");
}

static double amplgsl_sf_ellint_P(arglist *al) {
  double phi = al->ra[0], k = al->ra[1], n = al->ra[2];
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al, gsl_sf_ellint_P(phi, k, n, GSL_PREC_DOUBLE),
      "gsl_sf_ellint_P");
}

static double amplgsl_sf_ellint_D(arglist *al) {
  double phi = al->ra[0], k = al->ra[1], n = al->ra[2];
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al, gsl_sf_ellint_D(phi, k, n, GSL_PREC_DOUBLE),
      "gsl_sf_ellint_D");
}

static double amplgsl_sf_ellint_RC(arglist *al) {
  double x = al->ra[0], y = al->ra[1];
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al, gsl_sf_ellint_RC(x, y, GSL_PREC_DOUBLE),
      "gsl_sf_ellint_RC");
}

static double amplgsl_sf_ellint_RD(arglist *al) {
  double x = al->ra[0], y = al->ra[1], z = al->ra[2];
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al, gsl_sf_ellint_RD(x, y, z, GSL_PREC_DOUBLE),
      "gsl_sf_ellint_RD");
}

static double amplgsl_sf_ellint_RF(arglist *al) {
  double x = al->ra[0], y = al->ra[1], z = al->ra[2];
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al, gsl_sf_ellint_RF(x, y, z, GSL_PREC_DOUBLE),
      "gsl_sf_ellint_RF");
}

static double amplgsl_sf_ellint_RJ(arglist *al) {
  double x = al->ra[0], y = al->ra[1], z = al->ra[2], p = al->ra[3];
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al, gsl_sf_ellint_RJ(x, y, z, p, GSL_PREC_DOUBLE),
      "gsl_sf_ellint_RJ");
}

static double amplgsl_sf_erf(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = 2 * exp(-x * x) / sqrt(M_PI);
    if (al->hes)
      *al->hes = -2 * x * *al->derivs;
  }
  return check_result(al, gsl_sf_erf(x), "gsl_sf_erf");
}

static double amplgsl_sf_erfc(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = -2 * exp(-x * x) / sqrt(M_PI);
    if (al->hes)
      *al->hes = -2 * x * *al->derivs;
  }
  return check_result(al, gsl_sf_erfc(x), "gsl_sf_erfc");
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
  return check_result(al, gsl_sf_log_erfc(x), "gsl_sf_log_erfc");
}

static double amplgsl_sf_erf_Z(arglist *al) {
  double x = al->ra[0];
  double z = gsl_sf_erf_Z(x);
  if (al->derivs) {
    *al->derivs = -x * z;
    if (al->hes)
      *al->hes = -(z + x * *al->derivs);
  }
  return check_result(al, z, "gsl_sf_erf_Z");
}

static double amplgsl_sf_erf_Q(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    double deriv = *al->derivs = -gsl_sf_erf_Z(x);
    if (al->hes)
      *al->hes = -x * deriv;
  }
  return check_result(al, gsl_sf_erf_Q(x), "gsl_sf_erf_Q");
}

static double amplgsl_sf_hazard(arglist *al) {
  double x = al->ra[0];
  double hazard = gsl_sf_hazard(x);
  if (al->derivs) {
    *al->derivs = (hazard - x) * hazard;
    if (al->hes)
      *al->hes = hazard * (hazard * (2 * hazard - 3 * x) + (x * x - 1));
  }
  return check_result(al, hazard, "gsl_sf_hazard");
}

static double amplgsl_sf_expint_E1(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = -exp(-x) / x;
    if (al->hes)
      *al->hes = -*al->derivs * (1 / x + 1);
  }
  return check_result(al, gsl_sf_expint_E1(x), "gsl_sf_expint_E1");
}

static double amplgsl_sf_expint_E2(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = -gsl_sf_expint_E1(x);
    if (al->hes)
      *al->hes = exp(-x) / x;
  }
  return check_result(al, gsl_sf_expint_E2(x), "gsl_sf_expint_E2");
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
  return check_result(al, gsl_sf_expint_En(n, x), "gsl_sf_expint_En");
}

static double amplgsl_sf_expint_Ei(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = exp(x) / x;
    if (al->hes)
      *al->hes = *al->derivs * (1 - 1 / x);
  }
  return check_result(al, gsl_sf_expint_Ei(x), "gsl_sf_expint_Ei");
}

static double amplgsl_sf_Shi(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = x != 0 ? sinh(x) / x : 1;
    if (al->hes)
      *al->hes = x != 0 ? (cosh(x) - *al->derivs) / x : 0;
  }
  return check_result(al, gsl_sf_Shi(x), "gsl_sf_Shi");
}

static double amplgsl_sf_Chi(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = cosh(x) / x;
    if (al->hes)
      *al->hes = (sinh(x) - *al->derivs) / x;
  }
  return check_result(al, gsl_sf_Chi(x), "gsl_sf_Chi");
}

static double amplgsl_sf_expint_3(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = exp(-gsl_pow_3(x));
    if (al->hes)
      *al->hes = -3 * x * x * *al->derivs;
  }
  return check_result(al, gsl_sf_expint_3(x), "gsl_sf_expint_3");
}

static double amplgsl_sf_Si(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = x != 0 ? sin(x) / x : 1;
    if (al->hes)
      *al->hes = x != 0 ? (cos(x) - *al->derivs) / x : 0;
  }
  return check_result(al, gsl_sf_Si(x), "gsl_sf_Si");
}

static double amplgsl_sf_Ci(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = cos(x) / x;
    if (al->hes)
      *al->hes = -(sin(x) + *al->derivs) / x;
  }
  return check_result(al, gsl_sf_Ci(x), "gsl_sf_Ci");
}

static double amplgsl_sf_atanint(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = x != 0 ? atan(x) / x : 1;
    if (al->hes)
      *al->hes = x != 0 ? (1 / (x * x + 1) - *al->derivs) / x : 0;
  }
  return check_result(al, gsl_sf_atanint(x), "gsl_sf_atanint");
}

static double amplgsl_sf_fermi_dirac_m1(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = exp(x) / gsl_pow_2(exp(x) + 1);
    if (al->hes)
      *al->hes = -(exp(x) * (exp(x) - 1)) / gsl_pow_3(exp(x) + 1);
  }
  return check_result(al, gsl_sf_fermi_dirac_m1(x), "gsl_sf_fermi_dirac_m1");
}

static double amplgsl_sf_fermi_dirac_0(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = gsl_sf_fermi_dirac_m1(x);
    if (al->hes)
      *al->hes = exp(x) / gsl_pow_2(exp(x) + 1);
  }
  return check_result(al, gsl_sf_fermi_dirac_0(x), "gsl_sf_fermi_dirac_0");
}

static double amplgsl_sf_fermi_dirac_1(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = gsl_sf_fermi_dirac_0(x);
    if (al->hes)
      *al->hes = gsl_sf_fermi_dirac_m1(x);
  }
  return check_result(al, gsl_sf_fermi_dirac_1(x), "gsl_sf_fermi_dirac_1");
}

static double amplgsl_sf_fermi_dirac_2(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = gsl_sf_fermi_dirac_1(x);
    if (al->hes)
      *al->hes = gsl_sf_fermi_dirac_0(x);
  }
  return check_result(al, gsl_sf_fermi_dirac_2(x), "gsl_sf_fermi_dirac_2");
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
  return check_result(al, gsl_sf_fermi_dirac_int(j, x),
      "gsl_sf_fermi_dirac_int");
}

static double amplgsl_sf_fermi_dirac_mhalf(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al, gsl_sf_fermi_dirac_mhalf(x),
      "gsl_sf_fermi_dirac_mhalf");
}

static double amplgsl_sf_fermi_dirac_half(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al, gsl_sf_fermi_dirac_half(x),
      "gsl_sf_fermi_dirac_half");
}

static double amplgsl_sf_fermi_dirac_3half(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = gsl_sf_fermi_dirac_half(x);
    if (al->hes)
      *al->hes = gsl_sf_fermi_dirac_mhalf(x);
  }
  return check_result(al, gsl_sf_fermi_dirac_3half(x),
      "gsl_sf_fermi_dirac_3half");
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
  return check_result(al, gsl_sf_fermi_dirac_inc_0(x, b),
      "gsl_sf_fermi_dirac_inc_0");
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
  return check_result(al, gamma, "gsl_sf_gamma");
}

static double amplgsl_sf_lngamma(arglist *al) {
  double x = al->ra[0];
  if (al->derivs) {
    *al->derivs = x >= 0 || round(x) != x ? gsl_sf_psi(x) : GSL_NAN;
    if (al->hes)
      *al->hes = gsl_sf_psi_1(x);
  }
  return check_result(al, gsl_sf_lngamma(x), "gsl_sf_lngamma");
}

static double amplgsl_sf_gammastar(arglist *al) {
  double x = al->ra[0];
  double gammastar = gsl_sf_gammastar(x);
  if (al->derivs) {
    double coef = (0.5 / x - log(x) + gsl_sf_psi(x));
    *al->derivs = coef * gammastar;
    if (al->hes) {
      *al->hes = coef * *al->derivs +
          (gsl_sf_psi_1(x) - (1 + 0.5 / x) / x ) * gammastar;
    }
  }
  return check_result(al, gammastar, "gsl_sf_gammastar");
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
  return check_result(al, gammainv, "gsl_sf_gammainv");
}

static double amplgsl_sf_poch(arglist *al) {
  double a = al->ra[0], x = al->ra[1];
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al, gsl_sf_poch(a, x), "gsl_sf_poch");
}

static double amplgsl_sf_lnpoch(arglist *al) {
  double a = al->ra[0], x = al->ra[1];
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al, gsl_sf_lnpoch(a, x), "gsl_sf_lnpoch");
}

static double amplgsl_sf_pochrel(arglist *al) {
  double a = al->ra[0], x = al->ra[1];
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al, gsl_sf_pochrel(a, x), "gsl_sf_pochrel");
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
  return check_result(al, gsl_sf_gamma_inc(a, x), "gsl_sf_gamma_inc");
}

static double amplgsl_sf_gamma_inc_Q(arglist *al) {
  double a = al->ra[0], x = al->ra[1];
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al, gsl_sf_gamma_inc_Q(a, x), "gsl_sf_gamma_inc_Q");
}

static double amplgsl_sf_gamma_inc_P(arglist *al) {
  double a = al->ra[0], x = al->ra[1];
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al, gsl_sf_gamma_inc_P(a, x), "gsl_sf_gamma_inc_P");
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
  return check_result(al, beta, "gsl_sf_beta");
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
  return check_result(al, gsl_sf_lnbeta(a, b), "gsl_sf_lnbeta");
}

static double amplgsl_sf_beta_inc(arglist *al) {
  double a = al->ra[0], b = al->ra[1], x = al->ra[2];
  if (al->derivs) {
    error(al, DERIVS_NOT_PROVIDED);
    return 0;
  }
  return check_result(al, gsl_sf_beta_inc(a, b, x), "gsl_sf_beta_inc");
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
  return check_result(al, gsl_sf_gegenpoly_1(lambda, x),
      "gsl_sf_gegenpoly_1");
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
  return check_result(al, gsl_sf_gegenpoly_2(lambda, x),
      "gsl_sf_gegenpoly_2");
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
  return check_result(al, gsl_sf_gegenpoly_3(lambda, x),
      "gsl_sf_gegenpoly_3");
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
  return check_result(al, gsl_sf_gegenpoly_n(n, lambda, x),
      "gsl_sf_gegenpoly_n");
}

static double amplgsl_sf_hyperg_0F1(arglist *al) {
  double c = al->ra[0], x = al->ra[1];
  if (al->derivs) {
    if (!check_const_arg(al, 0, "c"))
      return 0;
    al->derivs[1] = gsl_sf_hyperg_0F1(c + 1, x) / c;
    if (al->hes)
      al->hes[2] = gsl_sf_hyperg_0F1(c + 2, x) / (c * (c + 1));
  }
  return check_result(al, gsl_sf_hyperg_0F1(c, x), "gsl_sf_hyperg_0F1");
}

static double amplgsl_sf_hyperg_1F1_int(arglist *al) {
  int m = (int)al->ra[0], n = (int)al->ra[1];
  double x = al->ra[2];
  if (!check_int_arg(al, 0, "m") || !check_int_arg(al, 1, "n"))
    return 0;
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
  return check_result(al, gsl_sf_hyperg_1F1_int(m, n, x),
      "gsl_sf_hyperg_1F1_int");
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

  /* Legendre Form of Complete Elliptic Integrals */
  addfunc("gsl_sf_ellint_Kcomp", amplgsl_sf_ellint_Kcomp,
      FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_ellint_Ecomp", amplgsl_sf_ellint_Ecomp,
      FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_ellint_Pcomp", amplgsl_sf_ellint_Pcomp,
      FUNCADD_REAL_VALUED, 2, 0);

  /* Legendre Form of Incomplete Elliptic Integrals */
  addfunc("gsl_sf_ellint_F", amplgsl_sf_ellint_F, FUNCADD_REAL_VALUED, 2, 0);
  addfunc("gsl_sf_ellint_E", amplgsl_sf_ellint_E, FUNCADD_REAL_VALUED, 2, 0);
  addfunc("gsl_sf_ellint_P", amplgsl_sf_ellint_P, FUNCADD_REAL_VALUED, 3, 0);
  addfunc("gsl_sf_ellint_D", amplgsl_sf_ellint_D, FUNCADD_REAL_VALUED, 3, 0);

  /* Carlson Forms */
  addfunc("gsl_sf_ellint_RC", amplgsl_sf_ellint_RC, FUNCADD_REAL_VALUED, 2, 0);
  addfunc("gsl_sf_ellint_RD", amplgsl_sf_ellint_RD, FUNCADD_REAL_VALUED, 3, 0);
  addfunc("gsl_sf_ellint_RF", amplgsl_sf_ellint_RF, FUNCADD_REAL_VALUED, 3, 0);
  addfunc("gsl_sf_ellint_RJ", amplgsl_sf_ellint_RJ, FUNCADD_REAL_VALUED, 4, 0);

  /* Elliptic Functions (Jacobi) */
  /* Wrapper for gsl_sf_elljac_e is not provided since the latter produces
     multiple values (through output parameters). */

  /* Error Functions */
  addfunc("gsl_sf_erf", amplgsl_sf_erf, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_erfc", amplgsl_sf_erfc, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_log_erfc", amplgsl_sf_log_erfc, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_erf_Z", amplgsl_sf_erf_Z, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_erf_Q", amplgsl_sf_erf_Q, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_hazard", amplgsl_sf_hazard, FUNCADD_REAL_VALUED, 1, 0);

  /* Exponential Integrals */
  addfunc("gsl_sf_expint_E1", amplgsl_sf_expint_E1, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_expint_E2", amplgsl_sf_expint_E2, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_expint_En", amplgsl_sf_expint_En, FUNCADD_REAL_VALUED, 2, 0);
  addfunc("gsl_sf_expint_Ei", amplgsl_sf_expint_Ei, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_Shi", amplgsl_sf_Shi, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_Chi", amplgsl_sf_Chi, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_expint_3", amplgsl_sf_expint_3, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_Si", amplgsl_sf_Si, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_Ci", amplgsl_sf_Ci, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_atanint", amplgsl_sf_atanint, FUNCADD_REAL_VALUED, 1, 0);

  /* Fermi-Dirac Function */
  addfunc("gsl_sf_fermi_dirac_m1", amplgsl_sf_fermi_dirac_m1,
      FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_fermi_dirac_0", amplgsl_sf_fermi_dirac_0,
      FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_fermi_dirac_1", amplgsl_sf_fermi_dirac_1,
      FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_fermi_dirac_2", amplgsl_sf_fermi_dirac_2,
      FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_fermi_dirac_int", amplgsl_sf_fermi_dirac_int,
      FUNCADD_REAL_VALUED, 2, 0);
  addfunc("gsl_sf_fermi_dirac_mhalf", amplgsl_sf_fermi_dirac_mhalf,
      FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_fermi_dirac_half", amplgsl_sf_fermi_dirac_half,
      FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_fermi_dirac_3half", amplgsl_sf_fermi_dirac_3half,
      FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_fermi_dirac_inc_0", amplgsl_sf_fermi_dirac_inc_0,
      FUNCADD_REAL_VALUED, 2, 0);

  /* Gamma Functions */
  addfunc("gsl_sf_gamma", amplgsl_sf_gamma, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_lngamma", amplgsl_sf_lngamma, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_gammastar", amplgsl_sf_gammastar, FUNCADD_REAL_VALUED, 1, 0);
  addfunc("gsl_sf_gammainv", amplgsl_sf_gammainv, FUNCADD_REAL_VALUED, 1, 0);

  /* Wrapper for factorials are not provided since these are easily
     implemented using built-in AMPL features like the prod operator. */

  /* Pochhammer Symbol */
  addfunc("gsl_sf_poch", amplgsl_sf_poch, FUNCADD_REAL_VALUED, 2, 0);
  addfunc("gsl_sf_lnpoch", amplgsl_sf_lnpoch, FUNCADD_REAL_VALUED, 2, 0);
  addfunc("gsl_sf_pochrel", amplgsl_sf_pochrel, FUNCADD_REAL_VALUED, 2, 0);

  /* Incomplete Gamma Functions */
  addfunc("gsl_sf_gamma_inc", amplgsl_sf_gamma_inc, FUNCADD_REAL_VALUED, 2, 0);
  addfunc("gsl_sf_gamma_inc_Q", amplgsl_sf_gamma_inc_Q,
      FUNCADD_REAL_VALUED, 2, 0);
  addfunc("gsl_sf_gamma_inc_P", amplgsl_sf_gamma_inc_P,
      FUNCADD_REAL_VALUED, 2, 0);

  /* Beta Functions */
  addfunc("gsl_sf_beta", amplgsl_sf_beta, FUNCADD_REAL_VALUED, 2, 0);
  addfunc("gsl_sf_lnbeta", amplgsl_sf_lnbeta, FUNCADD_REAL_VALUED, 2, 0);

  /* Incomplete Beta Function */
  addfunc("gsl_sf_beta_inc", amplgsl_sf_beta_inc, FUNCADD_REAL_VALUED, 3, 0);

  /* Gegenbauer Functions */
  addfunc("gsl_sf_gegenpoly_1", amplgsl_sf_gegenpoly_1,
      FUNCADD_REAL_VALUED, 2, 0);
  addfunc("gsl_sf_gegenpoly_2", amplgsl_sf_gegenpoly_2,
      FUNCADD_REAL_VALUED, 2, 0);
  addfunc("gsl_sf_gegenpoly_3", amplgsl_sf_gegenpoly_3,
      FUNCADD_REAL_VALUED, 2, 0);
  addfunc("gsl_sf_gegenpoly_n", amplgsl_sf_gegenpoly_n,
      FUNCADD_REAL_VALUED, 3, 0);

  /* Hypergeometric Functions */
  addfunc("gsl_sf_hyperg_0F1", amplgsl_sf_hyperg_0F1,
      FUNCADD_REAL_VALUED, 2, 0);
  addfunc("gsl_sf_hyperg_1F1_int", amplgsl_sf_hyperg_1F1_int,
      FUNCADD_REAL_VALUED, 3, 0);
  // TODO
}
