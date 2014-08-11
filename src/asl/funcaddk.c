/****************************************************************
Copyright (C) 1997-1998 Lucent Technologies
All Rights Reserved

Permission to use, copy, modify, and distribute this software and
its documentation for any purpose and without fee is hereby
granted, provided that the above copyright notice appear in all
copies and that both that the copyright notice and this
permission notice and warranty disclaimer appear in supporting
documentation, and that the name of Lucent or any of its entities
not be used in advertising or publicity pertaining to
distribution of the software without specific, written prior
permission.

LUCENT DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.
IN NO EVENT SHALL LUCENT OR ANY OF ITS ENTITIES BE LIABLE FOR ANY
SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER
IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,
ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
THIS SOFTWARE.
****************************************************************/

/* sample funcadd, K&R C variant */

#ifndef KR_headers
#define KR_headers
#endif

#include "math.h"	/* for sqrt */
#include "funcadd.h"	/* includes stdio1.h */

 char *ix_details_ASL[] = {0};	/* no -i command-line option */

 static AmplExports *AMPL_Exp;	/* for fprintf and strtod */

 static real
ginv(al) arglist *al;	/* generalized inverse of a single argument */
{
	real x = al->ra[0];
	x = x ? 1./x : 0.;
	if (al->derivs)
		*(real *)al->derivs = -x*x;
	return x;
	}

 static real
myhypot(al) arglist *al;	/* sqrt(x*x + y*y) */
{
	real *d, rv, x, x0, y, y0;

	x = x0 = al->ra[0];
	y = y0 = al->ra[1];

	if (x < 0.)
		x = -x;
	if (y < 0.)
		y = -y;
	rv = x;
	if (x < y) {
		rv = y;
		y = x;
		x = rv;
		}
	if (rv) {
		y /= x;
		rv = x * sqrt(1. + y*y);
		if (d = (real *)al->derivs) {
			d[0] = x0 / rv;
			d[1] = y0 / rv;
			}
		}
	else if (d = (real *)al->derivs)
		d[0] = d[1] = 0;
	return rv;
	}

 static real
mean(al) arglist *al;	/* mean of arbitrarily many arguments */
{
	real x, z;
	real *d, *de, *ra;
	int *at, i, j, n;
	char *se, *sym;
	AmplExports *ae = al->AE; /* for fprintf and strtod */

	if ((n = al->n) <= 0)
		return 0;
	at = al->at;
	ra = al->ra;
	d = de = (real *)al->derivs;
	x = 0.;
	for(i = 0; i < n;)
		if ((j = at[i++]) >= 0) {
			x += ra[j];
			++de;
			}
		else {
			x += strtod(sym = al->sa[-(j+1)], &se);
			if (*se) {
				fprintf(Stderr,
				"mean treating arg %d = \"%s\" as 0\n",
					i, sym);
				fflush(Stderr);
				}
			}
	if (d) {
		z = 1. / n;
		while(d < de)
			*d++ = z;
		}
	return x / n;
	}

 void
funcadd(ae) AmplExports *ae;
{
	AMPL_Exp = ae;
	/* Insert calls on addfunc here... */
	/* Arg 3, called argtype, can be 0 or 1:
	 *	0 ==> force all arguments to be numeric
	 *	1 ==> pass both symbolic and numeric arguments.
	 *
	 * Arg 4, called nargs, is interpretted as follows:
	 *	>=  0 ==> the function has exactly nargs arguments
	 *	<= -1 ==> the function has >= -(nargs+1) arguments.
	 *
	 * Arg 5, funcinfo, is passed to the functions in struct arglist;
	 *	it is not used in these examples, so we just pass 0.
	 */

	addfunc("ginv", (ufunc*)ginv, 0, 1, 0);
	addfunc("hypot", (ufunc*)myhypot, 0, 2, 0);
	addfunc("mean", (ufunc*)mean, 1, -1, 0);
	}
