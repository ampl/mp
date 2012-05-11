/****************************************************************
Copyright (C) 1997, 2000 Lucent Technologies
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

/* sample funcadd */

#include "stdlib.h"	/* for atoi */
#include "math.h"	/* for sqrt */
#include "funcadd.h"	/* includes "stdio1.h" */

 static real
ginv(arglist *al)	/* generalized inverse of a single argument */
{
	real x = al->ra[0];
	x = x ? 1./x : 0.;
	if (al->derivs) {
		*al->derivs = -x*x;
		if (al->hes)
			*al->hes = -2.*x * *al->derivs;
		}
	return x;
	}

 static char *
sginv(arglist *al)	/* character-valued version of ginv */
{
	AmplExports *ae = al->AE;	/* for sprintf and TempMem */
	real x = al->ra[0];
	/* Get memory for the return value from TempMem. */
	/* This memory is automatically freed after the return */
	/* value has been processed. */
	/* A nonreentrant alternative would be to declare */
	/*	static char buf[32];	*/
	char *buf = (char*)TempMem(al->TMI, 32);
	sprintf(buf, "x%.g", x ? 1./x : 0);
	return buf;
	}

 static real
myhypot(arglist *al)	/* myhypot(x,y) = sqrt(x*x + y*y) */
{
	real *d, *h, rv, x, x0, y, y0;

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
		if (d = al->derivs) {
			d[0] = x0 / rv;
			d[1] = y0 / rv;
			if (h = al->hes) {
				h[0] =  d[1]*d[1] / rv;
				h[1] = -d[0]*d[1] / rv;
				h[2] =  d[0]*d[0] / rv;
				}
			}
		}
	else if (d = al->derivs) {
		d[0] = d[1] = 0;
		if (h = al->hes)
			h[0] = h[1] = h[2] = 0;
		}
	return rv;
	}

 static real
ncall(void)		/* returns its invocation count */
{ static real x; return ++x; }

 static real
mean(arglist *al)	/* mean of arbitrarily many arguments */
{
	real x, z;
	real *d, *de, *ra;
	int *at, i, j, n;
	char *se;
	Const char *sym;
	AmplExports *ae = al->AE; /* for fprintf and strtod */

	if ((n = al->n) <= 0)
		return 0;
	at = al->at;
	ra = al->ra;
	d = de = al->derivs;
	x = 0.;
	for(i = 0; i < n;)
		if ((j = at[i++]) >= 0) {
			x += ra[j];
			++de;
			}
		else {
			x += z = strtod(sym = al->sa[-(j+1)], &se);
			if (*se) {
				fprintf(Stderr,
				"mean treating arg %d = \"%s\" as %.g\n",
					i, sym, z);
				/* Stderr may be stdout on some systems, */
				/* so flushing it is recommended. */
				fflush(Stderr);
				}
			}
	if (d) {
		z = 1. / n;
		while(d < de)
			*d++ = z;
		/* The Hessian is == 0 and, if needed, has been */
		/* initialized to 0. */
		}
	return x / n;
	}

/* Sample function kth optionally prints its arguments using
 * a variant of stdio.h supplied by funcadd.h. */

 static Const char *
kth(arglist *al)	/* kth(k,a1,a2,...,an) return ak */
{
	int j, k, n;
	char *buf;
	AmplExports *ae = al->AE;

	k = al->at[0] ? atoi(al->sa[0]) : (int)al->ra[0];
	n = al->n;
	if (k < 0) {
		al->Errmsg = buf = (char*)TempMem(al->TMI, 64);
		sprintf(buf, "kth(k,...) has k = %d < 0", k);
		return buf;
		}
	if (n <= 1 || k <= 0 || k >= n)
		return "";
	if ((j = al->at[k]) >= 0) {
		buf = (char*)TempMem(al->TMI, 32);
		sprintf(buf, "%.g", al->ra[j]);
		return buf;
		}
	return al->sa[-(j+1)];
	}

 static char *
get_env(arglist *al)	/* test/demonstrate ae->Getenv */
{
	AmplExports *ae = al->AE;

	if (al->at[0])
		return getenv(al->sa[0]);
	return "";
	}

/* Illustration of at_exit() and at_reset() processing */

 typedef struct Aeinfo
{
	AmplExports *ae;
	int n;
	} Aeinfo;

static Aeinfo AEI[10], *AEInext = AEI, *AEIlast = AEI + 10;

 static void
At_end(void *v)
{
	Aeinfo *aei = (Aeinfo *)v;
	AmplExports *ae = aei->ae;
	printf("Got to At_end: n = %d\n", aei->n);
	}

 static void
At_reset(void *v)
{
	Aeinfo *aei = (Aeinfo *)v;
	AmplExports *ae = aei->ae;
	printf("Got to At_reset: n = %d\n", aei->n);
	}

static real
ginvae(arglist *al)	/* like ginv, but enrolling At_reset and At_exit */
{
	static int nginv;
	AmplExports *ae = al->AE;

	if (AEInext < AEIlast) {
		AEInext->n = ++nginv;
		AEInext->ae = ae;
		at_reset(At_reset, AEInext);
		at_exit(At_end, AEInext++);
		}
	return ginv(al);
	}

 void
funcadd(AmplExports *ae){
	/* Insert calls on addfunc here... */

/* Arg 3, called type, must satisfy 0 <= type <= 6:
 * type&1 == 0:	0,2,4,6	==> force all arguments to be numeric.
 * type&1 == 1:	1,3,5	==> pass both symbolic and numeric arguments.
 * type&6 == 0:	0,1	==> the function is real valued.
 * type&6 == 2:	2,3	==> the function is char * valued; static storage
			    suffices: AMPL copies the return value.
 * type&6 == 4:	4,5	==> the function is random (real valued).
 * type&6 == 6: 6	==> random, real valued, pass nargs real args,
 *				0 <= nargs <= 2.
 *
 * Arg 4, called nargs, is interpretted as follows:
 *	>=  0 ==> the function has exactly nargs arguments
 *	<= -1 ==> the function has >= -(nargs+1) arguments.
 *
 * Arg 5, called funcinfo, is copied without change to the arglist
 *	structure passed to the function; funcinfo is for the
 *	function to use or ignore as it sees fit.
 */

	/* Solvers quietly ignore kth, sginv, and rncall, since */
	/* kth and sginv are symbolic (i.e., char* valued) and  */
	/* rncall is specified as random.  Thus kth, sginv, and */
	/* rncall may not appear nonlinearly in declarations in */
	/* an AMPL model. */

	addfunc("ginv", (rfunc)ginv, 0, 1, 0);
	addfunc("sginv", (rfunc)sginv, 2, 1, 0);
	addfunc("hypot", (rfunc)myhypot, 0, 2, 0);
	addfunc("ncall", (rfunc)ncall, 0, 0, 0);
	addfunc("rncall", (rfunc)ncall, 4, 0, 0);    /* could change 4 to 6 */
	addfunc("mean0", (rfunc)mean, 0, -1, 0);
	addfunc("mean", (rfunc)mean, 1, -1, 0);
	addfunc("kth", (rfunc)kth, 3, -2, 0);
	if (ae->ASLdate >= 20000608)
		addfunc("getenv", (rfunc)get_env, 3, 1, 0);
	addfunc("ginvae", (rfunc)ginvae, 0, 1, 0); /* demo at_exit, at_reset */
	/* at_end() and at_reset() calls could appear here, too. */
	}
