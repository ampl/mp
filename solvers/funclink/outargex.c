/****************************************************************
Copyright (C) 2002 AMPL Optimization LLC
All Rights Reserved

Permission to use, copy, modify, and distribute this software and its
documentation for any purpose and without fee is hereby granted,
provided that the above copyright notice appear in all copies and that
both that the copyright notice and this permission notice and warranty
disclaimer appear in supporting documentation.

AMPL Optimization LLC DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS
SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS.  IN NO EVENT SHALL AMPL Optimization LLC BE LIABLE FOR ANY
SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER
RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
****************************************************************/


/* Illustration of OUT args in an AMPL imported function */

/* Used with

	function outargex(INOUT ...);
	param nr integer > 0 default 4;
	param ns integer > 0 default 3;
	param p{1..nr};
	param s{i in 1..ns} symbolic;
	param sout symbolic;

	# outargex(nr, ns, {i in 1..nr} p[i], {i in 1..ns} s[i], sout)
	# returns the sum of the incoming p[i] values,
	# increments the p[i], sets sout to the concatenation of the
	# incoming s[i] and the return value, and modifies the s[i]
	# by "incrementing" non-numeric string values and interchanging
	# the numeric and string representations of ones that look
	# like numbers.

	# ... assign values to p, s; then

	call outargex(nr, ns, {i in 1..nr} p[i], {i in 1..ns} s[i], sout);
	display p, s, sout;
	display outargex(nr, ns, {i in 1..nr} p[i], {i in 1..ns} s[i], sout);
	display p, s, sout;
*/

#include "funcadd.h"

 static real
outargex(arglist *al)
{
	AmplExports *ae = al->AE;
	int *at, i, len, m, ms, n, nr, ns;
	real *ra, rv, t;
	char *buf, *se;
	const char *s, **sa;

	n = al->n;
	if (n < 3) {
		al->Errmsg = "outargex: too few args";
		return 0;
		}
	at = al->at;
	ra = al->ra;
	sa = al->sa;
	nr = ra[0];
	ns = ra[1];
	rv = 0;
	for(i = 2, m = nr + 2; i < m; i++) {
		if ((at[i] & (AMPLFUNC_INARG|AMPLFUNC_OUTARG|AMPLFUNC_STRING))
			!=   (AMPLFUNC_INARG|AMPLFUNC_OUTARG)) {
			al->Errmsg = buf = (char*)TempMem(al->TMI,60);
			sprintf(buf, "outargex: arg %d is not a real INOUT arg",i);
			return 0;
			}
		rv += ra[i]++;
		}
	/* i == m */
	len = 1; /* for trailing 0 */
	for(ms = m + ns; i < ms; i++) {
		if (!(at[i] & AMPLFUNC_INARG)) {
			al->Errmsg = buf = (char*)TempMem(al->TMI,60);
			sprintf(buf, "outargex: arg %d is not an IN arg",i);
			return 0;
			}
		len += at[i] & AMPLFUNC_STRING ? strlen(sa[i]) + 1: 20;
		/* 20 = upper bound on length contributed by ra[i] */
		}
	if ((at[ms] & (AMPLFUNC_OUTARG|AMPLFUNC_STROUT))
		   != (AMPLFUNC_OUTARG|AMPLFUNC_STROUT)) {
		al->Errmsg = "outargex: final arg is not symbolic OUT";
		return 0;
		}
	len += len;	/* storage for "incrementing" the s[i] */
	len += 20;	/* for concatenating the return value to sout */
	al->sa[ms] = buf = (char*)TempMem(al->TMI, len);
	for(i = m; i < ms; i++) {
		if (at[i] & AMPLFUNC_STRING)
			buf += sprintf(buf, "%s,", sa[i]);
		else
			buf += sprintf(buf, "%.g,", ra[i]);
		}
	buf += sprintf(buf, "%.g", rv); /* concatenate return value */

	/* "Increment" the other incoming char* values */
	/* just so they will be seen to change. */

	for(i = m; i < ms; i++) {
		if (s = sa[i]) {
			t = strtod(s,&se);
			if (!*se) {
				/* A string that looks like a number: */
				/* turn it into a real number. */
				ra[i] = t + 1.;
				sa[i] = 0;
				}
			else {
				strcpy(++buf, s);
				if (*buf >= ' ' && *buf < 0x7e)
					++*buf; /* increment first character */
				sa[i] = buf;
				buf += strlen(buf);
				}
			}
		else {
			/* replace numeric arg with incremented string */
			t = ra[i] + 1.;
			sa[i] = (const char*)buf;
			buf += sprintf(buf, "%.g", t) + 1;
			}
		}
	return rv;
	}

 void
funcadd(AmplExports *ae){
	addfunc("outargex", (rfunc)outargex,
		FUNCADD_STRING_ARGS|FUNCADD_OUTPUT_ARGS, -4, 0);
	}
