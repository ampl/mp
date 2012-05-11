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

#include "math.h"
#include "errno.h"
#include "limits.h"

#ifdef __cplusplus
extern "C" {
#endif

 double
#ifdef KR_headers
mypow_ASL(x,y) double x, y;
#else
mypow_ASL(double x, double y) /* return x ^ y (exponentiation) */
#endif
{
	double xy, y1, ye;
	unsigned long i;
	int ex, ey, flip;

	if (!y)
		return 1.0;	/* Kahan advocates this */

	flip = 0;
	if (y < 0.)
		{ y = -y; flip = 1; }
	if (y1 = modf(y, &ye)) {
		if (x <= 0.)
			goto zreturn;
		if (y1 > 0.5) {
			y1 -= 1.;
			ye += 1.;
			}
		xy = exp(y1 * log(x));
		}
	else {
#ifdef __TURBOC__
		if (!x)
			goto zreturn;
#endif
		xy = 1.0;
		}
	if (ye > (unsigned long)ULONG_MAX) {
		if (x <= 0.) {
 zreturn:
			if (x || flip)
				errno = EDOM;
			return 0.;
			}
		if (flip)
			y = -y;
		return exp(y * log(x));
		}
	x = frexp(x, &ex);
	ey = 0;
	if (i = (unsigned long)ye) for(;;) {
		if (i & 1) { xy *= x; ey += ex; }
		if (!(i >>= 1))
			break;
		x *= x;
		ex <<= 1;
		if (x < .5)
			{ x += x; ex -= 1; }
		}
	if (flip)
		{ xy = 1. / xy; ey = -ey; }
	errno = 0;
	x = ldexp(xy, ey);
	if (errno && ey < 0) {
		errno = 0;
		x = 0.;
		}
	return x;
	}
#ifdef __cplusplus
}
#endif
