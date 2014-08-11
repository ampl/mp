/****************************************************************
Copyright (C) 1997 Lucent Technologies
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

/* g_fmt(buf,x) stores the closest decimal approximation to x in buf;
 * it suffices to declare buf
 *	char buf[32];
 */

/* For systems (e.g. Crays) with such poor arithmetic that
 * dtoa does not work on them, #define No_dtoa (or compile with
 * -DNo_dtoa in CFLAGS).
 */

#include "arith.h"
#include "asl.h"

 char g_fmt_E = 'e';
 int g_fmt_decpt;

 int
#ifdef KR_headers
g_fmtp(b, x, prec) register char *b; double x; int prec;
#else
g_fmtp(register char *b, double x, int prec)
#endif
{
	char *b0 = b;
#ifdef No_dtoa
#ifndef g_fmt_prec
#define g_fmt_prec 15
#endif
	sprintf(b, "%.*g", prec ? prec : g_fmt_prec, x);
	while(*b)
		b++;
#else
	register int i, k;
	register char *s;
	int decpt, j, sign;
	char *s0, *se;

	if (!x) {
		*b++ = '0';
		if (g_fmt_decpt) {
			*b++ = '.';
			if (g_fmt_decpt == 2) {
				b[0] = g_fmt_E;
				b[1] = '+';
				b[2] = '0';
				b[3] = '0';
				b += 4;
				}
			}
		goto done;
		}
	s = s0 = dtoa(x, prec ? 2 : 0, prec, &decpt, &sign, &se);
	if (sign)
		*b++ = '-';
	if (decpt == 9999) /* Infinity or Nan */ {
		while(*b = *s++)
			b++;
		goto done0;
		}
	if (decpt <= -4 || decpt > se - s + 5 || g_fmt_decpt == 2) {
		*b++ = *s++;
		if (*s || g_fmt_decpt) {
			*b++ = '.';
			while(*b = *s++)
				b++;
			}
		*b++ = g_fmt_E;
		/* sprintf(b, "%+.2d", decpt - 1); */
		if (--decpt < 0) {
			*b++ = '-';
			decpt = -decpt;
			}
		else
			*b++ = '+';
		for(j = 2, k = 10; 10*k <= decpt; j++, k *= 10);
		for(;;) {
			i = decpt / k;
			*b++ = i + '0';
			if (--j <= 0)
				break;
			decpt -= i*k;
			decpt *= 10;
			}
		}
	else if (decpt <= 0) {
		*b++ = '0';
		*b++ = '.';
		for(; decpt < 0; decpt++)
			*b++ = '0';
		while(*b = *s++)
			b++;
		}
	else {
		while(*b = *s++) {
			b++;
			if (--decpt == 0 && (g_fmt_decpt || *s))
				*b++ = '.';
			}
		if (decpt > 0) {
			do *b++ = '0';
				while(--decpt);
			if (g_fmt_decpt)
				*b++ = '.';
			}
		}
 done0:
	freedtoa(s0);
 done:
	*b = 0;
#endif
	return b - b0;
	}

 int
#ifdef KR_headers
g_fmt(b, x) register char *b; double x;
#else
g_fmt(register char *b, double x)
#endif
{ return g_fmtp(b, x, 0); }

#ifdef KR_headers
 extern int obj_prec();

 int
g_fmtop(b, x) register char *b; double x;
#else
 extern int obj_prec(void);

 int
g_fmtop(register char *b, double x)
#endif
{
	return g_fmtp(b, x, obj_prec());
	}
