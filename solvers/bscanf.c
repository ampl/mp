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

#ifdef KR_headers
#include "varargs.h"
#else
#include "stdarg.h"
#endif

#include "asl.h"	/* for Long */

#ifndef Int
#define Int Long
#endif

 static char *
#ifdef KR_headers
Ladvance(s, Lp)
 char *s;
 Long *Lp;
#else
Ladvance(register char *s, Long *Lp)
#endif
{
	register Long x;
	register int c;

	int sign;
	c = *s++;
	if (c == '-') {
		sign = 1;
		c = *s++;
		}
	else
		sign = 0;
	if (c < '0' || c > '9')
		return 0;
	for(x = c - '0';; s++) {
		c = *s;
		if (c < '0' || c > '9')
			break;
		x = 10*x + c - '0';
		}
	*Lp = sign ? -x : x;
	return s;
	}

 static void
#ifdef KR_headers
badfmt(R, fmt0) EdRead *R; char *fmt0;
#else
badfmt(EdRead *R, const char *fmt0)
#endif
{
	badread(R);
	fprintf(Stderr, "bad format %s\n", fmt0);
	mainexit_ASL(1);
	}

 int
#ifdef KR_headers
ascanf(va_alist)
 va_dcl
#else
ascanf(EdRead *R, const char *fmt, ...)
#endif
{
	Const char *fmt0;
	int *ip;
	Long L, *Lp;
	double x, *xp;
	va_list ap;
	char *s, *s1;
	int rc = 0;

#ifdef KR_headers
	EdRead *R;
	char *fmt;
	va_start(ap);
	R = va_arg(ap, EdRead*);
	fmt = va_arg(ap, char*);
#else
	va_start(ap, fmt);
#endif
	s = read_line(R);
	if (!s)
		return 0;
	for(;;) {
		fmt0 = fmt;
		if (*fmt++ != '%')
			break;
		while(*s == ' ')
			s++;
		switch(*fmt++) {
			case 'd':
				if (!(s = Ladvance(s,&L)))
					return rc;
				ip = va_arg(ap, int*);
				*ip = (int)L;
				break;
			case 'l':
				switch(*fmt++) {
					case 'd':
						if (!(s = Ladvance(s,&L)))
							return rc;
						Lp = va_arg(ap, Long*);
						*Lp = (int)L;
						break;
					case 'f':
						x = strtod(s, &s);
						if (!s)
							return rc;
						xp = va_arg(ap, double*);
						*xp = x;
						break;
					default:
						badfmt(R, fmt0);
					}
				break;
			default:
				if (!(fmt = Ladvance((char *)fmt-1, &L))
				 || *fmt++ != 's')
					badfmt(R, fmt0);
				/* %127s */
				s1 = va_arg(ap, char*);
				while(--L > 0 && (*s1 = *s++))
					s1++;
				*s1 = 0;
			}
		rc++;
		while(*fmt == ' ')
			fmt++;
		}
	return rc;
	}

 int
#ifdef KR_headers
bscanf(va_alist)
 va_dcl
#else
bscanf(EdRead *R, const char *fmt, ...)
#endif
{
	Const char *fmt0;
	int *ip;
	short *shp;
	Long L, *Lp;
	Int I;
	double *xp;
	va_list ap;
	char *s1;
	int len;
	FILE *fd;
	int rc = 0;

#ifdef KR_headers
	EdRead *R;
	char *fmt;
	va_start(ap);
	R = va_arg(ap, EdRead*);
	fmt = va_arg(ap, char*);
#else
	va_start(ap, fmt);
#endif
	fd = R->nl;
	R->Line += R->lineinc;
	R->lineinc = 1;
	for(;;) {
		fmt0 = fmt;
		if (*fmt++ != '%')
			break;
		switch(*fmt++) {
			case 'd':
				ip = va_arg(ap, int*);
				if (!fread(&I, sizeof(Int), 1, fd))
					return rc;
				*ip = (int)I;
				if (R->iadjfcn)
					(*R->iadjfcn)(ip, sizeof(int));
				break;
			case 'h':
				if (*fmt == 'd')
					fmt++;
				shp = va_arg(ap, short*);
				if (!fread(shp, sizeof(short), 1, fd))
					return rc;
				if (R->iadjfcn)
					(*R->iadjfcn)(shp, sizeof(short));
				break;
			case 'l':
				switch(*fmt++) {
					case 'd':
						Lp = va_arg(ap, Long*);
						if (!fread(Lp, sizeof(Long),
								1, fd))
							return rc;
						if (R->iadjfcn)
							(*R->iadjfcn)(Lp,
								sizeof(Long));

						break;
					case 'f':
						xp = va_arg(ap, double*);
						if (!fread(xp, sizeof(double), 1, fd))
							return rc;
						if (R->dadjfcn)
							(*R->dadjfcn)(xp,
								sizeof(double));
						break;
					default:
						badfmt(R, fmt0);
					}
				break;
			default:
				if (!(fmt = Ladvance((char *)fmt-1, &L))
				 || *fmt++ != 's')
					badfmt(R, fmt0);
				/* %127s */
				s1 = va_arg(ap, char*);
				if (!fread(&len, sizeof(int), 1, fd))
					return rc;
				if (R->iadjfcn)
					(*R->iadjfcn)(&len, sizeof(int));
				if (len >= L || !fread(s1, len, 1, fd))
					return rc;
				s1[len] = 0;
				break;
			}
		rc++;
		while(*fmt == ' ')
			fmt++;
		}
	return rc;
	}
