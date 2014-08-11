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

#include "stdarg.h"

#include "asl.h"	/* for Long */

#ifndef Int
#define Int Long
#endif

 static char *
Ladvance(register char *s, Long *Lp)
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
badfmt(EdRead *R, const char *fmt0)
{
	badread(R);
	fprintf(Stderr, "bad format %s\n", fmt0);
	mainexit_ASL(1);
	}

 int
ascanf(EdRead *R, const char *fmt, ...)
{
	Const char *fmt0;
	Long L, *Lp;
	char *s, *s1;
	double x, *xp;
	int *ip, rc;
	va_list ap;

	va_start(ap, fmt);
	rc = 0;
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
				*ip = L;
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
						x = strtod(s, &s1);
						if (s1 <= s)
							return rc;
						xp = va_arg(ap, double*);
						*xp = x;
						s = s1;
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
bscanf(EdRead *R, const char *fmt, ...)
{
	Const char *fmt0;
	FILE *fd;
	Int I;
	Long L, *Lp;
	char *s1;
	double *xp;
	int *ip, len, rc;
	short sh;
	va_list ap;

	va_start(ap, fmt);
	rc = 0;
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
				if (!fread(&I, sizeof(I), 1, fd))
					return rc;
				if (R->iadjfcn)
					(*R->iadjfcn)(&I, sizeof(I));
				*ip = I;
				break;
			case 'h':
				if (*fmt == 'd')
					fmt++;
				Lp = va_arg(ap, Long*);
				if (!fread(&sh, sizeof(short), 1, fd))
					return rc;
				if (R->iadjfcn)
					(*R->iadjfcn)(&sh, sizeof(short));
				*Lp = sh;
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

 int
hscanf(EdRead *R, const char *fmt, ...)
{
	Const char *fmt0;
	FILE *fd;
#ifdef NO_LONG_LONG
	Long I[2];
#else
	long long I;
#endif
	Long L, *Lp;
	char *s1;
	double *xp;
	int *ip, len, rc;
	short sh;
	va_list ap;

	va_start(ap, fmt);
	rc = 0;
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
#ifdef NO_LONG_LONG /*{{*/
				if (!fread(I, sizeof(I), 1, fd))
					return rc;
				if (R->iadjfcn)
					(*R->iadjfcn)(I, sizeof(I));
#ifdef IEEE_8087
				*ip = I[0];
#else
				*ip = I[1];
#endif
#else /*}{*/
				if (!fread(&I, sizeof(I), 1, fd))
					return rc;
				if (R->iadjfcn)
					(*R->iadjfcn)(&I, sizeof(I));
				*ip = (int)I;
#endif /*}}*/
				break;
			case 'h':
				if (*fmt == 'd')
					fmt++;
				Lp = va_arg(ap, Long*);
				if (!fread(&sh, sizeof(short), 1, fd))
					return rc;
				if (R->iadjfcn)
					(*R->iadjfcn)(&sh, sizeof(short));
				*Lp = sh;
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
