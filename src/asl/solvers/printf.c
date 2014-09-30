/****************************************************************
Copyright (C) 1997, 1999, 2001 Lucent Technologies
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

/* This implements most of ANSI C's printf, fprintf, and sprintf,
 * omitting L, with %.0g and %.0G giving the shortest decimal string
 * that rounds to the number being converted, and with negative
 * precisions allowed for %f.
 */

#include "stddef.h"
#include "stdarg.h"
#include "stdlib.h"

#ifndef NO_PRINTF_A_FMT
#include "arith.h"
#endif

#include "stdio1.h"
#include "string.h"
#include "errno.h"

#ifdef KR_headers
#define Const /* const */
#define Voidptr char*
#ifndef size_t__
#define size_t int
#define size_t__
#endif

#else

#define Const const
#define Voidptr void*

#endif

#undef MESS
#ifndef Stderr
#define Stderr stderr
#endif

#ifdef _windows_
#undef PF_BUF
#define MESS
#include "mux0.h"
#define stdout_or_err(f) (f == stdout)
#else
#define stdout_or_err(f) (f == Stderr || f == stdout)
#endif

#ifdef __cplusplus
extern "C" {
#endif

 extern char *dtoa ANSI((double, int, int, int*, int*, char **));
 extern void freedtoa ANSI((char*));



#ifdef USE_ULDIV
/* This is for avoiding 64-bit divisions on the DEC Alpha, since */
/* they are not portable among variants of OSF1 (DEC's Unix). */

#define ULDIV(a,b) uldiv_ASL(a,(unsigned long)(b))

#ifndef LLBITS
#define LLBITS 6
#endif
#ifndef ULONG
#define ULONG unsigned long
#endif

 static int
klog(ULONG x)
{
	int k, rv = 0;

	if (x > 1L)
	    for(k = 1 << LLBITS-1;;) {
		if (x >= (1L << k)) {
			rv |= k;
			x >>= k;
			}
		if (!(k >>= 1))
			break;
		}
	return rv;
	}

 ULONG
uldiv_ASL(ULONG a, ULONG b)
{
	int ka;
	ULONG c, k;
	static ULONG b0;
	static int kb;

	if (a < b)
		return 0;
	if (b != b0) {
		b0 = b;
		kb = klog(b);
		}
	k = 1;
	if ((ka = klog(a) - kb) > 0) {
		k <<= ka;
		b <<= ka;
		}
	c = 0;
	for(;;) {
		if (a >= b) {
			a -= b;
			c |= k;
			}
		if (!(k >>= 1))
			break;
		a <<= 1;
		}
	return c;
	}

#else
#define ULDIV(a,b) a / b
#endif /* USE_ULDIV */

 typedef struct
Finfo {
	union { FILE *cf; char *sf; } u;
	char *ob0, *obe1;
	size_t lastlen;
	} Finfo;

 typedef char *(*Putfunc) ANSI((Finfo*, int*));

#ifdef PF_BUF
FILE *stderr_ASL = (FILE*)&stderr_ASL;
void (*pfbuf_print_ASL) ANSI((char*));
char *pfbuf_ASL;
static char *pfbuf_next;
static size_t pfbuf_len;
extern Char *mymalloc_ASL ANSI((size_t));
extern Char *myralloc_ASL ANSI((void *, size_t));

#undef fflush
#ifdef old_fflush_ASL
#define fflush old_fflush_ASL
#endif

 void
fflush_ASL(FILE *f)
{
	if (f == stderr_ASL) {
		if (pfbuf_ASL && pfbuf_print_ASL) {
			(*pfbuf_print_ASL)(pfbuf_ASL);
			free(pfbuf_ASL);
			pfbuf_ASL = 0;
			}
		}
	else
		fflush(f);
	}

 static void
pf_put(char *buf, int len)
{
	size_t x, y;
	if (!pfbuf_ASL) {
		x = len + 256;
		if (x < 512)
			x = 512;
		pfbuf_ASL = pfbuf_next = (char*)mymalloc_ASL(pfbuf_len = x);
		}
	else if ((y = (pfbuf_next - pfbuf_ASL) + len) >= pfbuf_len) {
		x = pfbuf_len;
		while((x <<= 1) <= y);
		y = pfbuf_next - pfbuf_ASL;
		pfbuf_ASL = (char*)myralloc_ASL(pfbuf_ASL, x);
		pfbuf_next = pfbuf_ASL + y;
		pfbuf_len = x;
		}
	memcpy(pfbuf_next, buf, len);
	pfbuf_next += len;
	*pfbuf_next = 0;
	}

 static char *
pfput(Finfo *f, int *rvp)
{
	int n;
	char *ob0 = f->ob0;
	*rvp += n = (int)(f->obe1 - ob0);
	pf_put(ob0, n);
	return ob0;
	}
#endif /* PF_BUF */

 static char *
Fput ( Finfo *f, int *rvp)
{
	char *ob0 = f->ob0;

	*rvp += f->obe1 - ob0;
	*f->obe1 = 0;
	fputs(ob0, f->u.cf);
	return ob0;
	}


#ifdef _windows_
int stdout_fileno_ASL = 1;

 static char *
Wput(Finfo *f, int *rvp)
{
	char *ob0 = f->ob0;

	*rvp += f->obe1 - ob0;
	*f->obe1 = 0;
	mwrite(ob0, f->obe1 - ob0);
	return ob0;
	}
#endif /*_windows_*/

#ifdef QUOTIFY /*{*/
 static char qtype[256];
 static char dig[256] = {
	13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
	13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
	13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 10, 13, 10, 13, 13,
	 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 13, 13, 13, 13, 13, 13,
	13, 13, 13, 13, 11, 11, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
	13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
	13, 13, 13, 13, 11, 11, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
	13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
	13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
	13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
	13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
	13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
	13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
	13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
	13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
	13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13
	};

 static int
lcstrcmp(const char *a, const char *b)
{
	int c, d;

	for(;;) {
		if ((c = *a++) >= 'A' && c <= 'Z')
			c += 'a' - 'A';
		if (d = c - *b++)
			return d;
		if (!c)
			break;
		}
	return 0;
	}

 static int
valid_param(char *t, int ekeep)
	/* return 0 unless valid decimal fp number */
	/* also adjust Fortran exponent notations to C style */
{
	char *Dig = dig, *e = 0, *se;
	int dc = 0;
	/* digit types: 0 <= d <= 9 ==> d, +- = 10, e,E,d,D = 11 */

#define D(n) Dig[*(unsigned char *)n]

	if (D(t) == 10)
		t++;
	switch(*t) {
	 case '0':
		if ((t[1] == 'x' || t[1] == 'X')) {
			strtod(t,&se);
			return *se == 0;
			}
		break;
	 case 'I':
	 case 'i':
		if (!lcstrcmp(t,"infinity"))
			return 1;
		break;
	 case 'N':
	 case 'n':
		if (Arith_Kind_ASL < 3 && !lcstrcmp(t,"nan"))
			return 1;
	 }
	while(D(t) < 10)
		{ t++; dc++; }
	if (*t == '.')
		while(D(++t) < 10)
			dc++;
	if (D(t) == 11) {
		e = t;
		if (D(++t) == 10)
			t++;
		if (D(t) >= 10)
			return 0;
		while(D(t) < 10)
			t++;
		}
	if (dc && !*t) {
		if (e && !ekeep)
			*e = 'e';
		return 1;
		}
	return 0;
	}
#undef D

 void
qt_init(void)
{
	int n1;
	const char *t;
	for(t = "\"'\n.+-_"; *t; t++)
		qtype[(int)*t] = *t;
	for(n1 = 'A'; n1 <= 'Z'; n1++)
		qtype[n1] = qtype[n1+'a'-'A'] = 1;
	for(n1 = '0'; n1 <= '9'; n1++)
		qtype[n1] = 2;
	}

 static void
qkind(char *s, int prec, int *quote, int *widthp)
{
	int n0, n1, n2, ne;
	char *s0;

	if (!*s)
		goto quote_it;
	if (!qtype['"'])
		qt_init();
	n0 = n1 = n2 = ne = 0;
	s0 = s;
	while(prec-- > 0)
		switch(qtype[*(unsigned char *)s++]) {
			case '\'': n1++; break;
			case '"':  n2++; break;
			case '\n': ne++; break;
			case 0:	   n0++;
			}
	if ((ne + n1 + n2) || n0) {
		if (n1 > n2) {
			ne += n2;
			*quote = '"';
			}
		else {
			ne += n1;
			*quote = '\'';
			}
		*widthp -= ne + 2;
		}
	else if (valid_param(s0,1)) {
 quote_it:
		*widthp -= 2;
		*quote = '\'';
		}
	}
#endif /*QUOTIFY*/ /*}*/

#define put(x) { *outbuf++ = x; if (outbuf == obe) outbuf = (*fput)(f,&rv); }

 static int
x_sprintf(char *obe, Putfunc fput, Finfo *f, const char *fmt, va_list ap)
{
	char *digits, *ob0, *outbuf, *s, *s0, *se;
	Const char *fmt0;
	char buf[32];
	double x;
	int alt, base, c, decpt, dot, conv, i1, k, lead0, left,
		len, prec, prec1, psign, rv, sgn, sign, width;
#ifdef QUOTIFY
	int quote;
#endif /*QUOTIFY*/
	long sl;
	short sh;
	size_t Ltmp, *ip, j, u;
	ssize_t i;
	unsigned int ui;
	unsigned long ul;
	unsigned short us;
#ifndef NO_PRINTF_A_FMT /*{*/
#ifdef IEEE_8087 /*{{*/
#define I0 1
#define I1 0
#else /*}{*/
#ifdef IEEE_MC68k /*{{*/
#define I0 0
#define I1 1
#else /*}{*/
#define NO_PRINTF_A_FMT
#endif /*}} IEEE_MC68k */
#endif /*}} IEEE_8087  */
#ifndef NO_PRINTF_A_FMT
	typedef union U { double d; unsigned int u[2]; } U;
	U uu;
	int bex, bw;
#endif
#endif /*} NO_PRINTF_A_FMT */
	static char hex[] = "0123456789abcdefpx";
	static char Hex[] = "0123456789ABCDEFPX";
	static char NullStr[] = "<NULL>";

	ob0 = outbuf = f->ob0;
	i = rv = 0;
	for(;;) {
		for(;;) {
			switch(c = *fmt++) {
				case 0:
					goto done;
				case '%':
					break;
				default:
					put(c)
					continue;
				}
			break;
			}
#ifdef QUOTIFY
		quote =
#endif /*QUOTIFY*/
		alt=dot=lead0=left=len=prec=psign=sign=width=0;
		fmt0 = fmt;
 fmtloop:
		switch(conv = *fmt++) {
			case ' ':
			case '+':
				sign = conv;
				goto fmtloop;
			case '-':
				if (dot)
					psign = 1;
				else
					left = 1;
				goto fmtloop;
			case '#':
				alt = 1;
				goto fmtloop;
			case '0':
				if (!lead0 && !dot) {
					lead0 = 1;
					goto fmtloop;
					}
			case '1':
			case '2':
			case '3':
			case '4':
			case '5':
			case '6':
			case '7':
			case '8':
			case '9':
				k = conv - '0';
				while((c = *fmt) >= '0' && c <= '9') {
					k = 10*k + c - '0';
					fmt++;
					}
				if (dot)
					prec = psign ? -k : k;
				else
					width = k;
				goto fmtloop;
			case 'h':
				len = 2;
				goto fmtloop;
			case 'l':
				len = 1;
				goto fmtloop;
			case 'z':
				len = 3;
				goto fmtloop;
			case '.':
				dot = 1;
				goto fmtloop;
			case '*':
				k = va_arg(ap, int);
				if (dot)
					prec = k;
				else {
					if (k < 0) {
						sign = '-';
						k = -k;
						}
					width = k;
					}
				goto fmtloop;
			case 'c':
				c = va_arg(ap, int);
				put(c)
				continue;
			case '%':
				put(conv)
				continue;
			case 'u':
				switch(len) {
				  case 0:
					ui = va_arg(ap, int);
					i = ui;
					break;
				  case 1:
					sl = va_arg(ap, long);
					i = sl;
					break;
				  case 2:
					us = va_arg(ap, int);
					i = us;
					break;
				  case 3:
					i = va_arg(ap, size_t);
				  }
				sign = 0;
				goto have_i;
			case 'i':
			case 'd':
				switch(len) {
				  case 0:
					k = va_arg(ap, int);
					i = k;
					break;
				  case 1:
					sl = va_arg(ap, long);
					i = sl;
					break;
				  case 2:
					sh = va_arg(ap, int);
					i = sh;
					break;
				  case 3:
					i = va_arg(ap, ssize_t);
				  }
				if (i < 0) {
					sign = '-';
					i = -i;
					}
 have_i:
				base = 10;
				u = i;
				digits = hex;
 baseloop:
				if (dot)
					lead0 = 0;
				s = buf;
				if (!u)
					alt = 0;
				do {
					j = ULDIV(u, base);
					*s++ = digits[u - base*j];
					}
					while((u = j));
				prec -= c = s - buf;
				if (alt && conv == 'o' && prec <= 0)
					prec = 1;
				if ((width -= c) > 0) {
					if (prec > 0)
						width -= prec;
					if (sign)
						width--;
					if (alt == 2)
						width--;
					}
				if (left) {
					if (alt == 2)
						put('0') /* for 0x */
					if (sign)
						put(sign)
					while(--prec >= 0)
						put('0')
					do put(*--s)
						while(s > buf);
					while(--width >= 0)
						put(' ')
					continue;
					}
				if (width > 0) {
					if (lead0) {
						if (alt == 2)
							put('0')
						if (sign)
							put(sign)
						while(--width >= 0)
							put('0')
						goto s_loop;
						}
					else
						while(--width >= 0)
							put(' ')
					}
				if (alt == 2)
					put('0')
				if (sign)
					put(sign)
 s_loop:
				while(--prec >= 0)
					put('0')
				do put(*--s)
					while(s > buf);
				continue;
			case 'n':
				ip = va_arg(ap, size_t*);
				if (!ip)
					ip = &Ltmp;
				c = outbuf - ob0 + rv;
				switch(len) {
				  case 0:
					*(int*)ip = c;
					break;
				  case 1:
					*(long*)ip = c;
					break;
				  case 2:
					*(short*)ip = c;
					break;
				  case 3:
					*ip = c;
				  }
				break;
			case 'p':
				len = alt = 1;
				/* no break */
			case 'x':
				digits = hex;
				goto more_x;
			case 'X':
				digits = Hex;
 more_x:
				if (alt) {
					alt = 2;
					sign = conv;
					}
				else
					sign = 0;
				base = 16;
 get_u:
				switch(len) {
				  case 0:
					ui = va_arg(ap, int);
					u = ui;
					break;
				  case 1:
					ul = va_arg(ap, long);
					u = ul;
					break;
				  case 2:
					us = va_arg(ap, int);
					u = us;
					break;
				  case 3:
					u = va_arg(ap, size_t);
				  }
				if (!u)
					sign = alt = 0;
				goto baseloop;
			case 'o':
				base = 8;
				digits = hex;
				goto get_u;
#ifdef QUOTIFY
			case 'Q':
				conv = 'q';
				quote = '\'';
			case 'q':
#endif /*QUOTIFY*/
			case 's':
				s0 = 0;
				s = va_arg(ap, char*);
				if (!s)
					s = NullStr;
				if (prec < 0)
					prec = 0;
 have_s:
				if (dot) {
					for(c = 0; c < prec; c++)
						if (!s[c])
							break;
					prec = c;
					}
				else
					prec = strlen(s);
				width -= prec;
#ifdef QUOTIFY
				if (conv == 'q')
					qkind(s,prec,&quote,&width);
#endif /*QUOTIFY*/
				if (!left)
					while(--width >= 0)
						put(' ')
#ifdef QUOTIFY
				if (quote) {
					put(quote)
					while(--prec >= 0) {
						if ((c = *s++) == quote)
							put(quote)
						else if (c == '\n')
							put('\\')
						put(c)
						}
					put(quote)
					}
				else
#endif /*QUOTIFY*/
				while(--prec >= 0)
					put(*s++)
				while(--width >= 0)
					put(' ')
				if (s0)
					freedtoa(s0);
				continue;
			case 'f':
				if (!dot)
					prec = 6;
				x = va_arg(ap, double);
 infnan:
				s = s0 = dtoa(x, 3, prec, &decpt, &sgn, &se);
				if (decpt == 9999) {
 fmt9999:
					dot = prec = alt = 0;
					if (*s == 'N')
						goto have_s;
					decpt = strlen(s);
					}
 f_fmt:
				if (sgn && (x||sign))
					sign = '-';
				if (prec > 0)
					width -= prec;
				if (width > 0) {
					if (sign)
						--width;
					if (decpt <= 0) {
						--width;
						if (prec > 0)
							--width;
						}
					else {
						if (s == se)
							decpt = 1;
						width -= decpt;
						if (prec > 0 || alt)
							--width;
						}
					}
				if (width > 0 && !left) {
					if (lead0) {
						if (sign)
							put(sign)
						sign = 0;
						do put('0')
							while(--width > 0);
						}
					else do put(' ')
						while(--width > 0);
					}
				if (sign)
					put(sign)
				if (decpt <= 0) {
					put('0')
					if (prec > 0 || alt)
						put('.')
					while(decpt < 0) {
						put('0')
						prec--;
						decpt++;
						}
					}
				else {
					do {
						if ((c = *s))
							s++;
						else
							c = '0';
						put(c)
						}
						while(--decpt > 0);
					if (prec > 0 || alt)
						put('.')
					}
				while(--prec >= 0) {
					if ((c = *s))
						s++;
					else
						c = '0';
					put(c)
					}
				while(--width >= 0)
					put(' ')
				freedtoa(s0);
				continue;
			case 'G':
			case 'g':
				if (!dot)
					prec = 6;
				x = va_arg(ap, double);
				if (prec < 0)
					prec = 0;
				s = s0 = dtoa(x, prec ? 2 : 0, prec, &decpt,
					&sgn, &se);
				if (decpt == 9999)
					goto fmt9999;
				c = se - s;
				prec1 = prec;
				if (!prec) {
					prec = c;
					prec1 = c + (s[1] || alt ? 5 : 4);
					/* %.0g gives 10 rather than 1e1 */
					}
				if (decpt > -4 && decpt <= prec1) {
					if (alt)
						prec -= decpt;
					else
						prec = c - decpt;
					if (prec < 0)
						prec = 0;
					goto f_fmt;
					}
				conv -= 2;
				if (!alt && prec > c)
					prec = c;
				--prec;
				goto e_fmt;
			case 'e':
			case 'E':
				if (!dot)
					prec = 6;
				x = va_arg(ap, double);
				if (prec < 0)
					prec = 0;
				s = s0 = dtoa(x, prec ? 2 : 0, prec+1, &decpt,
					&sgn, &se);
				if (decpt == 9999)
					goto fmt9999;
 e_fmt:
				if (sgn && (x||sign))
					sign = '-';
				if ((width -= prec + 5) > 0) {
					if (sign)
						--width;
					if (prec || alt)
						--width;
					}
				if ((c = --decpt) < 0)
					c = -c;
				while(c >= 100) {
					--width;
					c /= 10;
					}
				if (width > 0 && !left) {
					if (lead0) {
						if (sign)
							put(sign)
						sign = 0;
						do put('0')
							while(--width > 0);
						}
					else do put(' ')
						while(--width > 0);
					}
				if (sign)
					put(sign)
				put(*s++)
				if (prec || alt)
					put('.')
				while(--prec >= 0) {
					if ((c = *s))
						s++;
					else
						c = '0';
					put(c)
					}
				put(conv)
				if (decpt < 0) {
					put('-')
					decpt = -decpt;
					}
				else
					put('+')
				for(c = 2, k = 10; 10*k <= decpt; c++, k *= 10);
				for(;;) {
					i1 = decpt / k;
					put(i1 + '0')
					if (--c <= 0)
						break;
					decpt -= i1*k;
					decpt *= 10;
					}
				while(--width >= 0)
					put(' ')
				freedtoa(s0);
				continue;
#ifndef NO_PRINTF_A_FMT
			case 'a':
				digits = hex;
				goto more_a;
			case 'A':
				digits = Hex;
 more_a:
				uu.d = va_arg(ap, double);
				if ((uu.u[I0] & 0x7ff00000) == 0x7ff00000) {
					x = uu.d;
					goto infnan;
					}
				if (uu.d) {
					c = '1';
					if (uu.u[I0] & 0x80000000) {
						sign = '-';
						uu.u[I0] &= 0x7fffffff;
						}
					bex = (uu.u[I0] >> 20) - 1023;
					uu.u[I0] &= 0xfffff;
					if (bex == -1023) {
						++bex;
						if (uu.u[I0])
							do {
								--bex;
								uu.u[I0] <<= 1;
								if (uu.u[I1] & 0x80000000)
									uu.u[I0] |= 1;
								uu.u[I1] <<= 1;
								} while (uu.u[I0] < 0x100000);
						else {
							while(!(uu.u[I1] & 0x80000000)) {
								--bex;
								uu.u[I1] <<= 1;
								}
							bex -= 21;
							uu.u[I0] = uu.u[I1] >> 11;
							uu.u[I1] <<= 21;
							}
						}
					}
				else {
					c = '0';
					bex = 0;
					}
				if (dot) {
					if (prec > 13)
						prec = 13;
					if (uu.d && prec < 13) {
						uu.u[I0] |= 0x100000;
						if (prec < 5) {
						    ui = 1 << ((5-prec)*4 - 1);
						    if (uu.u[I0] & ui) {
							if (uu.u[I0] & ((ui-1) | (ui << 1))
							 || uu.u[I1]) {
							    uu.u[I0] += ui;
 bex_check:
							    if (uu.u[I0] & 0x200000) {
								++bex;
								uu.u[I0] >>= 1;
								}
							    }
							}
						    }
						else if (prec == 5) {
						    if (uu.u[I1] & 0x80000000) {
 u0_bump:
							++uu.u[I0];
							goto bex_check;
							}
						    }
						else {
						    i1 = (13 - prec) * 4;
						    ui = 1 << (i1 - 1);
						    if (uu.u[I1] & ui
						     && uu.u[I1] & ((ui-1) | (ui << 1))) {
							uu.u[I1] += ui;
							if (!(uu.u[I1] >> i1))
								goto u0_bump;
							}
						    }
						}
					}
				else {
					if ((ui = uu.u[I1]))
						for(prec = 6;
							(ui = (ui << 4) & 0xffffffff);
							++prec);
					else
						for(prec = 0, ui = uu.u[I0] & 0xfffff;
							ui;
							++prec, ui = (ui << 4) & 0xfffff);
					}
				bw = 1;
				if (bex) {
					if ((i1 = bex) < 0)
						i1 = -i1;
					while(i1 >= 10) {
						++bw;
						i1 /= 10;
						}
					}
				if ((sgn = uu.u[I0] & 0x80000000)) {
					uu.u[I0] &= 0x7fffffff;
					if (uu.d||sign)
						sign = '-';
					}
				if ((width -= bw + 5) > 0) {
					if (sign)
						--width;
					if (prec || alt)
						--width;
					}
				if (width > 0 && !left) {
					if (lead0) {
						if (sign) {
							put(sign)
							sign = 0;
							}
						do put('0')
							while(--width > 0);
						}
					else do put(' ')
						while(--width > 0);
					}
				if (sign)
					put(sign)
				put('0')
				put(digits[17])
				put(c)
				if (prec > 0 || alt)
					put('.')
				if (prec > 0) {
					if ((i1 = prec) > 5)
						i1 = 5;
					prec -= i1;
					do {
						put(digits[(uu.u[I0] >> 16) & 0xf])
						uu.u[I0] <<= 4;
						}
						while(--i1 > 0);
					while(prec > 0) {
						--prec;
						put(digits[(uu.u[I1] >> 28) & 0xf])
						uu.u[I1] <<= 4;
						}
					}
				put(digits[16])
				if (bex < 0) {
					put('-')
					bex = -bex;
					}
				else
					put('+')
				for(c = 1; 10*c <= bex; c *= 10);
				for(;;) {
					i1 = bex / c;
					put('0' + i1)
					if (!--bw)
						break;
					bex -= i1 * c;
					bex *= 10;
					}
				continue;
#endif /* NO_PRINTF_A_FMT */
			default:
				put('%')
				while(fmt0 < fmt)
					put(*fmt0++)
				continue;
			}
		}
 done:
	*outbuf = 0;
	return (f->lastlen = outbuf - ob0) + rv;
	}

#define Bsize 4096

 int
Printf(const char *fmt, ...)
{
	va_list ap;
	int rv;
	Finfo f;
	char buf[Bsize];

	va_start(ap, fmt);
	f.u.cf = stdout;
	f.ob0 = buf;
	f.obe1 = buf + Bsize - 1;
#ifdef _windows_
	if (fileno(stdout) == stdout_fileno_ASL) {
		rv = x_sprintf(f.obe1, Wput, &f, fmt, ap);
		mwrite(buf, f.lastlen);
		}
	else
#endif
#ifdef PF_BUF
	if (stdout == stderr_ASL) {
		rv = x_sprintf(f.obe1, pfput, &f, fmt, ap);
		pf_put(buf, f.lastlen);
		}
	else
#endif
		{
		rv = x_sprintf(f.obe1, Fput, &f, fmt, ap);
		fputs(buf, stdout);
		}
	va_end(ap);
	return rv;
	}

 static char *
Sput(Finfo *f, int *rvp)
{
	if (Printf("\nBUG! Sput called!\n", f, rvp))
		/* pass vp, rvp and return 0 to shut diagnostics off */
		exit(250);
	return 0;
	}

 int
Sprintf(char *s, const char *fmt, ...)
{
	va_list ap;
	int rv;
	Finfo f;

	va_start(ap, fmt);
	f.ob0 = s;
	rv = x_sprintf(s, Sput, &f, fmt, ap);
	va_end(ap);
	return rv;
	}

 int
Fprintf(FILE *F, const char *fmt, ...)
{
	va_list ap;
	int rv;
	Finfo f;
	char buf[Bsize];

	va_start(ap, fmt);
	f.u.cf = F;
	f.ob0 = buf;
	f.obe1 = buf + Bsize - 1;
#ifdef MESS
	if (stdout_or_err(F)) {
#ifdef _windows_
		if (fileno(stdout) == stdout_fileno_ASL) {
			rv = x_sprintf(f.obe1, Wput, &f, fmt, ap);
			mwrite(buf, f.lastlen);
			}
		else
#endif
#ifdef PF_BUF
		if (F == stderr_ASL) {
			rv = x_sprintf(f.obe1, pfput, &f, fmt, ap);
			pf_put(buf, f.lastlen);
			}
		else
#endif
			{
			rv = x_sprintf(f.obe1, Fput, &f, fmt, ap);
			fputs(buf, F);
			}
		}
	else
#endif /*MESS*/
		{
#ifdef PF_BUF
		if (F == stderr_ASL) {
			rv = x_sprintf(f.obe1, pfput, &f, fmt, ap);
			pf_put(buf, f.lastlen);
			}
	else
#endif
			{
			rv = x_sprintf(f.obe1, Fput, &f, fmt, ap);
			fputs(buf, F);
			}
		}
	va_end(ap);
	return rv;
	}

 int
Vsprintf(char *s, const char *fmt, va_list ap)
{
	Finfo f;
	return x_sprintf(f.ob0 = s, Sput, &f, fmt, ap);
	}

 int
Vfprintf(FILE *F, const char *fmt, va_list ap)
{
	char buf[Bsize];
	int rv;
	Finfo f;

	f.u.cf = F;
	f.ob0 = buf;
	f.obe1 = buf + Bsize - 1;
#ifdef MESS
	if (stdout_or_err(F)) {
#ifdef _windows_
		if (fileno(stdout) == stdout_fileno_ASL) {
			rv = x_sprintf(f.obe1, Wput, &f, fmt, ap);
			mwrite(buf, f.lastlen);
			}
		else
#endif
#ifdef PF_BUF
		if (F == stderr_ASL) {
			rv = x_sprintf(f.obe1, pfput, &f, fmt, ap);
			pf_put(buf, f.lastlen);
			}
		else
#endif
			{
			rv = x_sprintf(f.obe1, Fput, &f, fmt, ap);
			fputs(buf, F);
			}
		}
	else
#endif /*MESS*/
		{
#ifdef PF_BUF
		if (F == stderr_ASL) {
			rv = x_sprintf(f.obe1, pfput, &f, fmt, ap);
			pf_put(buf, f.lastlen);
			}
		else
#endif
			{
			rv = x_sprintf(f.obe1, Fput, &f, fmt, ap);
			fputs(buf, F);
			}
		}
	va_end(ap);
	return rv;
	}

 void
Perror(const char *s)
{
	if (s && *s)
		Fprintf(Stderr, "%s: ", s);
	Fprintf(Stderr, "%s\n", strerror(errno));
	}

 static char *
Snput(Finfo *f, int *rvp)
{
	char *s, *s0;
	size_t L;

	*rvp += Bsize;
	s0 = f->ob0;
	s = f->u.sf;
	if ((L = f->obe1 - s) > Bsize) {
		L = Bsize;
		goto copy;
		}
	if (L > 0) {
 copy:
		memcpy(s, s0, L);
		f->u.sf = s + L;
		}
	return s0;
	}

 int
Vsnprintf(char *s, size_t n, const char *fmt, va_list ap)
{
	Finfo f;
	char buf[Bsize];
	int L, rv;

	if (n <= 0 || !s) {
		n = 1;
		s = buf;
		}
	f.u.sf = s;
	f.ob0 = buf;
	f.obe1 = s + n - 1;
	rv = x_sprintf(buf + Bsize, Snput, &f, fmt, ap);
	if (f.lastlen > (L = f.obe1 - f.u.sf))
		f.lastlen = L;
	if (f.lastlen > 0) {
		memcpy(f.u.sf, buf, f.lastlen);
		f.u.sf += f.lastlen;
		}
	*f.u.sf = 0;
	return rv;
	}
 int
Snprintf(char *s, size_t n, const char *fmt, ...)
{
	int rv;
	va_list ap;

	va_start(ap, fmt);
	rv = Vsnprintf(s, n, fmt, ap);
	va_end(ap);
	return rv;
	}


#ifdef __cplusplus
}
#endif
