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

/* This is for compilers (like MetaWare High C) whose sscanf is broken.
 * It implements only the relevant subset of sscanf.
 * With sensible compilers, you can omit sscanf.o
 * if you add -DSscanf=sscanf to CFLAGS (in the makefile).
 */

#ifdef KR_headers
#include "varargs.h"
#else
#include "stddef.h"
#include "stdarg.h"
#include "stdlib.h"
#endif

#include "arith.h" /* for NO_SSIZE_T */
#include "stdio1.h"
#include "string.h"

#ifndef Stderr
extern FILE *Stderr;
#endif

#ifdef KR_headers
#ifndef size_t__
#define size_t int
#define size_t__
#endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

 static void
bad
#ifdef KR_headers
(fmt) char *fmt;
#else
(const char *fmt)
#endif
{
	fprintf(Stderr, "bad fmt in Sscanf, starting with \"%s\"\n", fmt);
	exit(1);
	}

 int
Sscanf
#ifdef KR_headers
	(va_alist)
 va_dcl
#else
	(char *s, const char *fmt, ...)
#endif
{
	char *s0;
	va_list ap;
	long L, *Lp;
	int c, i, *ip, rc, sgn;
	ssize_t Ll, *Llp;

#ifdef KR_headers
	char *fmt, *s;
	va_start(ap);
	s = va_arg(ap, char*);
	fmt = va_arg(ap, char*);
#else
	va_start(ap, fmt);
#endif
	rc = 0;
	for(;;) {
		for(;;) {
			switch(i = *(unsigned char *)fmt++) {
				case 0:
					goto done;
				case '%':
					break;
				default:
					if (i <= ' ') {
						while(*s <= ' ')
							if (!*s++)
								return rc;
						}
					else if (*s++ != i)
						return rc;
					continue;
				}
			break;
			}
		switch(*fmt++) {
			case 'l':
				if (*fmt != 'd')
					bad(fmt);
				fmt++;
				Lp = va_arg(ap, long*);
				L = strtol(s0 = s, &s, 10);
				if (s > s0) {
					rc++;
					*Lp = L;
					continue;
					}
				return rc;
			case 'd':
				ip = va_arg(ap, int*);
				L = strtol(s0 = s, &s, 10);
				if (s > s0) {
					rc++;
					*ip = (int)L;
					continue;
					}
				return rc;
			case 'D':
				Llp = va_arg(ap, ssize_t*);
				sgn = 0;
				c = *s;
				if (c == '-') {
					sgn = 1;
					c = *++s;
					}
				if (c < '0' || c > '9')
					return rc;
				Ll = c - '0';
				while((c = *++s) >= '0' && c <= '9')
					Ll = 10*Ll + c - '0';
				++rc;
				if (sgn)
					Ll = -Ll;
				*Llp = Ll;
				continue;
			default:
				bad(fmt);
			}
		}
 done:
	return rc;
	}
#ifdef __cplusplus
	}
#endif
