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

/* This is for Sun (and other systems with brain-damaged sprintf that
 * doesn't return the count of characters transmitted).
 * It implements only the relevant subset of sprintf.
 * With ANSI C (and other sensible systems), you can omit Sprintf.o
 * if you add -DSprintf=sprintf to CFLAGS (in the makefile).
 */

#include "stdio.h"
#include "string.h"

#ifdef KR_headers
#ifndef size_t__
#define size_t int
#define size_t__
#endif

#include "varargs.h"

#else

#include "stddef.h"
#include "stdarg.h"
#include "stdlib.h"

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
	printf("bad fmt in Sprintf, starting with \"%s\"\n", fmt);
	exit(1);
	}

#define put(x) *outbuf++ = x

 int
Sprintf
#ifdef KR_headers
	(va_alist)
 va_dcl
#else
	(char *outbuf, const char *fmt, ...)
#endif
{
	char *ob0, *s;
	char buf[32];
	va_list ap;
	long i, j;

#ifdef KR_headers
	char *outbuf, *fmt;
	va_start(ap);
	outbuf = va_arg(ap, char*);
	fmt = va_arg(ap, char*);
#else
	va_start(ap, fmt);
#endif
	ob0 = outbuf;
	for(;;) {
		for(;;) {
			switch(i = *fmt++) {
				case 0:
					goto done;
				case '%':
					break;
				default:
					put(i);
					continue;
				}
			break;
			}
		switch(*fmt++) {
			case 'c':
				i = va_arg(ap, int);
				put(i);
				continue;
			case 'l':
				if (*fmt != 'd')
					bad(fmt);
				fmt++;
				i = va_arg(ap, long);
				goto have_i;
			case 'd':
				i = va_arg(ap, int);
 have_i:
				if (i < 0) {
					put('-');
					i = -i;
					}
				s = buf;
				do {
					j = i / 10;
					*s++ = i - 10*j + '0';
					}
					while(i = j);
				do {
					i = *--s;
					put(i);
					}
					while(s > buf);
				continue;
			case 's':
				s = va_arg(ap, char*);
				while(i = *s++)
					{ put(i); }
				continue;
			default:
				bad(fmt);
			}
		}
 done:
	*outbuf = 0;
	return outbuf - ob0;
	}
#ifdef __cplusplus
	}
#endif
