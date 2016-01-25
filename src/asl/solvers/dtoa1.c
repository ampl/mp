/****************************************************************
Copyright (C) 1997, 2001 Lucent Technologies
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

#ifndef No_dtoa /*{{*/
#ifdef __cplusplus
#include "memory.h"
#endif
#include "arith.h"
#include "math.h"
#undef  __MATH_H__
#define __MATH_H__

#include "stdlib.h"	/* A recent (2008) version of cygwin gratuitously
			   #includes stdlib.h.  To prevent confusion, any
			   such #include must come before the next line. */
#define strtod strtod_ASL
#define No_leftright
#ifdef ALLOW_OPENMP /*{*/
#define MULTIPLE_THREADS
#include <omp.h>

/* This assumes use of separate ASL* values for each thread, so we */
/* can ignore ACQUIRE_DTOA_LOCK(n) and FREE_DTOA_LOCK(n) for n >= 2. */

#define dtoa_get_threadno omp_get_thread_num
 static omp_lock_t Locks[2];

 void
init_dtoa_locks(void)
{
	int i;
	static int need_init = 1;

	if (need_init) {
		for(i = 0; i < 2; ++i)
			omp_init_lock(&Locks[i]);
		need_init = 0;
		}
	}

 void
ACQUIRE_DTOA_LOCK(unsigned int n)
{
	if (n < 2)
		omp_set_lock(&Locks[n]);
	}

 void
FREE_DTOA_LOCK(unsigned int n)
{
	if (n < 2)
		omp_unset_lock(&Locks[n]);
	}
#endif /*}*/

#ifndef MALLOC /*{{*/
#define MALLOC dtoamalloc
#include "dtoa.c" /* >= 19991215 */

#include "stdio1.h"

#ifdef __cplusplus
extern "C" {
#endif
#ifndef Stderr
extern FILE *Stderr;
#endif

 void *
dtoamalloc(size_t x)
{
	void *rv = malloc(x);
	if (!rv) {
		fprintf(Stderr, "\nmalloc failure in dtoa!\n");
		exit(63);
		}
	return rv;
	}
#ifdef __cplusplus
	}
#endif
#else /*}{MALLOC*/
#include "dtoa.c"
#endif /*}}MALLOC*/

#else /*}{*/

#ifdef __cplusplus
extern "C" {

extern char *ecvt(double value, int ndigit, int *decpt, int *sign);
extern char *fcvt(double value, int ndigit, int *decpt, int *sign);

 char *
dtoa(double d, int mode, int ndigits, int *decpt, int *sign, char **rve)
#endif
{
	char *s, *s1;

	if (!d) {
		*decpt = 1;
		s = "0";
		if (rve)
			*rve = s + 1;
		return s;
		}
	if (mode > 3)
		mode = 2 + (mode & 1);
	if (mode == 3)
		s = fcvt(d, ndigits, decpt, sign);
	else {
		if (mode <= 1)
			ndigits = 15;
		s = ecvt(d, ndigits, decpt, sign);
		}
	if (*(s1 = s)) {
		while(*++s1);
		while(--s1 > s && *s1 == '0');
		*++s1 = 0;
		}
	if (rve)
		*rve = s1;
	return s;
	}
#endif
#ifdef __cplusplus
}
#endif /*}}*/
/* 20070913:  dtoa.c INFNAN_CHECK := default */
/* 20151020:  sync with /netlib/fp */
