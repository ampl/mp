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

#include "nlp.h"
#include "nlc.h"
#define asl cur_ASL
#include <errno.h>
#include <setjmp.h>
#ifdef __cplusplus
extern "C" {
#endif
 real acosh_ ANSI((real *));
 real asinh_ ANSI((real *));
 real acoshd_ ANSI((real *, real *));
 real asinhd_ ANSI((real *, real *));
 void in_trouble ANSI((char *, real));
 void in_trouble2 ANSI((char *, real, real));
 int zerset ANSI((real *, int));
 int zerset_ ANSI((real *, fint *));
 void domain_ ANSI((char *, real *, fint));
 void zerdiv_ ANSI((real *));
 real feval_ ANSI((fint *, fint *, real *, real *));
 int ceval_ ANSI((fint *, double *, double *, double *));
#ifdef __cplusplus
	}
#endif

 void
#ifdef KR_headers
in_trouble(who, a) char *who; real a;
#else
in_trouble(char *who, real a)
#endif
{
	char buf[64];
	if (err_jmp)
		longjmp(err_jmp->jb, errno ? errno : 1);
	sprintf(buf, "can't evaluate %s(%g)", who, a);
	perror(buf);
	exit(1);
	}

 void
#ifdef KR_headers
in_trouble2(who, a, b) char *who; real a, b;
#else
in_trouble2(char *who, real a, real b)
#endif
{
	char buf[64];
	if (err_jmp)
		longjmp(err_jmp->jb, errno ? errno : 1);
	sprintf(buf, "can't evaluate %s(%g,%g)", who, a, b);
	perror(buf);
	exit(1);
	}

 void
#ifdef KR_headers
domain_(who, x, wholen) char *who; real *x; fint wholen;
#else
domain_(char *who, real *x, fint wholen)
#endif
{
	if (err_jmp)
		longjmp(err_jmp->jb, errno = EDOM);
	fprintf(Stderr, "can't compute %.*s(%g)\n", wholen, who, *x);
	exit(1);
	}

 void
#ifdef KR_headers
zerdiv_(L) real *L;
#else
zerdiv_(real *L)
#endif
{
	if (err_jmp)
		longjmp(err_jmp->jb, errno = EDOM);
	fprintf(Stderr, "can't compute %g/0\n", *L);
	exit(1);
	}

 real
#ifdef KR_headers
 asinh_(x) real *x;
#else
 asinh_(real *x)
#endif
{
	real rv, t;
	int sign;

	if (sign = (t = *x) < 0.)
		t = -t;
	rv = log(t + sqrt(t*t + 1.));
	if (sign)
		rv = -rv;
	if (errno)
		in_trouble("asinh",t);
	return rv;
	}

 real
#ifdef KR_headers
 asinhd_(x, dL) real *x, *dL;
#else
 asinhd_(real *x, real *dL)
#endif
{
	real rv, t, t1;
	int sign;

	if (sign = (t = *x) < 0.)
		t = -t;
	rv = log(t + (t1 = sqrt(t*t + 1.)));
	if (sign)
		rv = -rv;
	if (errno)
		in_trouble("asinh",t);
	*dL = 1. / t1;
	return rv;
	}

 real
#ifdef KR_headers
acosh_(x) real *x;
#else
acosh_(real *x)
#endif
{
	real rv, t;

	if ((t = *x) < 1.) {
		errno = EDOM;
		rv = 0.;
		}
	else
		rv = log(t + sqrt(t*t - 1.));
	if (errno)
		in_trouble("acosh",t);
	return rv;
	}

 real
#ifdef KR_headers
acoshd_(x, dL) real *x, *dL;
#else
acoshd_(real *x, real *dL)
#endif
{
	real rv, t, t1;

	if ((t = *x) < 1.) {
		errno = EDOM;
		rv = 0.;
		}
	else
		rv = log(t + (t1 = sqrt(t*t - 1.)));
	if (errno)
		in_trouble("acosh",t);
	*dL = 1. / t1;
	return rv;
	}
