/****************************************************************
Copyright (C) 1992, 1998 Lucent Technologies
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

/* Many Unix systems have a getrusage() that returns time with */
/* a resolution of microseconds.  Compile with -DNO_RUSAGE to */
/* avoid using this routine. */

#ifdef KR_headers
#define Void /*void*/
#else
#define Void void
#endif
#ifdef __cplusplus
extern "C" double xectim_(Void);
#endif

#ifdef __STDC__
#define ASSUME_STDC
#endif
#ifdef _WIN32
#undef  NO_RUSAGE
#define NO_RUSAGE
#undef  ASSUME_STDC
#define ASSUME_STDC
#include <windows.h>
#endif

#ifndef NO_RUSAGE /*{{*/
#include <sys/time.h>
#include <sys/resource.h>

 double
xectim_(Void)
{
	struct rusage R;
	getrusage(RUSAGE_SELF, &R);
	return R.ru_utime.tv_sec + R.ru_stime.tv_sec
		+ 1e-6*(R.ru_utime.tv_usec + R.ru_stime.tv_usec);
	}

#else /* }NO_RUSAGE{ */

#if defined(ASSUME_STDC) && !defined(Clocks_to_seconds) /*{{*/
#include "time.h"

 double
xectim_(Void)
{
#ifdef _WIN32
	FILETIME PrCT, PrET, PrKT, PrUT;
	double st, ut;
	static HANDLE mypid;
	static int first = 1;

	if (first) {
		first = 0;
		mypid = GetCurrentProcess();
		}
	if (GetProcessTimes(mypid, &PrCT, &PrET, &PrKT, &PrUT)) {
		ut = (PrUT.dwLowDateTime
				+ 4294967296.*PrUT.dwHighDateTime)*1e-7;
		st = (PrKT.dwLowDateTime
				+ 4294967296.*PrKT.dwHighDateTime)*1e-7;
		return st + ut;
		}
#endif
	return (double)clock() / CLOCKS_PER_SEC;
	}

#else /*}defined(__STDC__) && !defined(Clocks_to_seconds){*/
#include "sys/types.h"
#include "sys/times.h"
#ifndef Clocks_to_seconds /*{*/
#ifdef CRAY /*{{*/
#define Clocks_to_seconds * 9.5e-9
#else /*}CRAY{*/

#ifndef HZ /*{*/
#include "limits.h"
#ifdef CLOCKS_PER_SEC /*{{*/
#define HZ (double)(CLOCKS_PER_SEC)
#else /*}CLOCKS_PER_SEC{*/
#ifdef CLK_TCK /*{{*/
#define HZ (double)(CLK_TCK)
#else /*}CLK_TCK{*/
#ifdef mips /*{{*/
#define HZ 100.0
#else /*}mips{*/
#define HZ 60.0
#endif /*}}mips*/
#endif /*}}CLK_TCK*/
#endif /*}}CLOCKS_PER_SEC*/
#endif /*}HZ*/
#define Clocks_to_seconds / HZ
#endif /*}}CRAY*/
#endif /*}Clocks_to_seconds*/

double xectim_(Void)
{
	struct tms t;

	times(&t);
	return((t.tms_utime + t.tms_stime) Clocks_to_seconds);
	}
#endif /*}}defined(__STDC__) && !defined(Clocks_to_seconds)*/
#endif /*}}NO_RUSAGE*/
