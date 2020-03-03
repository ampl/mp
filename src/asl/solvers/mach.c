/*******************************************************************
Copyright (C) 2019 AMPL Optimization, Inc.; written by David M. Gay.

Permission to use, copy, modify, and distribute this software and its
documentation for any purpose and without fee is hereby granted,
provided that the above copyright notice appear in all copies and that
both that the copyright notice and this permission notice and warranty
disclaimer appear in supporting documentation.

The author and AMPL Optimization, Inc. disclaim all warranties with
regard to this software, including all implied warranties of
merchantability and fitness.  In no event shall the author be liable
for any special, indirect or consequential damages or any damages
whatsoever resulting from loss of use, data or profits, whether in an
action of contract, negligence or other tortious action, arising out
of or in connection with the use or performance of this software.
*******************************************************************/

/* machine-dependent initializations */

#include "arith.h"
#include "math.h"

#ifndef Long
#define Long int
#endif

double Infinity, negInfinity;
#ifdef __cplusplus
extern "C" void Mach_ASL();
#define Cextern "C"
#else
#define Cextern /*nothing*/
#endif
#ifdef KR_headers
#define Void /*void*/
#else
#define Void void
#endif

#if !defined(NO_arith_h_check) && (defined(IEEE_8087) || defined(IEEE_MC68k))
#include "stdio1.h"
#ifndef Stderr
extern FILE *Stderr;
#endif
#include "stdlib.h" /* for exit() */

extern Cextern int Same_Double_ASL(double, double);

 static int
NaN_test(Void)
{
	union { unsigned int u[2]; double d; } U;
	U.u[0] = QNaN0;
	U.u[1] = QNaN1;
	return Same_Double_ASL(U.d, U.d);
	}
#endif

 void
Mach_ASL(Void)
{
#ifdef IEEE_8087
	union {
		double d;
		Long L[2];
		} u;
	u.L[1] = 0x7ff00000;
	u.L[0] = 0;
	Infinity = u.d;
/* would use #elif defined, but MIPS doesn't understand #elif */
#else
#ifdef IEEE_MC68k
	union {
		double d;
		Long L[2];
		} u;
	u.L[0] = 0x7ff00000;
	u.L[1] = 0;
	Infinity = u.d;
#else
#ifdef HUGE_VAL
	Infinity = HUGE_VAL;
#else
#ifdef HUGE
#ifdef CRAY
	Infinity = 0.5 * HUGE;
	/* brain-damaged machine dies testing if (HUGE > -HUGE) */
#else
	Infinity = HUGE;
#endif
#else
	Cannot define Infinity!!!
#endif
#endif
#endif
#endif
#if !defined(NO_arith_h_check) && (defined(IEEE_8087) || defined(IEEE_MC68k))
	if (Infinity == 0. || !Same_Double_ASL(Infinity, Infinity + Infinity) || NaN_test()) {
		fprintf(Stderr, "Compiled with an incorrect \"arith.h\".\n"
			"Please invoke \"make clean; make\" in your solvers directory,\n"
			"then relink this program.\n");
		exit(1);
		}
#endif
	negInfinity = -Infinity;
	}
