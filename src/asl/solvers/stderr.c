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

#include <stdio.h>

FILE *Stderr = 0;

#ifdef __cplusplus
extern "C" void Stderr_init_ASL(void);
extern "C" int Same_Double_ASL(double, double);
#endif

 int
Same_Double_ASL(double a, double b)
/* Located here rather than in mach.c in hopes of */
/* avoiding defeat by compiler optimizations. */
{ return a == b; }

#ifdef _WIN32

/* In the MS Windows world, we must jump through */
/* extra hoops in case we're just in a .dll, e.g., */
/* if we're used in a MATLAB mex function. */

#include <windows.h>
#include <io.h>		/* for _open_osfhandle() */
#include <fcntl.h>	/* for _O_TEXT */
#include "arith.h"	/* for LONG_LONG_POINTERS */

 void
Stderr_init_ASL(void)
{
	HANDLE h;
	int ih;

	AllocConsole();	/* fails if unnecessary */
	h = GetStdHandle(STD_OUTPUT_HANDLE);
	ih = _open_osfhandle((ssize_t)h, _O_TEXT);
	if (ih == -1)
		Stderr = fopen("con","w");
	else
		Stderr = _fdopen(ih, "w");
	}

#else /*!_WIN32*/

#ifndef STDERR
#define STDERR stderr
#endif

#ifdef __cplusplus
extern "C" void Stderr_init_ASL(void);
#endif

 void
Stderr_init_ASL(void)
{ Stderr = STDERR; }

#endif /*WIN32*/
