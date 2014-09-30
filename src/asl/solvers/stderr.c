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

#include <stdio.h>

FILE *Stderr = 0;

#ifdef __cplusplus
extern "C" void Stderr_init_ASL(void);
#endif

#ifdef _WIN32

/* In the MS Windows world, we must jump through */
/* extra hoops in case we're just in a .dll, e.g., */
/* if we're used in a MATLAB mex function. */

#include <windows.h>
#include <io.h>		/* for _open_osfhandle() */
#include <fcntl.h>	/* for _O_TEXT */
#include "arith.h"	/* for LONG_LONG_POINTERS */

#ifdef LONG_LONG_POINTERS
#define Long long long
#else
#define Long long
#endif

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
