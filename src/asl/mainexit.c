/****************************************************************
Copyright (C) 1997-1999 Lucent Technologies
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

/* Replaceable exit routine */

#ifdef __cplusplus
extern "C" {
#endif

#ifdef KR_headers

extern void at_exit_ASL();

#ifdef PF_BUF
extern void (*pfbuf_print_ASL)();
extern char *pfbuf_ASL;
#endif

 void
mainexit_ASL(n) int n;

#else
#include "stdlib.h"

extern void at_exit_ASL(void);

#ifdef PF_BUF
extern void (*pfbuf_print_ASL)(char*);
extern char *pfbuf_ASL;
#endif

 void
mainexit_ASL(int n)
#endif
{
#ifdef PF_BUF
	if (pfbuf_ASL && pfbuf_print_ASL)
		(*pfbuf_print_ASL)(pfbuf_ASL);
#endif
	at_exit_ASL();
	exit(n);
	}

#ifdef __cplusplus
	}
#endif
