/****************************************************************
Copyright (C) 1999 Lucent Technologies
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

#include "asl.h"
#include "signal.h"

#ifdef KR_headers
#define ANSI(x) ()
#else
#define ANSI(x) x
#endif

#ifdef __cplusplus
extern "C" {
#endif

extern void mainexit_ASL ANSI((int));

#ifndef Sig_ret_type
#define Sig_ret_type void
#endif

#ifdef SIGHUP
 static Sig_ret_type
#ifdef KR_headers
hupcatch(n) int n;
#else
hupcatch(int n)
#endif
{
	mainexit_ASL(n);
	}
#endif

 void
sigcatch_ASL ANSI((void))
{
#ifdef SIGHUP
	int i;
	static int sig[] = { SIGABRT, SIGQUIT, SIGTERM, 0 };

	if (signal(SIGHUP, hupcatch) == SIG_IGN)
		signal(SIGHUP, SIG_IGN);
	for(i = 0; sig[i]; i++)
		signal(sig[i], hupcatch);
#endif
	}

#ifndef _WIN32

 extern void (*breakfunc_ASL) ANSI((int,void*)), *breakarg_ASL;

 static Sig_ret_type
#ifdef KR_headers
intcatch(n) int n;
#else
intcatch(int n)
#endif
{
	(*breakfunc_ASL)(n, breakarg_ASL);
	}

 void
#ifdef KR_headers
intcatch_ASL(a, f, v) ASL *a; void (*f)(); void *v;
#else
intcatch_ASL(ASL *a, void (*f)(int,void*), void *v)
#endif
{
	AmplExports *ae;
	if (f) {
		breakfunc_ASL = f;
		breakarg_ASL = v;
		signal(SIGINT, intcatch);
		}
	else
		signal(SIGINT, SIG_IGN);
	if (ae = a->i.ae) {
		ae->Breakfunc = f;
		ae->Breakarg = v;
		}
	}
#endif /*_WIN32*/

#ifdef __cplusplus
}
#endif
