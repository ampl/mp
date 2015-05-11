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

int isatty_ASL; /* for use with "sw" under NT */

#ifndef _MT
#define _MT 1
#endif

#include <errno.h>
#include <float.h>
#include <math.h>
#include <signal.h>
#include <windows.h>
#include <process.h>
#include "asl.h"

#undef Need_set_errno
#undef errno_defined
#ifdef errno
#ifndef Keep_errno
#define errno_defined
#endif
#endif /*errno*/
#ifdef errno_defined
#ifndef No_set_errno_ASL
#define Need_set_errno
#undef errno
int errno;
static int set_errno(int);
#endif /*No_set_errno_ASL*/
#else
#define set_errno(x) x
#endif /*errno_defined*/

#ifdef oldWATCOM
#define matherr_rettype double
#else
#ifndef matherr_rettype
#define matherr_rettype int
#endif
#endif

#ifndef _DOMAIN
#define _DOMAIN DOMAIN
#endif
#ifndef _SING
#define _SING SING
#endif
#ifndef _TLOSS
#define _TLOSS TLOSS
#endif
#ifndef _OVERFLOW
#define _OVERFLOW OVERFLOW
#endif

#ifdef __cplusplus
extern "C" {
 void fpinit_ASL(void);
 void catch_SIGINT_ASL(void (*)(int,void*), void*);
 }
#endif

 extern void (*breakfunc_ASL)(int,void*), *breakarg_ASL;

 void
intcatch_ASL(ASL *a, void (*f)(int,void*), void *v)
{
	AmplExports *ae;
	if (f) {
		breakfunc_ASL = f;
		breakarg_ASL = v;
		}
	else
		signal(SIGINT, SIG_DFL);
	if (ae = a->i.ae) {
		ae->Breakfunc = f;
		ae->Breakarg = v;
		}
	}

 static void
siglistener(void *arg)
{
	DWORD n;
	HANDLE *Sig = (HANDLE*)arg;
	int s[3], sig;
	void(*oldsig)(int);

	for(;;) {
		s[1] = 0;
		if (!ReadFile(Sig[0], s, sizeof(s), &n,0) || n <= 0)
			break;
		if (s[1] > 0 && s[1] != GetCurrentProcessId())
			s[0] = 0;
		/* Must use hard-wired numbers since the brain-dead compiler */
		/* vendors use inconsistent signal numbering... */
		switch(s[0]) {
		  case 2:
			sig = SIGINT;
			goto int_or_break;
		  case 21:
			sig = SIGBREAK;
 int_or_break:
			if ((oldsig = signal(sig,SIG_IGN)) == SIG_IGN) {
				sig = 0;
				break;
				}
			if (breakfunc_ASL) {
				(*breakfunc_ASL)(sig, breakarg_ASL);
				breakfunc_ASL = 0;
				}
			else
				signal(sig, oldsig);
			break;
		  case 15:
			sig = SIGTERM;
			if ((oldsig = signal(sig,SIG_DFL)) != SIG_IGN)
				signal(sig, oldsig);
			break;
		  default:
			sig = 0;
		  }
		if (sig) {
			s[2]++;
			if (s[1] >= 0)
				s[0] = 0;
			}
		WriteFile(Sig[1], s, n, &n, 0);
		if (sig)
			raise(sig);
		Sleep(50);
		}
	}

 static void
siglisten(void)
{
	static HANDLE Sig[2];
	HANDLE h;
	char *s;

#ifdef LONG_LONG_POINTERS
#define STRTOUL strtoull
#else
#define STRTOUL strtoul
#endif
	if (s = getenv("SW_sigpipe")) {
		if (!(Sig[0] = (HANDLE)STRTOUL(s,&s,10))
		 || *s != ','
		 || !(Sig[1] = (HANDLE)STRTOUL(s+1,&s,10)))
			return;
		if (*s == ',')
			isatty_ASL = (int)strtoul(s+1,&s,10);
		h = (HANDLE)_beginthread(siglistener, 0, Sig);
		SetThreadPriority(h, THREAD_PRIORITY_HIGHEST);
		}
	}

 void
fpinit_ASL(void)
{
	static int first = 1;

#ifndef No_Control87 /* for DEC Alpha */
#ifndef MCW_EM
#ifndef _MCW_EM	/* for cygwin with -mno-cygwin */
#define _MCW_EM 0x0008001F
#endif
#define MCW_EM _MCW_EM
#endif
#ifndef PC_53
#ifndef _PC_53
#define _PC_53 0x00010000
#endif
#define PC_53 _PC_53
#endif
#ifndef MCW_PC
#ifndef _MCW_PC
#define _MCW_PC 0x00030000
#endif
#define MCW_PC _MCW_PC
#endif
	_control87(MCW_EM | PC_53, MCW_EM | MCW_PC);
#endif
	if (first) {
		first = 0;
		siglisten();
		}
	}

#ifdef __GNUC__
#if (__GNUC__ >= 4 && __GNUC_MINOR__ >= 5) || __GNUC__ >= 5
#undef NO_matherr
#define NO_matherr
#endif
#endif

#ifndef NO_matherr
#ifdef __MINGW32__
#define matherr _matherr
#endif

 matherr_rettype
matherr( struct _exception *e )
{
	switch(e->type) {
	  case _DOMAIN:
	  case _SING:
		errno = set_errno(EDOM);
		break;
	  case _TLOSS:
	  case _OVERFLOW:
		errno = set_errno(ERANGE);
	  }
	return 0;
	}
#endif /*NO_matherr*/

#ifdef Need_set_errno

 static int
set_errno(int n)
{ return errno = n; }

#endif /*Need_set_errno*/
