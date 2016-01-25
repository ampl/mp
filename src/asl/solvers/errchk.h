/****************************************************************
Copyright (C) 1998 Lucent Technologies
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

#if defined(NO_ERRNO) || (!defined(CHECK_ERRNO) && !defined(ALSO_CHECK_ERRNO))
#define errno_set(x) /*nothing*/
#define ErrnoChk /* nothing */
#else
#define errno_set(x) errno = x
#define ErrnoChk errno ||
#endif

typedef union { real d; unsigned int u[2]; } U;

#ifdef IEEE_MC68k
#define errchk(x) ErrnoChk ((x.u[0] & 0x7ff00000) == 0x7ff00000)
#elif defined(IEEE_8087)
#define errchk(x) ErrnoChk ((x.u[1] & 0x7ff00000) == 0x7ff00000)
#elif !defined(NO_ERRNO)
#define errchk(x) errno
#else
!!! Cannot use -DNO_ERRNO on this machine!!!
#endif
