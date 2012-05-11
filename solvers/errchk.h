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

#ifdef NO_ERRNO
#define errno_set(x) /*nothing*/
#define ErrnoChk /* nothing */
#else
#define errno_set(x) errno = x
#define ErrnoChk errno ||
#endif

#ifdef NANCHECK
#ifdef IEEE_MC68k
#define errchk(x) ErrnoChk ((((Long *)&(x))[0] & 0x7ff00000L) == 0x7ff00000L)\
	&& (((Long *)&(x))[1] || ((Long *)&(x))[0] & 0xfffffL)
#else
#ifdef IEEE_8087
#define errchk(x) ErrnoChk ((((Long *)&(x))[1] & 0x7ff00000L) == 0x7ff00000L)\
	&& (((Long *)&(x))[0] || ((Long *)&(x))[1] & 0xfffffL)
#else
!!!! Cannot use -DNANCHECK on non-IEEE machines
#endif
#endif
#else
#define errchk(x) errno
#endif

#ifdef KR_headers
extern char *strerror();
#endif
