/****************************************************************
Copyright (C) 1993 Lucent Technologies
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


/* Activate one of the following four #defines, as approproate; */

/*#define IEEE_MC68k*/	/* for SPARC, 3B2, Aliant, Motorola and other
			   big-endian IEEE arithmetic machines. */

/*#define IEEE_8087*/	/* for Intel 80x87 arithmetic processors
			   (e.g., WG6386 machines) and other
			   little-endian IEEE arithmetic machines. */

/*#define VAX*/		/* VAX arithmetic */

/*#define IBM*/		/* IBM-mainframe arithmetic (see below for UTS) */

/*#define CRAY*/	/* Cray arithmetic */

#ifdef CRAY
#define No_dtoa
#endif

/* Note that MIPS chips can be either big-endian or little-endian;
 * for example, MIPS and SGI computers are big-endian (#define IEEE_MC68k)
 * whereas DECStations are little-endian (#define IEEE_8087).
 */

/* If your C >> operator does a logical rather than an arithmetic right
 * shift, you need to #define Unsigned_Shifts .  This affects the
 * conversion routines in dtoa.c (included by dtoa1.c).  In an AMPL
 * session, if
 *		print .1;
 * prints something other than .1, try #define Unsigned_Shifts .
 * The following #ifdef UTS is appropriate for the UTS on mhuxo...
 */

#ifdef UTS
#define IBM
#define Unsigned_Shifts
#endif
