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

/* As an alternative to "make arith.h", on systems with compilers sufficiently */
/* compatible with gcc or Microsoft C/C++, it may suffice to link or copy this */
/* file to arith.h . */

#ifdef _WIN64
#define X64_bit_pointers
#define LONG_LONG_POINTERS
#define ssize_t long long
#define IEEE_8087
#elif defined(_WIN32)
#define ssize_t long
#define NO_LONG_LONG
#define IEEE_8087
#elif defined(__GNUC__)
#if defined(__ORDER_LITTLE_ENDIAN__)
#define IEEE_8087
#elif defined(__ORDER_BIG_ENDIAN__)
#define IEEE_MC68k
#else
Unexpected byte ordering!!!!
#endif
#endif
#if __SIZEOF_POINTER__ == 8
#define X64_bit_pointers
#endif
#define Double_Align
#ifdef IEEE_8087
#define Arith_Kind_ASL 1
#define QNaN0 0x0
#define QNaN1 0xfff80000
#elif defined(IEEE_MC68k)
#define Arith_Kind_ASL 2
#define QNaN0 0xfff80000
#define QNaN1 0x0
#endif
#ifndef Long
#define Long int
#endif
