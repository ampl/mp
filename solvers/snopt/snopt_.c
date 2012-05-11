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

/* Variant of snopt.c for use with Fortran compilers that do not */
/* add a trailing underscore to external names, but otherwise */
/* follow the f2c calling conventions... */

/* The following #defines are for use with both SNOPT 6.* and SNOPT 7.* */

#define s1file_	s1file
#define s3opt_	s3opt
#define s8mem_	s8mem
#define sninit_	sninit
#define snmemb_ snmemb
#define snopt_	snopt
#define snset_	snset
#define sntitl_	sntitl
#define sqopt_	sqopt

#include "snopt.c"
