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

#include "asl_pfgh.h"

 int
#ifdef KR_headers
htcl_ASL(x) register unsigned int x;
#else
htcl_ASL(register unsigned int x)
#endif
{
	register int k = 0;
	register unsigned int L = sizeof(char *);

	while(L < x) {
		k++;
		if (!(L <<= 1))
			break;
		}
	return k;
	}

 Char *
#ifdef KR_headers
new_mblk_ASL(a, k) ASL *a; int k;
#else
new_mblk_ASL(ASL *a, int k)
#endif
{
	Char *rv, **t;
	switch(a->i.ASLtype) {
	  case ASL_read_pfg:
	  case ASL_read_pfgh:
		break;
	  default:
		badasl_ASL(a, ASL_read_pfgh, "new_mblk");
	  }
	ACQUIRE_DTOA_LOCK(MBLK_LOCK);
	t = ((ASL_pfgh*)a)->mblk_free + k;
	if (rv = *t) {
		*t = *(Char**)rv;
		FREE_DTOA_LOCK(MBLK_LOCK);
		}
	else {
		FREE_DTOA_LOCK(MBLK_LOCK);
		rv = mem_ASL(a, sizeof(char*)<<k);
		}
	return rv;
	}

#ifndef PSHVREAD
 static
#endif
 void
#ifdef KR_headers
Del_mblk_ASL(a, k, x) ASL *a; int k; Char *x;
#else
Del_mblk_ASL(ASL *a, int k, Char *x)
#endif
{
	Char **t;
	switch(a->i.ASLtype) {
	  case ASL_read_pfg:
	  case ASL_read_pfgh:
		break;
	  default:
		badasl_ASL(a, ASL_read_pfgh, "del_mblk");
	  }
	t = ((ASL_pfgh*)a)->mblk_free + k;
	*(Char**)x = *t;
	*t = x;
	}
