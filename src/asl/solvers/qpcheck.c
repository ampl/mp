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

#include "nlp.h"

 fint
#ifdef KR_headers
qpcheck_ASL(a, rowqp, colqp, delsqp)
	ASL *a; fint **rowqp,**colqp;double **delsqp;
#else
qpcheck_ASL(ASL *a, fint **rowqp, fint **colqp, double **delsqp)
#endif
{

	fint rv;
	int i;
	ASL_fg *asl;

	ASL_CHECK(a, ASL_read_fg, "qpcheck");
	asl = (ASL_fg*)a;
	i = obj_no;
	if (i < 0 || i >= n_obj)
		return 0;
	rv = nqpcheck(i, rowqp, colqp, delsqp);
	if (rv < 0) {
		if (rv == -2L)
			fprintf(Stderr,
			 "Quadratic objective involves division by 0.\n");
		else
			fprintf(Stderr,
			 "Sorry, %s can't handle nonlinearities.\n",
				progname ? progname : "");
		exit(1);
		}
	return rv;
	}
