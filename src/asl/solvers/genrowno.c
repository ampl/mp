/****************************************************************
Copyright (C) 2000, 2001 Lucent Technologies
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

/* Solvers that handle sparse nonlinear constraints and wish */
/* to have Jacobian nonzeros stored sparsely by columns in */
/* A_rownos may simply call gen_rownos() after calling the .nl */
/* reader (with A_vals left at its default value = NULL). */

 void
#ifdef KR_headers
gen_rownos_ASL(asl) ASL *asl;
#else
gen_rownos_ASL(ASL *asl)
#endif
{
	cgrad *cg, **cgp, **cgpe;
	int *a, n;

	if (n_con <= 0 || nzc <= 0)
		return;
	if (!(a = A_rownos))
		A_rownos = a = (int*)M1alloc(nzc*sizeof(int));
	n = Fortran;
	cgp = Cgrad;
	for(cgpe = cgp + n_con; cgp < cgpe; n++)
		for(cg = *cgp++; cg; cg = cg->next)
			a[cg->goff] = n;
	}
