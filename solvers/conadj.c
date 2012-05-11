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

/* Make adjustments for a nonlinear complementarity problem
 * expressed in AMPL via constraints of the form
 *	lb <= variable * expression <= ub
 * Set cv[i] = index of variable in i-th constraint (first var = 0),
 * and arrange for Jacobian computations to omit the "variable *".
 */

#include "nlp.h"
#include "r_opn.hd"

 int
#ifdef KR_headers
conadj_ASL(a, cv, bailout) ASL *a; int *cv, bailout;
#else
conadj_ASL(ASL *a, int *cv, int bailout)
#endif
{
	cde *d, *dend;
	expr *e;
	ASL_fg *asl;

	ASL_CHECK(a, ASL_read_fg, "conadj");
	asl = (ASL_fg*)a;
	for(d = con_de, dend = d + n_con; d < dend; d++) {
		e = d->e;
		if (e->op != f_OPMULT || e->L.e->op != f_OPVARVAL) {
			if (bailout) {
				fprintf(Stderr,
					"Not a complementarity problem!\n");
				exit(1);
				}
			return 1;
			}
		d->d = d->d->next->next;
		d->e = e->R.e;
		*cv++ = ((expr_v *)(e->L.e))->a;
		}
	return 0;
	}
