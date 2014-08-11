/****************************************************************
Copyright (C) 1997, 1999 Lucent Technologies
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

#include "jac2dim.h"

 int
x2_check_ASL(ASL_fgh *asl, real *X)
{
	expr_v *V;
	int *vm;
	real *Xe, *vscale;

	if (x0kind == ASL_first_x || memcmp(Lastx, X, x0len)) {
		if (asl->i.Derrs)
			deriv_errclear_ASL(&asl->i);
		want_deriv = want_derivs;
		memcpy(Lastx, X, x0len);
		asl->i.nxval++;
		V = var_e;
		Xe = X + n_var;
		vscale = asl->i.vscale;
		if ((vm = asl->i.vmap)) {
			if (vscale)
				while(X < Xe)
					V[*vm++].v = *vscale++ * *X++;
			else
				while(X < Xe)
					V[*vm++].v = *X++;
			}
		else {
			if (vscale)
				while(X < Xe)
					(V++)->v = *vscale++ * *X++;
			else
				while(X < Xe)
					(V++)->v = *X++;
			}
		x0kind = 0;
		if (comb)
			comeval(asl, 0, comb);
		return 1;
		}
	return 0;
	}

 void
x2known_ASL(ASL *a, real *X, fint *nerror)
{
	Jmp_buf err_jmp0;
	int ij;

	ASL_CHECK(a, ASL_read_fgh, "x2known");
	if (a->i.xknown_ignore)
		return;
	if (nerror && *nerror >= 0) {
		a->i.err_jmp_ = &err_jmp0;
		ij = setjmp(err_jmp0.jb);
		if ((*nerror = ij))
			goto done;
		}
	errno = 0;	/* in case f77 set errno opening files */
	x2_check_ASL((ASL_fgh*)a, X);
	a->i.x_known = 1;
 done:
	a->i.err_jmp_ = 0;
	}
