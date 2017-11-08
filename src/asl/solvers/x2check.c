/*******************************************************************
Copyright (C) 2017 AMPL Optimization, Inc.; written by David M. Gay.

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

#include "jac2dim.h"

 int
x2_check_ASL(ASL_fgh *asl, real *X)
{
	expr_v *V;
	int *vm;
	real *Xe, *vscale;

	if (x0len == 0) {
		x0kind = 0;
		return 0;
		}
	if (x0kind == ASL_first_x || memcmp(Lastx, X, x0len)) {
		if (asl->i.Derrs)
			deriv_errclear_ASL(&asl->i);
		want_deriv = want_derivs;
		memcpy(Lastx, X, x0len);
		asl->i.nxval++;
		V = var_e;
		Xe = (real*)((char*)X + x0len);
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

 int
x2known_ASL(ASL *asl, real *X, fint *nerror)
{
	Jmp_buf err_jmp0;
	int ij, rc;

	ASL_CHECK(asl, ASL_read_fgh, "x2known");
	rc = 1;
	if (asl->i.xknown_ignore)
		goto ret;
	if (nerror && *nerror >= 0) {
		asl->i.err_jmp_ = &err_jmp0;
		ij = setjmp(err_jmp0.jb);
		if ((*nerror = ij))
			goto done;
		}
	errno = 0;	/* in case f77 set errno opening files */
	co_index = nlo ? -1 : 0;
	rc = x2_check_ASL((ASL_fgh*)asl, X);
	asl->i.x_known = 1;
 done:
	asl->i.err_jmp_ = 0;
 ret:
	return rc;
	}
