/*******************************************************************
Copyright (C) 2017, 2018 AMPL Optimization, Inc.; written by David M. Gay.

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

 void
con2val_ASL(ASL *a, real *X, real *F, fint *nerror)
{
	ASL_fgh *asl;
	Jmp_buf err_jmp0;
	cde *d, *d0;
	cgrad *gr, **gr0;
	expr *e1;
	int *cm, i, j, j1, k, kv, *vmi;
	real *cscale, f, *vscale;

	ASL_CHECK(a, ASL_read_fgh, "con2val");
	asl = (ASL_fgh*)a;
	if (nerror && *nerror >= 0) {
		err_jmp = &err_jmp0;
		i = setjmp(err_jmp0.jb);
		if ((*nerror = i))
			goto done;
		}
	want_deriv = want_derivs;
	errno = 0;	/* in case f77 set errno opening files */
	j = n_conjac[0];
	if (!asl->i.x_known) {
		co_index = j;
		x2_check_ASL(asl,X);
		}
	if (!(x0kind & ASL_have_concom)) {
		if (comb < combc)
			com2eval_ASL(asl, comb, combc);
		if (comc1)
			com21eval_ASL(asl, 0,comc1);
		x0kind |= ASL_have_concom;
		}
	d0 = con_de;
	k = n_conjac[1];
	cscale = asl->i.cscale;
	cm = asl->i.cmap;
	kv = 0;
	vmi = 0;
	if ((vscale = asl->i.vscale))
		kv = 2;
	if (asl->i.vmap) {
		vmi = get_vminv_ASL(a);
		++kv;
		}
	if (!(gr0 = asl->i.Cgrad0))
		asl->i.Cgrad0 = gr0 = asl->i.Cgrad_;
	for(; j < k; ++j) {
		i = j;
		if (cm)
			i = cm[j];
		co_index = i;
		d = d0 + i;
		e1 = d->e;
		f = (*e1->op)(e1 C_ASL);
		if (!F)
			continue;
		gr = gr0[i];
		switch(kv) {
		  case 3:
			for(; gr; gr = gr->next) {
				j1 = vmi[gr->varno];
				f += X[j1] * vscale[j1] * gr->coef;
				}
			break;
		  case 2:
			for(; gr; gr = gr->next) {
				j1 = gr->varno;
				f += X[j1] * vscale[j1] * gr->coef;
				}
			break;
		  case 1:
			for(; gr; gr = gr->next)
				f += X[vmi[gr->varno]] * gr->coef;
			break;
		  case 0:
			for(; gr; gr = gr->next)
				f += X[gr->varno] * gr->coef;
		  }
		if (cscale)
			f *= cscale[j];
		*F++ = f;
		}
	x0kind |= ASL_have_conval;
 done:
	err_jmp = 0;
	}

 void
jac2val_ASL(ASL *a, real *X, real *G, fint *nerror)
{
	ASL_fgh *asl;
	Jmp_buf err_jmp0;
	cde *d, *d0;
	cgrad *gr, **gr0, *gr1;
	fint ne0;
	int *cm, i, i1, j, k, xksave, *vmi;
	real *Adjoints, *cscale, t, *vscale;
	size_t L;
	static char who[] = "jac2val";

	ASL_CHECK(a, ASL_read_fgh, who);
	asl = (ASL_fgh *)a;
	if (!want_derivs)
		No_derivs_ASL(who);
	ne0 = -1;
	if (nerror && (ne0 = *nerror) >= 0) {
		err_jmp = &err_jmp0;
		L = setjmp(err_jmp0.jb);
		if ((*nerror = L))
			goto done;
		}
	errno = 0;	/* in case f77 set errno opening files */
	co_index = j = n_conjac[0];
	if ((!asl->i.x_known && x2_check_ASL(asl,X))
	|| !(x0kind & ASL_have_conval)) {
		xksave = asl->i.x_known;
		asl->i.x_known = 1;
		con2val_ASL(a, X, 0, nerror);
		asl->i.x_known = xksave;
		if (ne0 >= 0 && *nerror)
			goto done;
		}
	k = n_conjac[1];
	if (asl->i.Derrs)
		deriv_errchk_ASL(a, nerror, j, k-j);
	Adjoints = adjoints;
	d0 = con_de;
	cscale = asl->i.cscale;
	vscale = asl->i.vscale;
	cm = asl->i.cmap;
	vmi = 0;
	if (asl->i.vmap)
		vmi = get_vminv_ASL(a);
	if (f_b)
		fun2set_ASL(asl, f_b);
	if (f_c)
		fun2set_ASL(asl, f_c);
	gr0 = asl->i.Cgrad0;
	for(; j < k; ++j) {
		i = j;
		if (cm)
			i = cm[j];
		d = d0 + i;
		for(gr = gr1 = gr0[i]; gr; gr = gr->next)
			Adjoints[gr->varno] = gr->coef;
		if ((L = d->zaplen)) {
			memset(adjoints_nv1, 0, L);
			derprop(d->d);
			}
		if (vscale) {
			gr = gr1;
			if (vmi)
				for(; gr; gr = gr->next) {
					i1 = gr->varno;
					Adjoints[i1] *= vscale[vmi[i1]];
					}
			else
				for(; gr; gr = gr->next) {
					i1 = gr->varno;
					Adjoints[i1] *= vscale[i1];
					}
			}
		gr = gr1;
		if (cscale)
			for(t = cscale[j]; gr; gr = gr->next)
				G[gr->goff] = t*Adjoints[gr->varno];
		else
			for(; gr; gr = gr->next)
				G[gr->goff] = Adjoints[gr->varno];
		}
 done:
	err_jmp = 0;
	}

 int
jac2dim_ASL(ASL *asl, const char *stub, fint *M, fint *N, fint *NO, fint *NZ,
	fint *MXROW, fint *MXCOL, ftnlen stub_len)
{
	FILE *nl;

	nl = jac_dim_ASL(asl, stub, M, N, NO, NZ, MXROW, MXCOL, stub_len);
	if (!nl)
		return ASL_readerr_nofile;
	X0 = (real *)M1alloc(n_var*sizeof(real));
	return pfgh_read_ASL(asl, nl, ASL_return_read_err);
	}
