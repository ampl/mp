/****************************************************************
Copyright (C) 1997, 1999, 2000 Lucent Technologies
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

 void
con1val_ASL(ASL *a, register real *X, real *F, fint *nerror)
{
	cde *d, *d0;
	expr *e1;
	real *cscale, f;
	int *cm, i, j, k, useV;
	cgrad *gr, **gr0;
	Jmp_buf err_jmp0;
	ASL_fg *asl;
	expr_v *V;

	ASL_CHECK(a, ASL_read_fg, "con1val");
	asl = (ASL_fg *)a;
	if (nerror && *nerror >= 0) {
		err_jmp = &err_jmp0;
		i = setjmp(err_jmp0.jb);
		if ((*nerror = i))
			goto done;
		}
	want_deriv = want_derivs;
	errno = 0;	/* in case f77 set errno opening files */
	if (!asl->i.x_known)
		x0_check_ASL(asl,X);
	if (!(x0kind & ASL_have_concom)) {
		if (comb < combc)
			comeval_ASL(asl, comb, combc);
		if (comc1)
			com1eval_ASL(asl, 0,comc1);
		x0kind |= ASL_have_concom;
		}
	x0kind |= ASL_have_conval;
	d0 = con_de;
	j = n_conjac[0];
	k = n_conjac[1];
	cscale = asl->i.cscale;
	cm = asl->i.cmap;
	V = var_e;
	useV = asl->i.vscale || asl->i.vmap;
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
		gr = gr0[i];
		if (useV)
			for(; gr; gr = gr->next)
				f += V[gr->varno].v * gr->coef;
		else
			for(; gr; gr = gr->next)
				f += X[gr->varno] * gr->coef;
		if (F) {
			if (cscale)
				f *= cscale[j];
			*F++ = f;
			}
		}
 done:
	err_jmp = 0;
	}

 void
jac1val_ASL(ASL *a, real *X, real *G, fint *nerror)
{
	ASL_fg *asl;
	Jmp_buf err_jmp0;
	cde *d, *d0;
	cgrad *gr, **gr0, *gr1;
	fint ne0;
	int *cm, i, i1, j, k, xksave, *vmi;
	real *Adjoints, *cscale, t, *vscale;
	size_t L;
	static char who[] = "jac1val";

	ASL_CHECK(a, ASL_read_fg, who);
	asl = (ASL_fg*)a;
	if (!want_derivs)
		No_derivs_ASL(who);
	ne0 = -1;
	if (nerror && (ne0 = *nerror) >= 0) {
		err_jmp = &err_jmp0;
		j = setjmp(err_jmp0.jb);
		if ((*nerror = j))
			goto done;
		}
	errno = 0;	/* in case f77 set errno opening files */
	if ((!asl->i.x_known && x0_check_ASL(asl,X))
	 || !(x0kind & ASL_have_conval)) {
		xksave = asl->i.x_known;
		asl->i.x_known = 1;
		con1val_ASL(a, X, 0, nerror);
		asl->i.x_known = xksave;
		if (ne0 >= 0 && *nerror)
			goto done;
		}
	j = n_conjac[0];
	k = n_conjac[1];
	if (asl->i.Derrs)
		deriv_errchk_ASL(a, nerror, j, k-j);
	if (asl->i.zap_J)
		memset(G, 0, asl->i.zap_J);
	Adjoints = adjoints;
	d0 = con_de;
	cscale = asl->i.cscale;
	vscale = asl->i.vscale;
	cm = asl->i.cmap;
	vmi = 0;
	if (asl->i.vmap)
		vmi = get_vminv_ASL(a);
	if (f_b)
		funnelset_ASL(asl, f_b);
	if (f_c)
		funnelset_ASL(asl, f_c);
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
