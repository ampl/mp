/****************************************************************
Copyright (C) 1997-2001 Lucent Technologies
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

 int
x0_check_ASL(ASL_fg *asl, real *X)
{
	expr_v *V;
	int *vm;
	real *vscale, *Xe;

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
			comeval_ASL(asl, 0, comb);
		return 1;
		}
	return 0;
	}

 void
x1known_ASL(ASL *asl, real *X, fint *nerror)
{
	Jmp_buf err_jmp0;
	int ij;

	ASL_CHECK(asl, ASL_read_fg, "x1known");
	if (asl->i.xknown_ignore)
		return;
	if (nerror && *nerror >= 0) {
		err_jmp = &err_jmp0;
		ij = setjmp(err_jmp0.jb);
		if ((*nerror = ij))
			goto done;
		}
	errno = 0;	/* in case f77 set errno opening files */
	x0_check_ASL((ASL_fg*)asl, X);
	asl->i.x_known = 1;
 done:
	err_jmp = 0;
	}

 static void
NNOBJ_chk(ASL *asl, int i, const char *who)
{
	ASL_CHECK(asl, ASL_read_fg, who);
	if (i < 0 || i >= n_obj) {
		fprintf(Stderr,
		"objval: got NOBJ = %d; expected 0 <= NOBJ < %d\n",
			i, n_obj);
		exit(1);
		}
	}

 real
obj1val_ASL(ASL *a, int i, real *X, fint *nerror)
{
	ASL_fg *asl;
	Jmp_buf err_jmp0;
	cde *d;
	expr *e1;
	expr_v *V;
	int ij;
	ograd *gr;
	real f;

	NNOBJ_chk(a, i, "obj1val");
	asl = (ASL_fg*)a;
	if (nerror && *nerror >= 0) {
		err_jmp = &err_jmp0;
		ij = setjmp(err_jmp0.jb);
		if ((*nerror = ij)) {
			f = 0.;
			goto done;
			}
		}
	want_deriv = want_derivs;
	errno = 0;	/* in case f77 set errno opening files */
	if (!asl->i.x_known)
		x0_check_ASL(asl,X);
	if (!asl->i.noxval)
		asl->i.noxval = (int*)M1zapalloc(n_obj*sizeof(int));
	co_index = -(i + 1);
	if (!(x0kind & ASL_have_objcom)) {
		if (ncom0 > combc)
			comeval_ASL(asl, combc, ncom0);
		if (comc1 < ncom1)
			com1eval_ASL(asl, comc1, ncom1);
		x0kind |= ASL_have_objcom;
		}
	d = obj_de + i;
	gr = Ograd[i];
	e1 = d->e;
	f = (*e1->op)(e1 C_ASL);
	asl->i.noxval[i] = asl->i.nxval;
	if (asl->i.vmap || asl->i.vscale)
		for(V = var_e; gr; gr = gr->next)
			f += gr->coef * V[gr->varno].v;
	else
		for(; gr; gr = gr->next)
			f += gr->coef * X[gr->varno];
 done:
	err_jmp = 0;
	return f;
	}

 void
obj1grd_ASL(ASL *a, int i, real *X, real *G, fint *nerror)
{
	ASL_fg *asl;
	Jmp_buf err_jmp0;
	cde *d;
	fint ne0;
	int ij, j, *vmi, xksave, *z;
	ograd *gr, **gr0;
	real *Adjoints, *vscale;
	size_t L;
	static char who[] = "obj1grd";

	NNOBJ_chk(a, i, who);
	asl = (ASL_fg*)a;
	if (!want_derivs)
		No_derivs_ASL(who);
	ne0 = -1;
	if (nerror && (ne0 = *nerror) >= 0) {
		err_jmp = &err_jmp0;
		ij = setjmp(err_jmp0.jb);
		if ((*nerror = ij))
			goto done;
		}
	errno = 0;	/* in case f77 set errno opening files */
	if (!asl->i.x_known)
		x0_check_ASL(asl,X);
	if (!asl->i.noxval || asl->i.noxval[i] != asl->i.nxval) {
		xksave = asl->i.x_known;
		asl->i.x_known = 1;
		obj1val_ASL(a, i, X, nerror);
		asl->i.x_known = xksave;
		if (ne0 >= 0 && *nerror)
			goto done;
		}
	if (asl->i.Derrs)
		deriv_errchk_ASL(a, nerror, -(i+1), 1);
	if (f_b)
		funnelset_ASL(asl, f_b);
	if (f_o)
		funnelset_ASL(asl, f_o);
	Adjoints = adjoints;
	d = obj_de + i;
	gr0 = Ograd + i;
	for(gr = *gr0; gr; gr = gr->next)
		Adjoints[gr->varno] = gr->coef;
	if ((L = d->zaplen)) {
		memset(adjoints_nv1, 0, L);
		derprop(d->d);
		}
	if (zerograds) {	/* sparse gradients */
		z = zerograds[i];
		while((i = *z++) >= 0)
			G[i] = 0;
		}
	gr = *gr0;
	vmi = 0;
	if (asl->i.vmap)
		vmi = get_vminv_ASL(a);
	if ((vscale = asl->i.vscale)) {
		if (vmi)
			for(; gr; gr = gr->next) {
				j = vmi[i = gr->varno];
				G[j] = Adjoints[i] * vscale[j];
				}
		else
			for(; gr; gr = gr->next) {
				i = gr->varno;
				G[i] = Adjoints[i] * vscale[i];
				}
		}
	else if (vmi)
		for(; gr; gr = gr->next) {
			i = gr->varno;
			G[vmi[i]] = Adjoints[i];
			}
	else
		for(; gr; gr = gr->next) {
			i = gr->varno;
			G[i] = Adjoints[i];
			}
 done:
	err_jmp = 0;
	}
