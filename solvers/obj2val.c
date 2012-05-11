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

#include "jac2dim.h"

 extern int x2_check_ASL ANSI((ASL_fgh*, real *));
#define x2_check(x) x2_check_ASL(asl,x)

 static void
NNOBJ_chk(ASL *asl, int i, const char *who)
{
	ASL_CHECK(asl, ASL_read_fgh, who);
	if (i < 0 || i >= n_obj) {
		fprintf(Stderr,
			"objval: got NOBJ = %d; expected 0 <= NOBJ < %d\n",
			i, n_obj);
		exit(1);
		}
	}

 real
obj2val_ASL(ASL *a, int i, real *X, fint *nerror)
{
	ASL_fgh *asl;
	Jmp_buf err_jmp0;
	cde *d;
	expr *e1;
	expr_v *V;
	int ij;
	ograd *gr, **gr0;
	real f;

	NNOBJ_chk(a, i, "obj2val");
	asl = (ASL_fgh*)a;
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
	x2_check(X);
	if (!asl->i.noxval)
		asl->i.noxval = (int*)M1zapalloc(n_obj*sizeof(int));
	co_index = -(i + 1);
	if (!(x0kind & ASL_have_objcom)) {
		if (ncom0 > combc)
			comeval(asl, combc, ncom0);
		x0kind |= ASL_have_objcom;
		}
	d = obj_de + i;
	if (d->n_com1)
		com1eval(asl, d->com11, d->n_com1);
	gr0 = Ograd + i;
	e1 = d->e;
	f = (*e1->op)(e1 C_ASL);
	asl->i.noxval[i] = asl->i.nxval;
	gr = *gr0;
	if (asl->i.vscale || asl->i.vmap)
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
obj2grd_ASL(ASL *a, int i, real *X, real *G, fint *nerror)
{
	ASL_fgh *asl;
	Jmp_buf err_jmp0;
	cde *d;
	fint ne0;
	int ij, j, *vmi, xksave, *z;
	ograd *gr, **gr0;
	real *Adjoints, *vscale;
	size_t L;
	static char who[] = "obj2grd";

	NNOBJ_chk(a, i, who);
	asl = (ASL_fgh*)a;
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
		x2_check(X);
	if (!asl->i.noxval || asl->i.noxval[i] != asl->i.nxval) {
		xksave = asl->i.x_known;
		asl->i.x_known = 1;
		obj2val_ASL(a, i, X, nerror);
		asl->i.x_known = xksave;
		if (ne0 >= 0 && *nerror)
			goto done;
		}
	if (asl->i.Derrs)
		deriv_errchk_ASL(a, nerror, -(i+1), 1);
	if (f_b)
		funnelset(asl, f_b);
	if (f_o)
		funnelset(asl, f_o);
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
	vmi = 0;
	if (a->i.vmap)
		vmi = get_vminv_ASL(a);
	gr = *gr0;
	if ((vscale = asl->i.vscale)) {
		if (vmi)
			for(; gr; gr = gr->next) {
				j = vmi[i = gr->varno];
				G[j] = vscale[j] * Adjoints[i];
				}
		else
			for(; gr; gr = gr->next) {
				i = gr->varno;
				G[i] = vscale[i] * Adjoints[i];
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
