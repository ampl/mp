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
	int ij, j1, kv, *vmi;
	ograd *gr, **gr0;
	real f, *vscale;

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
	co_index = -(i + 1);
	x2_check(X);
	if (!asl->i.noxval)
		asl->i.noxval = (int*)M1zapalloc(n_obj*sizeof(int));
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
	kv = 0;
	vmi = 0;
	if ((vscale = asl->i.vscale))
		kv = 2;
	if (asl->i.vmap) {
		vmi = get_vminv_ASL(a);
		++kv;
		}
	gr = *gr0;
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
	if (!asl->i.x_known) {
		co_index = -(i + 1);
		x2_check(X);
		}
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
