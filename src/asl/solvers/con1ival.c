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

#include "nlp.h"

#ifdef __cplusplus
 extern "C" {
 extern  real con1ival_nomap_ASL(ASL *a, int i, real *X, fint *nerror);
 extern void con1grd_nomap_ASL(ASL *a, int i, real *X, real *G, fint *nerror);
	}
#endif

 static void
INchk(ASL *asl, const char *who, int i, int ix)
{
	ASL_CHECK(asl, ASL_read_fg, who);
	if (i < 0 || i >= ix) {
		fprintf(Stderr, "%s: got I = %d; expected 0 <= I < %d\n",
			who, i, ix);
		exit(1);
		}
	}

 static real
cival(ASL_fg *asl, int i, real *X, fint *nerror)
{
	Jmp_buf err_jmp0;
	expr *e;
	int ij, nc;
	real f;

	if (nerror && *nerror >= 0) {
		err_jmp = &err_jmp0;
		ij = setjmp(err_jmp0.jb);
		if ((*nerror = ij))
			return 0.;
		}
	want_deriv = want_derivs;
	errno = 0;	/* in case f77 set errno opening files */
	co_index = i;
	if (!asl->i.x_known)
		x0_check_ASL(asl,X);
	if (!asl->i.ncxval)
		asl->i.ncxval = (int*)M1zapalloc(nclcon*sizeof(int));
	if (!(x0kind & ASL_have_concom)) {
		if (comb < combc)
			comeval_ASL(asl, comb, combc);
		if (comc1)
			com1eval_ASL(asl, 0,comc1);
		x0kind |= ASL_have_concom;
		}
	if (i >= (nc = asl->i.n_con0))
		e = lcon_de[i-nc].e;
	else
		e = con_de[i].e;
	f = (*e->op)(e C_ASL);
	asl->i.ncxval[i] = asl->i.nxval;
	err_jmp = 0;
	return f;
	}

 static real
Conival1(ASL_fg *asl, int i, real *X, fint *nerror)
{
	cgrad *gr;
	int j1, kv, *vmi;
	real f, *vscale;

	if (i < asl->i.n_con0)
		f = cival(asl, i, X, nerror);
	else
		f = 0.;
	kv = 0;
	vmi = 0;
	if ((vscale = asl->i.vscale))
		kv = 2;
	if (asl->i.vmap) {
		vmi = get_vminv_ASL((ASL*)asl);
		++kv;
		}
	gr = asl->i.Cgrad0[i];
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
	return f;
	}

 real
con1ival_nomap_ASL(ASL *a, int i, real *X, fint *nerror)
{
	INchk(a, "con1ival_nomap", i, a->i.n_con0);
	return  Conival1((ASL_fg*)a, i, X, nerror);
	}

 real
con1ival(ASL *a, int i, real *X, fint *nerror)
{
	int *cm;

	INchk(a, "con1ival", i, a->i.n_con_);
	if ((cm = a->i.cmap))
		i = cm[i];
	return  Conival1((ASL_fg*)a, i, X, nerror);
	}

 int
lcon1val(ASL *a, int i, real *X, fint *nerror)
{
	real f;

	INchk(a, "lcon1ival", i, a->i.n_lcon_);
	f = cival((ASL_fg*)a, i + a->i.n_con0, X, nerror);
	return f != 0.;
	}

 static void
Congrd1(ASL_fg *asl, int i, real *X, real *G, fint *nerror)
{
	Jmp_buf err_jmp0;
	cde *d;
	cgrad *gr, *gr1;
	int i0, ij, j, *vmi, xksave;
	real *Adjoints, *vscale;
	size_t L;

	if (nerror && *nerror >= 0) {
		err_jmp = &err_jmp0;
		ij = setjmp(err_jmp0.jb);
		if ((*nerror = ij))
			return;
		}
	errno = 0;	/* in case f77 set errno opening files */
	if (!asl->i.x_known) {
		co_index = i;
		x0_check_ASL(asl,X);
		}
	if ((!asl->i.ncxval || asl->i.ncxval[i] != asl->i.nxval)
	 && (!(x0kind & ASL_have_conval)
	     || i < n_conjac[0] || i >= n_conjac[1])) {
		xksave = asl->i.x_known;
		asl->i.x_known = 1;
		con1ival_nomap_ASL((ASL*)asl,i,X,nerror);
		asl->i.x_known = xksave;
		if (nerror && *nerror)
			return;
		}
	if (asl->i.Derrs)
		deriv_errchk_ASL((ASL*)asl, nerror, i, 1);
	if (!(x0kind & ASL_have_funnel)) {
		if (f_b)
			funnelset_ASL(asl, f_b);
		if (f_c)
			funnelset_ASL(asl, f_c);
		x0kind |= ASL_have_funnel;
		}
	Adjoints = adjoints;
	d = &con_de[i];
	gr1 = asl->i.Cgrad0[i];
	for(gr = gr1; gr; gr = gr->next)
		Adjoints[gr->varno] = gr->coef;
	if ((L = d->zaplen)) {
		memset(adjoints_nv1, 0, L);
		derprop(d->d);
		}
	vmi = 0;
	if (asl->i.vmap)
		vmi = get_vminv_ASL((ASL*)asl);
	if ((vscale = asl->i.vscale)) {
		if (vmi)
			for(gr = gr1; gr; gr = gr->next) {
				i0 = gr->varno;
				Adjoints[i0] *= vscale[vmi[i0]];
				}
		else
			for(gr = gr1; gr; gr = gr->next) {
				i0 = gr->varno;
				Adjoints[i0] *= vscale[i0];
				}
		}
	gr = gr1;
	i0 = 0;
	switch(asl->i.congrd_mode) {
	  case 1:
		for(; gr; gr = gr->next)
			G[i0++] = Adjoints[gr->varno];
		break;
	  case 2:
		for(; gr; gr = gr->next)
			G[gr->goff] = Adjoints[gr->varno];
		break;
	  default:
		if (vmi) {
			for(; gr; gr = gr->next) {
				i = vmi[j = gr->varno];
				while(i0 < i)
					G[i0++] = 0;
				G[i] = Adjoints[j];
				i0 = i + 1;
				}
			}
		else
			for(; gr; gr = gr->next) {
				i = gr->varno;
				while(i0 < i)
					G[i0++] = 0;
				G[i] = Adjoints[i];
				i0 = i + 1;
				}
		i = n_var;
		while(i0 < i)
			G[i0++] = 0;
	  }
	err_jmp = 0;
	}

 void
con1grd_nomap_ASL(ASL *a, int i, real *X, real *G, fint *nerror)
{
	ASL_fg *asl;
	static char who[] = "con1grd_nomap";

	INchk(a, who, i, a->i.n_con0);
	asl = (ASL_fg*)a;
	if (!want_derivs)
		No_derivs_ASL(who);
	Congrd1(asl, i, X, G, nerror);
	}

 void
con1grd_ASL(ASL *a, int i, real *X, real *G, fint *nerror)
{
	ASL_fg *asl;
	int *cm;
	static char who[] = "con1grd";

	INchk(a, who, i, a->i.n_con_);
	asl = (ASL_fg*)a;
	if (!want_derivs)
		No_derivs_ASL(who);
	if ((cm = a->i.cmap))
		i = cm[i];
	Congrd1(asl, i, X, G, nerror);
	}
