/****************************************************************
Copyright (C) 1997, 1999-2001 Lucent Technologies
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

#include "jacpdim.h"

#ifdef __cplusplus
extern "C" {
#endif

extern void funpset_ASL ANSI((ASL_pfgh*, funnel*));
#undef funnelset
#define funnelset(x) funpset_ASL(asl, x)

 void
ihd_clear_ASL(ASL_pfgh *asl)
{
	Ihinfo *ihi;
	int i = asl->P.ihdcur;
	asl->P.ihdcur = 0;
	for(ihi = asl->P.ihi1; ihi->ihd <= i; ihi = ihi->next) {
		del_mblk(ihi->k, ihi->hest);
		ihi->hest = 0;
		}
	}

 int
xp_check_ASL(ASL_pfgh *asl, real *x)
{
	cexp *c, *c1, *ce;
	expr *e;
	expr_v *v, *v0;
	int *dvsp0, i, i0, i1, *vm;
	linarg *la;
	linpart *L, *Le;
	ograd *og;
	real t, *vscale, *xe;

	if (x0kind != ASL_first_x && !memcmp(Lastx, x, x0len))
		return 0;

	if (asl->i.Derrs)
		deriv_errclear_ASL(&asl->i);
	want_deriv = want_derivs;
	memcpy(Lastx, x, x0len);
	asl->i.nxval++;
	if (asl->P.ihdcur)
		ihd_clear_ASL(asl);
	x0kind = asl->I.x0kind_init;
	xe = x + n_var;
	v = v0 = var_e;
	vscale = asl->i.vscale;
	if ((vm = asl->i.vmap)) {
		if (vscale)
			while(x < xe)
				v0[*vm++].v = *vscale++ * *x++;
		else
			while(x < xe)
				v0[*vm++].v = *x++;
		}
	else {
		if (vscale)
			while(x < xe)
				(v++)->v = *vscale++ * *x++;
		else
			while(x < xe)
				(v++)->v = *x++;
		}

	for(la = asl->P.lalist; la; la = la->lnext) {
		og = la->nz;
		t = var_e[og->varno].v*og->coef;
		while((og = og->next))
			t += var_e[og->varno].v*og->coef;
		la->v->v = t;
		}
	errno = 0;
	if (asl->P.ncom) {
		c = cexps;
		dvsp0 = asl->P.dvsp0;
		i0 = *dvsp0;
		i = 0;
		/* Normally v == var_ex here, but v < var_ex is possible if sos_add */
		/* overestimated the number of new variables that sos_finish added. */
		v = var_ex;
		for(ce = c1 = c + asl->P.ncom; c < ce; c++) {
			for(i1 = *++dvsp0; i0 < i1; i0++, c1++) {
				cv_index = i0 + 1;
				e = c1->e;
				asl->P.vp[i0]->v = (*e->op)(e C_ASL);
				if (c1->funneled)
					funnelset(c1->funneled);
				}
			e = c->e;
			cv_index = ++i;
			t = (*e->op)(e C_ASL);
			if ((L = c->L))
				for(Le = L + c->nlin; L < Le; L++)
					t += L->fac * ((expr_v*)L->v.vp)->v;
			else if (!c->d && (og = asl->P.dv[c-cexps].ll)) {
				if (og->varno < 0) {
					t += og->coef;
					og = og->next;
					}
				while(og) {
					t += og->coef*v0[og->varno].v;
					og = og->next;
					}
				}
			(v++)->v = t;
			if (c->funneled)
				funnelset(c->funneled);
			}
		cv_index = 0;
		}
	return 1;
	}

 void
xp2known_ASL(ASL* asl, real *X, fint *nerror)
{
	Jmp_buf err_jmp0;
	int ij;

	ASL_CHECK(asl, ASL_read_pfgh, "xp2known");
	if (asl->i.xknown_ignore)
		return;
	if (nerror && *nerror >= 0) {
		err_jmp = &err_jmp0;
		ij = setjmp(err_jmp0.jb);
		if ((*nerror = ij))
			goto done;
		}
	errno = 0;	/* in case f77 set errno opening files */
	xp_check_ASL((ASL_pfgh*)asl, X);
	asl->i.x_known = 1;
 done:
	err_jmp = 0;
	}

 static real *
bigUmult(ASL_pfgh *asl, real *h, range *r, int nobj, real *ow, real *y)
{
	real *s, t;
	Umultinfo *u, *u0, *u1, *ue, **utodo, **utodoi;
	int i, j, ku, n, nv;
	int *imap, *iv;
	linarg *la, **lap;
	ograd *og;

	s = asl->P.dOscratch;
	utodo = utodoi = (Umultinfo**)asl->P.utodo;
	n = r->n;
	ku = htcl(n*sizeof(Umultinfo) + n_var*sizeof(int));
	u = u0 = (Umultinfo*)new_mblk(ku);
	imap = (int*)(u + n);
	iv = r->ui;
	nv = r->nv;
	for(i = 0; i < nv; i++) {
		imap[j = *iv++] = i;
		utodo[j] = 0;
		}
	lap = r->lap;
	for(i = 0; i < n; i++) {
		la = *lap++;
		u->v = la->v;
		u->i = i;
		u->og = u->og0 = og = la->nz;
		utodoi = utodo + og->varno;
		u->next = *utodoi;
		*utodoi = u++;
		}
	ue = u;
	iv = r->ui;
	for(i = 0; i < nv; i++) {
		utodoi = utodo + *iv++;
		u1 = *utodoi;
		*utodoi = 0;
		for(u = u1; u; u = u->next)
			s[u->i] = u->og->coef;
		pshv_prod(r, nobj, ow, y);
		h += i;
		for(j = 0; j <= i; j++)
			h[j] = 0.;
		while((u = u1)) {
			u1 = u->next;
			s[u->i] = 0.;
			if ((og = u->og->next)) {
				u->og = og;
				utodoi = utodo + og->varno;
				u->next = *utodoi;
				*utodoi = u;
				}
			}
		for(u = u0; u < ue; u++)
			if ((t = u->v->aO))
				for(og = u->og0; og &&
						(j = imap[og->varno]) <= i;
						og = og->next)
					h[j] += t*og->coef;
		}
	del_mblk(ku, u0);
	return h + nv;
	}

 void
hvpinit_ASL(ASL *a, int ndhmax, int nobj, real *ow, real *y)
{
	ASL_pfgh *asl;
	Ihinfo *ihi;
	range *r;
	real *h, *s, *si;
	int i, ihc, n1;
	linarg **lap, **lap1, **lape;
	expr_v *v;

	ASL_CHECK(a, ASL_read_pfgh, "xvpinit");
	asl = (ASL_pfgh*)a;
	xpsg_check_ASL(asl, nobj, ow, y);
	asl->P.nhvprod = 0;
	if (!asl->P.hes_setup_called)
		(*asl->p.Hesset)(a, 1, 0, nlo, 0, nlc);
	ihc = 0;
	if (ndhmax > asl->P.ihdmax)
		ndhmax = asl->P.ihdmax;
	if ((asl->P.ndhmax = ndhmax) <= 0)
		goto done;
	if (!(ihi = asl->P.ihi1) || ndhmax < asl->P.ihdmin)
		return;
	if (nobj < 0 || nobj >= n_obj)
		nobj = -1;
	s = asl->P.dOscratch;
	for(ihc = 0; ihi->ihd <= ndhmax; ihi = ihi->next) {
		ihc = ihi->ihd;
		ihi->hest = h = (real *)new_mblk(ihi->k);
		for(r = ihi->r; r; r = r->rlist.prev) {
			r->hest = h;
			if ((n1 = r->n) < r->nv) {
				si = s;
				lape = lap = r->lap;
				for(i = 0; i < n1; i++) {
					*si = 1.;
					pshv_prod(r, nobj, ow, y);
					*si++ = 0;
					lape++;
					lap1 = lap;
					do {
						v = (*lap1++)->v;
						*h++ = v->aO;
						}
						while(lap1 < lape);
					}
				}
			else
				h = bigUmult(asl, h, r, nobj, ow, y);
			}
		}
 done:
	asl->P.ihdcur = ihc;
	}

#ifdef __cplusplus
}
#endif
