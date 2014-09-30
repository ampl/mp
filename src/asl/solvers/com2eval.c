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

#include "jac2dim.h"

 void
com2eval_ASL(ASL_fgh *asl, int i, int ie)
{
	cexp *c, *ce;
	expr *e;
	expr_v *V = var_ex + i;
	linpart *L, *Le;
	real t;

	c = cexps + i;
	ce = cexps + ie;
	do {
		e = c->e;
		t = (*e->op)(e C_ASL);
		if ((L = c->L))
			for(Le = L + c->nlin; L < Le; L++)
				t += L->fac * ((expr_v*)L->v.vp)->v;
		(*V++).v = t;
		}
		while(++c < ce);
	}

 void
com1eval(ASL_fgh *asl, int i, int n)
{
	cexp1 *c, *ce;
	expr *e;
	expr_v *V;
	linpart *L, *Le;
	real t;

	V = var_ex1 + i;
	c = cexps1 + i;
	ce = c + n;
	do {
		e = c->e;
		t = (*e->op)(e C_ASL);
		if ((L = c->L))
			for(Le = L + c->nlin; L < Le; L++)
				t += L->fac * ((expr_v*)L->v.vp)->v;
		(*V++).v = t;
		}
		while(++c < ce);
	}

 void
funnelset(ASL_fgh *asl, funnel *f)
{
	cplist	*cl;
	derp	*d;

	for(; f; f = f->next) {
		memset(adjoints_nv1, 0, f->fcde.zaplen);
		cl = f->cl;
		do *cl->ca.rp = 0;
			while((cl = cl->next));
		d = f->fcde.d;
		*d->b.rp = 1.;
		do *d->a.rp += *d->b.rp * *d->c.rp;
			while((d = d->next));
		cl = f->cl;
		do *cl->cfa = *cl->ca.rp;
			while((cl = cl->next));
		}
	}

 static void
hv_fwd(expr *e)
{
	argpair *da, *dae;
	expr *e1, **ep;
	real dO;

	for(; e; e = e->fwd) {
	    e->aO = e->adO = 0;
	    switch(e->a) {

		case Hv_timesR:
		case Hv_binaryR:
			e->dO.r = e->R.e->dO.r * e->dR;
			break;

		case Hv_timesLR:
		case Hv_binaryLR:
			e->dO.r = e->L.e->dO.r*e->dL + e->R.e->dO.r*e->dR;
			break;

		case Hv_timesL:
		case Hv_unary:
			e->dO.r = e->L.e->dO.r * e->dL;
			break;

		case Hv_vararg:
		case Hv_if:
			if ((e1 = ((expr_va *)e)->valf)) {
				hv_fwd(e1);
				e->dO.r = ((expr_va *)e)->vale->dO.r;
				}
			else if ((e1 = ((expr_va *)e)->val))
				e->dO.r = e1->dO.r;
			else
				e->dO.r = 0;
			break;

		case Hv_plterm:
			e->dO.r = e->dL * e->R.e->dO.r;
			break;

		case Hv_sumlist:
			ep = e->R.ep;
			for(dO = 0; (e1 = *ep); ep++)
				dO += e1->dO.r;
			e->dO.r = dO;
			break;

		case Hv_func:
			da = ((expr_f *)e)->da;
			dae = ((expr_f *)e)->dae;
			for(dO = 0; da < dae; da++)
				dO += da->e->dO.r * *da->u.v;
			e->dO.r = dO;
			break;

		case Hv_negate:
			e->dO.r = -e->L.e->dO.r;
			break;

		case Hv_plusR:
			e->dO.r = e->R.e->dO.r;
			break;

		case Hv_plusL:
			e->dO.r = e->L.e->dO.r;
			break;

		case Hv_plusLR:
			e->dO.r = e->L.e->dO.r + e->R.e->dO.r;
			break;

		case Hv_minusR:
			e->dO.r = -e->R.e->dO.r;
			break;

		case Hv_minusLR:
			e->dO.r = e->L.e->dO.r - e->R.e->dO.r;
			break;

		default:/*DEBUG*/
			fprintf(Stderr, "bad e->a = %d in hv_fwd\n", e->a);
			exit(1);
		}
	    }
	}

 static void
func_back(expr_f *f)
{
	argpair *da, *da1, *dae;
	expr *e;
	real **fh;
	real aO, adO, t;

	fh = f->fh;
	aO = f->aO;
	adO = f->adO;
	da = f->da;
	for(dae = f->dae; da < dae; da++) {
		e = da->e;
		t = *da->u.v;
		e->aO += t*aO;
		e->adO += t*adO;
		for(da1 = f->da; da1 < dae; da1++) {
			e = da1->e;
			e->aO += adO*e->dO.r * **fh++;
			}
		}
	}

 static void
hv_back(expr *e)
{
	expr *e1, **ep, *e2;
	real adO, t1, t2;

	if (!e || (!e->aO && !e->adO))
		return;
	for(; e; e = e->bak)
	    switch(e->a) {
		case Hv_binaryR:
			e1 = e->R.e;
			e1->adO += e->adO * e->dR;
			e1->aO += e->aO * e->dR  +  e->adO * e1->dO.r * e->dR2;
			break;

		case Hv_binaryLR:
			e1 = e->L.e;
			e2 = e->R.e;
			adO = e->adO;
			t1 = adO * e1->dO.r;
			t2 = adO * e2->dO.r;
			e1->aO  += e->aO*e->dL + t1*e->dL2 + t2*e->dLR;
			e2->aO  += e->aO*e->dR + t1*e->dLR + t2*e->dR2;
			e1->adO += adO * e->dL;
			e2->adO += adO * e->dR;
			break;

		case Hv_unary:
			e1 = e->L.e;
			e1->adO += e->adO * e->dL;
			e1->aO  += e->aO * e->dL  +  e->adO * e1->dO.r * e->dL2;
			break;

		case Hv_vararg:
		case Hv_if:
			if ((e1 = ((expr_va *)e)->vale)) {
				e1->aO = e->aO;
				e1->adO = e->adO;
				hv_back(e1);
				}
			else {
				e1 = ((expr_va *)e)->val;
				if (e1->op != f_OPNUM) {
					e1->aO = e->aO;
					e1->adO = e->adO;
					}
				}
			break;

		case Hv_plterm:
			e->R.e->aO += e->dL * e->aO;
			break;

		case Hv_sumlist:
			ep = e->R.ep;
			t1 = e->aO;
			t2 = e->adO;
			while((e1 = *ep++)) {
				e1->aO += t1;
				e1->adO += t2;
				}
			break;

		case Hv_func:
			func_back((expr_f *)e);
			break;

		case Hv_negate:
			e1 = e->L.e;
 neg_end:
			e1->aO -= e->aO;
			e1->adO -= e->adO;
			break;

		case Hv_plusR:
			e1 = e->R.e;
			goto plus_end;

		case Hv_plusL:
			e1 = e->L.e;
 plus_end:
			e1->aO += e->aO;
			e1->adO += e->adO;
			break;

		case Hv_plusLR:
			e1 = e->L.e;
			e1->aO += t1 = e->aO;
			e1->adO += t2 = e->adO;
			e1 = e->R.e;
			e1->aO += t1;
			e1->adO += t2;
			break;

		case Hv_minusR:
			e1 = e->R.e;
			goto neg_end;

		case Hv_minusLR:
			e1 = e->L.e;
			e1->aO += t1 = e->aO;
			e1->adO += t2 = e->adO;
			e1 = e->R.e;
			e1->aO -= t1;
			e1->adO -= t2;
			break;

		case Hv_timesR:
			e1 = e->R.e;
			e1->adO += e->adO * (t1 = e->dR);
			e1->aO += e->aO * t1;
			break;

		case Hv_timesL:
			e1 = e->L.e;
			e1->adO += e->adO * (t1 = e->dL);
			e1->aO  += e->aO * t1;
			break;

		case Hv_timesLR:
			e1 = e->L.e;
			e2 = e->R.e;
			adO = e->adO;
			e1->aO  += e->aO*e->dL  +  adO * e2->dO.r;
			e2->aO  += e->aO*e->dR  +  adO * e1->dO.r;
			e1->adO += adO * e->dL;
			e2->adO += adO * e->dR;
			break;

		default:/*DEBUG*/
			fprintf(Stderr, "bad e->a = %d in hv_back\n", e->a);
			exit(1);
		}
	}

 static void
hv_fwd0(cexp *c, expr_v *v)
{
	linpart *L, *Le;
	real x;

	v->aO = v->adO = 0;
	if (c->ef) {
		hv_fwd(c->ef);
		x = c->ee->dO.r;
		}
	else if (c->e->op != f_OPNUM)
		x = c->e->dO.r;
	else
		x = 0;
	if ((L = c->L))
		for(Le = L + c->nlin; L < Le; L++)
			x += L->fac * ((expr_v*)L->v.vp)->dO.r;
	v->dO.r = x;
	}

 static void
hv_fwd1(ASL_fgh *asl, cde *d)
{
	cexp1 *c, *ce;
	expr_v *v;
	int i;

	i = d->com11;
	v = var_ex1 + i;
	c = cexps1 + i;
	ce = c + d->n_com1;
	do hv_fwd0((cexp *)c, v++);
		while(++c < ce);
	}

 static void
hv_back0(ASL_fgh *asl, int i, int n)
{
	cexp *c, *ce;
	expr *e;
	expr_v *v;
	linpart *L, *Le;
	real t;

	v = var_ex + (i + n);
	ce = cexps + i;
	c = ce + n;
	do {
		--c;
		--v;
		if ((t = v->aO) && (L = c->L))
			for(Le = L + c->nlin; L < Le; L++)
				((expr_v*)L->v.vp)->aO += t * L->fac;
		if ((e = c->ee)) {
			e->aO = t;
			e->adO = v->adO;
			hv_back(e);
			}
		else if ((e = c->e)->op != f_OPNUM) {
			e->aO = t;
			e->adO = v->adO;
			}
		}
		while(c != ce);
	}

 static void
hv_back1(ASL_fgh *asl, cde *d)
{
	cexp1 *c, *ce;
	expr *e;
	expr_v *v;
	int i;
	linpart *L, *Le;
	real t;

	i = d->com11 + d->n_com1;
	v = var_ex1 + i;
	c = cexps1 + i;
	ce = cexps1 + d->com11;
	do {
		--c;
		--v;
		if ((t = v->aO) && (L = c->L))
			for(Le = L + c->nlin; L < Le; L++)
				((expr_v*)L->v.vp)->aO += t * L->fac;
		if ((e = c->ee)) {
			e->aO = t;
			e->adO = v->adO;
			hv_back(e);
			}
		else if ((e = c->e)->op != f_OPNUM) {
			e->aO = t;
			e->adO = v->adO;
			}
		}
		while(c != ce);
	}

#undef asl

 void
hv2comp_ASL(ASL *a, real *hv, real *p, int nobj, real *ow, real *y0)
	/* p = direction */
	/* y = Lagrange multipliers */
	/* hv = result */
{
	ASL_fgh *asl;
	cde *d, *d0;
	cexp *c, *ce;
	expr *e;
	expr_v *x, *xe;
	int n, no, noe;
	real *cscale, t, *vscale, *y, *ye, yi;

	ASL_CHECK(a, ASL_read_fgh, "hv2comp");
	asl = (ASL_fgh*)a;

	if (nobj >= 0 && nobj < n_obj) {
		no = nobj;
		noe = no + 1;
		ow = &edag_one_ASL;
		}
	else {
		nobj = -1;
		no = noe = 0;
		if (ow)
			noe = n_obj;
		}

	n = n_var;
	for(x = var_e, xe = x + n; x < xe; x++) {
		x->dO.r = *p++;
		x->aO = x->adO = 0;
		}
	if ((vscale = asl->i.vscale))
		for(x = var_e; x < xe; x++)
			x->dO.r *= *vscale++;
	x = var_ex;
	if (comb)
		for(c = cexps, ce = cexpsc; c < ce; c++)
			hv_fwd0(c, x++);
	if ((y = y0)) {
		ye = y + n_con;
		d = con_de;
		for(y0 = 0; y < ye; y++, d++)
			if ((e = d->ef) && *y) {
				if (!y0) {
					if (comc) {
						c = cexpsc;
						ce = cexpso;
						for(; c < ce; c++)
							hv_fwd0(c, x++);
						}
					y0 = y;
					d0 = d;
					}
				if (d->n_com1)
					hv_fwd1(asl, d);
				hv_fwd(e);
				}
		}
	for(; no < noe; no++)
	    if ((t = *ow++)) {
		d = obj_de + no;
		if (cexpso < cexpse) {
			x = var_ex + combc;
			c = cexpso;
			ce = cexpse;
			do hv_fwd0(c, x++);
				while(++c < ce);
			}
		if (d->n_com1)
			hv_fwd1(asl, d);
		if ((e = d->ef)) {
			hv_fwd(e);
			e = d->ee;
			e->aO = 0;
			e->adO = t;
			hv_back(e);
			}
		else if ((e = d->e)->op != f_OPNUM) {
			e->aO = 0;
			e->adO = t;
			}
		if (d->n_com1)
			hv_back1(asl, d);
		if (como)
			hv_back0(asl, combc, como);
		}
	if (y0) {
		if ((cscale = a->i.lscale))
			cscale += d0 - con_de;
		do {
			if ((yi = *y0++)) {
				if (cscale)
					yi *= *cscale++;
				if (!(e = d0->ee)) {
					if ((e = d0->e)->op != f_OPNUM) {
						e->aO = 0;
						e->adO = yi;
						}
					}
				else if ((e->adO = yi)) {
					e->aO = 0;
					hv_back(e);
					if (d0->n_com1)
						hv_back1(asl, d0);
					}
				}
			else if (cscale)
				++cscale;
			d0++;
			}
			while(y0 < ye);
		if (comc)
			hv_back0(asl, comb, comc);
		}
	if (comb)
		hv_back0(asl, 0, comb);
	x = var_e;
	if (vscale)
		for(vscale = asl->i.vscale; x < xe; x++)
			*hv++ = x->aO * *vscale++;
	else
		for(; x < xe; x++)
			*hv++ = x->aO;
	}

 void
hv2compd_ASL(ASL *a, real *hv, real *p, int co)
	/* p = direction */
	/* hv = result */
	/* co >= 0: behave like hv2comp_ASL with nobj = -1, ow = 0, y[i] = i == co ? 1. : 0. */
	/* co < 0: behave like hv2comp_ASL with nobj = -1 - co, ow = 0, y = 0 */
{
	ASL_fgh *asl;
	cde *d;
	cexp *c, *ce;
	cgrad *cg, *cg0;
	expr *e;
	expr_v *x, *x0;
	int n, no;
	ograd *og;
	real *cscale, t, *vscale;
	varno_t i;

	ASL_CHECK(a, ASL_read_fgh, "hv2comp");
	asl = (ASL_fgh*)a;

	n = n_var;
	memset(hv, 0, n*sizeof(real));

	no = -1 - co;
	if (co >= nlc || no >= nlo)
		return;
	cg0 = 0;
	vscale = asl->i.vscale;
	x0 = var_e;
	if (co >= 0) {
		cg = cg0 = Cgrad[co];
		if (vscale) {
			for(; cg; cg = cg->next) {
				i = cg->varno;
				x = x0 + i;
				x->dO.r = p[i]*vscale[i];
				x->aO = x->adO = 0.;
				}
			}
		else {
			for(; cg; cg = cg->next) {
				i = cg->varno;
				x = x0 + i;
				x->dO.r = p[i];
				x->aO = x->adO = 0.;
				}
			}
		}
	else {
		og = Ograd[no];
		if (vscale) {
			for(; og; og = og->next) {
				i = og->varno;
				x = x0 + i;
				x->dO.r = p[i]*vscale[i];
				x->aO = x->adO = 0.;
				}
			}
		else {
			for(; og; og = og->next) {
				i = og->varno;
				x = x0 + i;
				x->dO.r = p[i];
				x->aO = x->adO = 0.;
				}
			}
		}
	x = var_ex;
	if (comb)
		for(c = cexps, ce = cexpsc; c < ce; c++)
			hv_fwd0(c, x++);
	if (co >= 0) {
		if (comc) {
			c = cexpsc;
			ce = cexpso;
			for(; c < ce; c++)
				hv_fwd0(c, x++);
			}
		d = con_de + co;
		t = 1.;
		if ((cscale = a->i.lscale))
			t = cscale[co];
		}
	else {
		d = obj_de + no;
		if (cexpso < cexpse) {
			x = var_ex + combc;
			c = cexpso;
			ce = cexpse;
			do hv_fwd0(c, x++);
				while(++c < ce);
			}
		t = 1.;
		}
	if (d->n_com1)
		hv_fwd1(asl, d);
	if ((e = d->ef)) {
		hv_fwd(e);
		e = d->ee;
		e->aO = 0;
		e->adO = t;
		hv_back(e);
		}
	else if ((e = d->e)->op != f_OPNUM) {
		e->aO = 0;
		e->adO = t;
		}
	if (d->n_com1)
		hv_back1(asl, d);
	if (co >= 0) {
		if (comc)
			hv_back0(asl, comb, comc);
		}
	else {
		if (como)
			hv_back0(asl, combc, como);
		}
	if (comb)
		hv_back0(asl, 0, comb);
	if ((cg = cg0)) {
		if (vscale) {
			while(cg) {
				i = cg->varno;
				hv[i] = vscale[i]*x0[i].aO;
				cg = cg->next;
				}
			}
		else {
			while(cg) {
				i = cg->varno;
				hv[i] = x0[i].aO;
				cg = cg->next;
				}
			}
		}
	else {
		og = Ograd[no];
		if (vscale) {
			while(og) {
				i = og->varno;
				hv[i] = vscale[i]*x0[i].aO;
				og = og->next;
				}
			}
		else {
			while(og) {
				i = og->varno;
				hv[i] = x0[i].aO;
				og = og->next;
				}
			}
		}
	}

 varno_t
hv2comps_ASL(ASL *a, real *hv, real *p, int co, varno_t nz, varno_t *z)
	/* p = direction */
	/* hv = result */
	/* co >= 0: behave like hv2comp_ASL with nobj = -1, ow = 0, y[i] = i == co ? 1. : 0. */
	/* co < 0: behave like hv2comp_ASL with nobj = -1 - co, ow = 0, y = 0 */
{
	ASL_fgh *asl;
	cde *d;
	cexp *c, *ce;
	cgrad *cg, *cg0;
	expr *e;
	expr_v *x, *x0;
	int n, no;
	ograd *og;
	real *cscale, *hve, t, *vscale;
	varno_t f, i, rv, *ze;

	ASL_CHECK(a, ASL_read_fgh, "hv2comp");
	asl = (ASL_fgh*)a;

	n = n_var;
	memset(hv, 0, n*sizeof(real));

	no = -1 - co;
	if (co >= nlc || no >= nlo)
		return 0;
	cg0 = 0;
	vscale = asl->i.vscale;
	x0 = var_e;
	if (co >= 0) {
		cg = cg0 = Cgrad[co];
		if (vscale) {
			for(; cg; cg = cg->next) {
				i = cg->varno;
				x = x0 + i;
				x->dO.r = p[i]*vscale[i];
				x->aO = x->adO = 0.;
				}
			}
		else {
			for(; cg; cg = cg->next) {
				i = cg->varno;
				x = x0 + i;
				x->dO.r = p[i];
				x->aO = x->adO = 0.;
				}
			}
		}
	else {
		og = Ograd[no];
		if (vscale) {
			for(; og; og = og->next) {
				i = og->varno;
				x = x0 + i;
				x->dO.r = p[i]*vscale[i];
				x->aO = x->adO = 0.;
				}
			}
		else {
			for(; og; og = og->next) {
				i = og->varno;
				x = x0 + i;
				x->dO.r = p[i];
				x->aO = x->adO = 0.;
				}
			}
		}
	x = var_ex;
	if (comb)
		for(c = cexps, ce = cexpsc; c < ce; c++)
			hv_fwd0(c, x++);
	if (co >= 0) {
		if (comc) {
			c = cexpsc;
			ce = cexpso;
			for(; c < ce; c++)
				hv_fwd0(c, x++);
			}
		d = con_de + co;
		t = 1.;
		if ((cscale = a->i.lscale))
			t = cscale[co];
		}
	else {
		d = obj_de + no;
		if (cexpso < cexpse) {
			x = var_ex + combc;
			c = cexpso;
			ce = cexpse;
			do hv_fwd0(c, x++);
				while(++c < ce);
			}
		t = 1.;
		}
	if (d->n_com1)
		hv_fwd1(asl, d);
	if ((e = d->ef)) {
		hv_fwd(e);
		e = d->ee;
		e->aO = 0;
		e->adO = t;
		hv_back(e);
		}
	else if ((e = d->e)->op != f_OPNUM) {
		e->aO = 0;
		e->adO = t;
		}
	if (d->n_com1)
		hv_back1(asl, d);
	if (co >= 0) {
		if (comc)
			hv_back0(asl, comb, comc);
		}
	else {
		if (como)
			hv_back0(asl, combc, como);
		}
	if (comb)
		hv_back0(asl, 0, comb);
	rv = 0;
	f = Fortran;
	if ((ze = z))
		ze += nz;
	if ((hve = hv))
		hve += nz;
	if ((cg = cg0)) {
		if (!hv) {
			while(cg) {
				++rv;
				if (z < ze)
					*z++ = cg->varno;
				cg = cg->next;
				}
			}
		else if (vscale) {
			while(cg) {
				++rv;
				i = cg->varno;
				if (z < ze)
					*z++ = f + i;
				if (hv < hve)
					*hv++ = vscale[i]*x0[i].aO;
				cg = cg->next;
				}
			}
		else {
			while(cg) {
				++rv;
				i = cg->varno;
				if (z < ze)
					*z++ = f + i;
				if (hv < hve)
					*hv++ = x0[i].aO;
				cg = cg->next;
				}
			}
		}
	else {
		og = Ograd[no];
		if (!hv) {
			while(og) {
				++rv;
				if (z < ze)
					*z++ = og->varno;
				og = og->next;
				}
			}
		else if (vscale) {
			while(og) {
				++rv;
				i = og->varno;
				if (z < ze)
					*z++ = f + i;
				if (hv < hve)
					*hv++ = vscale[i]*x0[i].aO;
				og = og->next;
				}
			}
		else {
			while(og) {
				++rv;
				i = og->varno;
				if (z < ze)
					*z++ = f + i;
				if (hv < hve)
					*hv++ = x0[i].aO;
				og = og->next;
				}
			}
		}
	return rv;
	}
