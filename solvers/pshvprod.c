/****************************************************************
Copyright (C) 1997, 1998, 2000, 2001 Lucent Technologies
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

 static void
#ifdef KR_headers
hv_fwd(e) register expr *e;
#else
hv_fwd(register expr *e)
#endif
{
	register expr *e1, **ep;
	argpair *da, *dae;
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
			if (e1 = ((expr_va *)e)->valf) {
				hv_fwd(e1);
				e->dO.r = ((expr_va *)e)->vale->dO.r;
				}
			else if ((e1 = ((expr_va *)e)->val) && e1->op != f_OPNUM)
				e->dO.r = e1->dO.r;
			else
				e->dO.r = 0;
			break;

		case Hv_plterm:
			e->dO.r = e->dL * e->R.e->dO.r;
			break;

		case Hv_sumlist:
			ep = e->R.ep;
			for(dO = 0; e1 = *ep; ep++)
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
#ifdef KR_headers
func_back(f) expr_f *f;
#else
func_back(expr_f *f)
#endif
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
		e->adO += (t = *da->u.v) * adO;
		e->aO += t * aO;
		t = adO*e->dO.r;
		for(da1 = f->da; da1 < dae; da1++) {
			e = da1->e;
			e->aO += t * **fh++;
			}
		}
	}

 static void
#ifdef KR_headers
funnel_back(asl, c, v, t) ASL_pfgh *asl; cexp *c; expr_v *v; real t;
#else
funnel_back(ASL_pfgh *asl, cexp *c, expr_v *v, real t)
#endif
{
	real aO, adO;
	real *g, *h;
	hes_fun *hf;
	expr_v **vp, **vp1, **vpe;
	ograd *og;

	aO = v->aO = t;
	adO = v->adO;
	hf = c->hfun;
	if (og = hf->og) {
		do {
			v = var_e + og->varno;
			v->adO += (t = og->coef) * adO;
			v->aO += t*aO;
			}
			while(og = og->next);
		return;
		}
	g = hf->grdhes;
	h = g + hf->n;
	vp = hf->vp;
	vpe = vp + hf->n;
	do {
		v = *vp++;
		v->adO += (t = *g++) * adO;
		v->aO += t*aO;
		t = adO * v->dO.r;
		vp1 = hf->vp;
		do (*vp1++)->aO += t * *h++;
		   while(vp1 < vpe);
		}
		while(vp < vpe);
	}

 static void
#ifdef KR_headers
hv_back(e) register expr *e;
#else
hv_back(register expr *e)
#endif
{
	register expr *e1, **ep, *e2;
	real adO, t1, t2;

	if (!e || !e->aO && !e->adO)
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
			if (e1 = ((expr_va *)e)->vale) {
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
			while(e1 = *ep++) {
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
#ifdef KR_headers
hv_fwd0(asl, c, v) ASL_pfgh *asl; cexp *c; expr_v *v;
#else
hv_fwd0(ASL_pfgh *asl, cexp *c, expr_v *v)
#endif
{
	register linpart *L, *Le;
	real *g, x;
	hes_fun *hf;
	ograd *og;
	expr_v **vp, **vpe;

	v->aO = v->adO = 0;
	if (hf = c->hfun) {
		x = 0;
		if (og = hf->og)
			do x += og->coef * var_e[og->varno].dO.r;
			while(og = og->next);
		else {
			g = hf->grdhes;
			vp = hf->vp;
			vpe = vp + hf->n;
			do x += *g++ * (*vp++)->dO.r;
			   while(vp < vpe);
			}
		}
	else if (c->ef) {
		hv_fwd(c->ef);
		x = c->ee->dO.r;
		}
	else if (c->e->op != f_OPNUM)
		x = c->e->dO.r;
	else
		x = 0;
	if (L = c->L)
		for(Le = L + c->nlin; L < Le; L++)
			x += L->fac * ((expr_v*)L->v.vp)->dO.r;
	v->dO.r = x;
	}

 static void
#ifdef KR_headers
hfg_fwd(e) register expr *e;
#else
hfg_fwd(register expr *e)
#endif
{
	expr *e1;

	for(; e; e = e->fwd) {
		e->aO = 0;
		switch(e->a) {
		  case Hv_vararg:
		  case Hv_if:
			if (e1 = ((expr_va *)e)->valf)
				hfg_fwd(e1);
		  }
		}
	}

 static void
#ifdef KR_headers
hfg_back(e) register expr *e;
#else
hfg_back(register expr *e)
#endif
{
	register expr *e1, **ep;
	real aO;

	if (!e || !e->aO && !e->adO)
		return;
	for(; e; e = e->bak)
	    switch(e->a) {
		case Hv_timesR:
		case Hv_binaryR:
			e->R.e->aO += e->aO * e->dR;
			break;

		case Hv_binaryLR:
		case Hv_timesLR:
			aO = e->aO;
			e->L.e->aO += aO * e->dL;
			e->R.e->aO += aO * e->dR;
			break;

		case Hv_unary:
		case Hv_timesL:
			e->L.e->aO += e->aO * e->dL;
			break;

		case Hv_vararg:
		case Hv_if:
			if (e1 = ((expr_va *)e)->vale) {
				e1->aO = e->aO;
				hfg_back(e1);
				}
			else {
				e1 = ((expr_va *)e)->val;
				if (e1->op != f_OPNUM)
					e1->aO = e->aO;
				}
			break;

		case Hv_plterm:
			e->R.e->aO += e->dL * e->aO;
			break;

		case Hv_sumlist:
			ep = e->R.ep;
			aO = e->aO;
			while(e1 = *ep++)
				e1->aO += aO;
			break;

		case Hv_func:
			func_back((expr_f *)e);
			break;

		case Hv_negate:
			e->L.e->aO -= e->aO;
			break;

		case Hv_plusR:
			e->R.e->aO += e->aO;
			break;

		case Hv_plusL:
			e->L.e->aO += e->aO;
			break;

		case Hv_plusLR:
			e->L.e->aO += aO = e->aO;
			e->R.e->aO += aO;
			break;

		case Hv_minusR:
			e->R.e->aO -= e->aO;
			break;

		case Hv_minusLR:
			e->L.e->aO += aO = e->aO;
			e->R.e->aO -= aO;
			break;

		default:/*DEBUG*/
			fprintf(Stderr, "bad e->a = %d in hfg_back\n", e->a);
			exit(1);
		}
	}

 static void
#ifdef KR_headers
funnelhes(asl) ASL_pfgh *asl;
#else
funnelhes(ASL_pfgh *asl)
#endif
{
	cexp *c;
	int n;
	hes_fun *hf;
	expr_v *v, **vp, **vp1, **vpe;
	real *g, *h;
	expr *e;

	x0kind &= ~ASL_need_funnel;
	for(hf = asl->I.hesthread; hf; hf = hf->hfthread) {
		if (hf->og)
			continue;
		n = hf->n;
		g = hf->grdhes;
		h = g + n;
		c = hf->c;
		vp = hf->vp;
		vpe = vp + n;

		do (*vp++)->aO = 0;
		   while(vp < vpe);

		hfg_fwd(c->ef);
		e = c->ee;
		e->aO = 1;
		hfg_back(e);

		vp = hf->vp;
		do {
			v = *vp++;
			*g++ = v->aO;
			v->dO.r = v->aO = v->adO = 0;
			}
			while(vp < vpe);

		vp = hf->vp;
		do {
			v = *vp++;
			v->dO.r = 1;
			if (e = c->ef)
				hv_fwd(e);
			if (e = c->ee) {
				e->aO = 0;
				e->adO = 1;
				hv_back(e);
				}
			else if ((e = c->e)->op != f_OPNUM) {
				e->aO = 0;
				e->adO = 1;
				}
			v->dO.r = 0;
			vp1 = hf->vp;
			do {
				v = *vp1++;
				*h++ = v->aO;
				v->aO = v->adO = 0;
				}
				while(vp1 < vpe);
			}
			while(vp < vpe);
		}
	}

 static void
#ifdef KR_headers
hvp0comp_ASL(a, hv, p, nobj, ow, y0)
	ASL *a; real *hv, *p, *ow, *y0; int nobj;
#else
hvp0comp_ASL(ASL *a, real *hv, real *p, int nobj, real *ow, real *y0)
#endif
	/* p = direction */
	/* y = Lagrange multipliers */
	/* hv = result */
{
	expr_v *v, *x, *x0, *xe;
	expr *e;
	real *cscale, t, t1, t2, *p1, *y, *ye, yi;
	cexp *c, *c1, *ce;
	int *dvsp0, i0, i1, n, no, noe;
	linpart *L, *Le;
	linarg *la;
	ograd *og;
	ps_func *f, *f0;
	psb_elem *b, *be;
	psg_elem *g, *ge;

#define asl ((ASL_pfgh*)a)
	if (x0kind & ASL_need_funnel)
		funnelhes(asl);
	if (nobj >= 0 && nobj < n_obj) {
		ow = ow ? ow + nobj : &edag_one_ASL;
		no = nobj;
		noe = no + 1;
		}
	else {
		no = noe = 0;
		if (ow)
			noe = n_obj;
		}
	n = c_vars >= o_vars ? c_vars : o_vars;
	for(la = asl->P.lalist; la; la = la->lnext) {
		og = la->nz;
		t = p[og->varno]*og->coef;
		while(og = og->next)
			t += p[og->varno]*og->coef;
		x = la->v;
		x->dO.r = t;
		x->aO = x->adO = 0;
		}
	p1 = p;
	x = x0 = var_e;
	for(xe = x + n; x < xe; x++) {
		x->dO.r = *p1++;
		x->aO = x->adO = 0;
		}
	if (asl->P.ncom) {
		x = var_ex;
		dvsp0 = asl->P.dvsp0;
		c = cexps;
		i0 = *dvsp0;
		for(ce = c1 = c + asl->P.ncom; c < ce; c++) {
			for(i1 = *++dvsp0; i0 < i1; i0++)
				hv_fwd0(asl, c1++, asl->P.vp[i0]);
			hv_fwd0(asl, c, x++);
			}
		}
	if (y = y0) {
		ye = y + n_con;
		f = asl->P.cps;
		for(y0 = 0; y < ye; y++, f++)
			if (*y) {
			    for(b = f->b, be = b + f->nb; b < be; b++)
				if (e = b->D.ef) {
					if (!y0) {
						y0 = y;
						f0 = f;
						}
					hv_fwd(e);
					}
			    for(g = f->g, ge = g + f->ng; g < ge; g++)
				for(b = g->E, be = b + g->ns; b < be; b++)
					if (e = b->D.ef) {
						if (!y0) {
							y0 = y;
							f0 = f;
							}
						hv_fwd(e);
						}
			    }
		}
	for(; no < noe; no++)
	    if (t2 = *ow++) {
		f = asl->P.ops + no;
		for(b = f->b, be = b + f->nb; b < be; b++)
			if (e = b->D.ef) {
				hv_fwd(e);
				e = b->D.ee;
				e->aO = 0;
				e->adO = t2;
				hv_back(e);
				}
			else if ((e = b->D.e)->op != f_OPNUM) {
				e->aO = 0;
				e->adO = t2;
				}
		for(g = f->g, ge = g + f->ng; g < ge; g++) {
			for(b = g->E, be = b + g->ns; b < be; b++)
				if (e = b->D.ef) {
					hv_fwd(e);
					e = b->D.ee;
					e->aO = 0;
					e->adO = t2*g->g1;
					hv_back(e);
					}
				else if ((e = b->D.e)->op != f_OPNUM) {
					e->aO = 0;
					e->adO = t2*g->g1;
					}
			if (g->g2) {
				t = 0.;
				for(og = g->og; og; og = og->next)
					t += og->coef * p[og->varno];
				t *= t2*g->g2;
				for(og = g->og; og; og = og->next)
					x0[og->varno].aO += t*og->coef;
				}
			}
		}
	if (y0) {
		if (cscale = a->i.lscale)
			cscale += f - asl->P.cps;
		do {
			yi = *y0++;
			if (cscale)
				yi *= *cscale++;
			for(b = f0->b, be = b + f0->nb; b < be; b++) {
				if (!(e = b->D.ee)) {
					if ((e = b->D.e)->op != f_OPNUM) {
						e->aO = 0;
						e->adO = yi;
						}
					}
				else if (e->adO = yi) {
					e->aO = 0;
					hv_back(e);
					}
				}
			for(g = f0->g, ge = g + f0->ng; g < ge; g++) {
				t = g->g1 * yi;
				for(b = g->E, be = b + g->ns; b < be; b++) {
				    if (!(e = b->D.ee)) {
					if ((e = b->D.e)->op != f_OPNUM) {
						e->aO = 0;
						e->adO = t;
						}
					}
				    else if (e->adO = t) {
					e->aO = 0;
					hv_back(e);
					}
				    }
				if (!(t = g->g2 * yi))
					continue;
				t1 = 0;
				for(og = g->og; og; og = og->next)
					t1 += og->coef * p[og->varno];
				t *= t1;
				for(og = g->og; og; og = og->next)
					x0[og->varno].aO += t*og->coef;
				}
			f0++;
			}
			while(y0 < ye);
		}
	if (asl->P.ncom)
	    for(c = cexps; c < ce--; ) {
		for(i0 = *--dvsp0; i0 < i1; ) {
			v = asl->P.vp[--i1];
			--c1;
			if ((t = v->aO) && (L = c1->L))
				for(Le = L + c1->nlin; L < Le; L++)
					((expr_v*)L->v.vp)->aO += t * L->fac;
			if (c1->hfun)
				funnel_back(asl, c1, v, t);
			else if (e = c1->ee) {
				e->aO = t;
				e->adO = v->adO;
				hv_back(e);
				}
			else if ((e = c1->e)->op != f_OPNUM) {
				e->aO = t;
				e->adO = v->adO;
				}
			}
		if ((t = (--x)->aO) && (L = ce->L))
			for(Le = L + ce->nlin; L < Le; L++)
				((expr_v*)L->v.vp)->aO += t * L->fac;
		if (ce->hfun)
			funnel_back(asl, ce, x, t);
		else if (e = ce->ee) {
			e->aO = t;
			e->adO = x->adO;
			hv_back(e);
			}
		else if ((e = ce->e)->op != f_OPNUM) {
			e->aO = t;
			e->adO = v->adO;
			}
		}
	x = var_e;
	for(la = asl->P.lalist; la; la = la->lnext)
		if (t = la->v->aO) {
			og = la->nz;
			do x[og->varno].aO += t*og->coef;
				while(og = og->next);
			}
	while(x < xe)
		*hv++ = (x++)->aO;
	}
#undef asl

 static real *	/* Compute vector x0 = mat(h)*y0,	*/
		/* where h = upper triang of mat(h).	*/
#ifdef KR_headers
dtmul(n, x0, h, y0) int n; real *x0; real *h; real *y0;
#else
dtmul(int n, real *x0, real *h, real *y0)
#endif
{
	int i;
	real *hi, t, *x, *y, *y1, yi;

	y1 = y0;
	--h;
	for(i = 0; i < n; i++) {
		hi = ++h + i;
		yi = *y1++;
		t = yi**hi;
		x = x0;
		y = y0;
		while(h < hi) {
			t += *y++**h;
			*x++ += yi**h++;
			}
		*x = t;
		}
	return x0;
	}

 void
#ifdef KR_headers
hvpcomp_ASL(a, hv, p, nobj, ow, y)
	ASL *a; real *hv, *p, *ow, *y; int nobj;
#else
hvpcomp_ASL(ASL *a, real *hv, real *p, int nobj, real *ow, real *y)
#endif
	/* p = direction */
	/* y = Lagrange multipliers */
	/* hv = result */
{
	int kp, kw, n, no, noe, ns, nv, *ui, *uie;
	Ihinfo *ihi;
	range *r;
	linarg *la, **lap, **lape;
	ograd *og;
	real *cscale, *owi, t, t1, t2, *p0, *s, *w, *wi, *x;
	psg_elem *g, *ge;
	ps_func *ps, *pe;

	ASL_CHECK(a, ASL_read_pfgh, "hvpcomp");
#define asl ((ASL_pfgh*)a)
	if (a->i.x_known != 2)
		xpsg_check_ASL(asl, nobj, ow, y);
	nv = n_var;
	kp = htcl(nv*sizeof(real));
	p0 = 0;
	if (s = asl->i.vscale) {
		p0 = (real*)new_mblk(kp);
		for(n = 0; n < nv; n++)
			p0[n] = s[n] * p[n];
		p = p0;
		}
	if (!asl->P.ihdcur) {
		if (asl->P.ihdmin <= 0) {
			hvp0comp_ASL(a,hv,p,nobj,ow,y);
			goto done;
			}
		if (!(n = asl->P.nhvprod))
			n = asl->P.ihdmin;
		if (n >= asl->P.ihdmin)
			hvpinit_ASL(a, ihd_limit, nobj, ow, y);
		}
	asl->P.nhvprod++;
	memset(hv, 0, nv*sizeof(real));
	kw = kp + 1;
	w = (real*)new_mblk(kw);
	x = w + n_var;
	s = asl->P.dOscratch;
	ns = 0;
	for(ihi = asl->P.ihi1; ihi; ihi = ihi->next) {
		r = ihi->r;
		if (ihi->hest)
		    for(; r; r = r->rlist.prev) {
			n = r->n;
			nv = r->nv;
			wi = w;
			if (n < nv) {
				lap = r->lap;
				lape = lap + n;
				do {
					og = (*lap++)->nz;
					t = p[og->varno]*og->coef;
					while(og = og->next)
						t += p[og->varno]*og->coef;
					*wi++ = t;
					}
					while(lap < lape);
				wi = dtmul(n, x, r->hest, w);
				lap = r->lap;
				do if (t = *wi++) {
					og = (*lap)->nz;
					do hv[og->varno] +=
						t*og->coef;
					   while(og = og->next);
					}
					while(++lap < lape);
				}
			else {
				ui = r->ui;
				uie = ui + nv;
				do *wi++ = p[*ui++];
					while(ui < uie);
				wi = dtmul(nv, x, r->hest, w);
				ui = r->ui;
				do hv[*ui++] += *wi++;
					while(ui < uie);
				}
			}
		else
		    for(; r; r = r->rlist.prev) {
			n = r->n;
			if (ns < n)
				ns = n;
			wi = s;
			lap = r->lap;
			lape = lap + n;
			do {
				og = (*lap++)->nz;
				t = p[og->varno]*og->coef;
				while(og = og->next)
					t += p[og->varno]*og->coef;
				*wi++ = t;
				}
				while(lap < lape);
			pshv_prod_ASL(asl, r, nobj, ow, y);
			lap = r->lap;
			do {
				la = *lap++;
				if (t = la->v->aO) {
					og = la->nz;
					do hv[og->varno] += t*og->coef;
						while(og = og->next);
					}
				}
				while(lap < lape);
			}
		}
	del_mblk(kw,w);
	wi = s + ns;
	while(wi > s)
		*--wi = 0.;
	if (asl->P.nobjgroups) {
	    if (nobj >= 0 && nobj < n_obj) {
		owi = ow ? ow + nobj : &edag_one_ASL;
		no = nobj;
		noe = no + 1;
		}
	    else {
		nobj = -1;
		no = noe = 0;
		if (owi = ow)
			noe = n_obj;
		}
	    for(; no < noe; no++)
		if (t = *owi++) {
		    ps = asl->P.ops + no;
		    g = ps->g;
		    for(ge = g + ps->ng; g < ge; g++)
			if ((t2 = g->g2) && (og = g->og)) {
				t1 = p[og->varno]*og->coef;
				while(og = og->next)
					t1 += p[og->varno]*og->coef;
				t2 *= t*t1;
				og = g->og;
				do hv[og->varno] += t2*og->coef;
					while(og = og->next);
				}
		}
	    }
	if (asl->P.ncongroups && y) {
		cscale = a->i.lscale;
		ps = asl->P.cps;
		for(pe = ps + n_con; ps < pe; ps++, y++)
		    if (t = cscale ? *cscale++ * *y : *y)
			for(g = ps->g, ge = g + ps->ng; g < ge; g++)
			    if ((t2 = g->g2) && (og = g->og)) {
				t1 = p[og->varno]*og->coef;
				while(og = og->next)
					t1 += p[og->varno]*og->coef;
				t2 *= t*t1;
				og = g->og;
				do hv[og->varno] += t2*og->coef;
					while(og = og->next);
				}
		}
 done:
	if (p0) {
		del_mblk(kp, p0);
		s = asl->i.vscale;
		w = hv + n_var;
		while(hv < w)
			*hv++ *= *s++;
		}
	}
#undef asl

 void
#ifdef KR_headers
pshv_prod_ASL(a, r, nobj, ow, y) ASL_pfgh *a; range *r; int nobj; real *ow, *y;
#else
pshv_prod_ASL(ASL_pfgh *a, range *r, int nobj, real *ow, real *y)
#endif
{
	int *cei, *cei0, *ceie, i;
	linarg *la, **lap, **lape;
	linpart *L, *Le;
	expr_v *v;
	ps_func *p;
	cexp *c;
	expr *e;
	real *cscale, *s, owi, t;
	psb_elem *b;
	psg_elem *g;

	cscale = a->i.lscale;
#define asl a
	if (nobj >= 0 && nobj < n_obj) {
		if (ow) {
			if ((owi = ow[nobj]) == 0.)
				nobj = -1;
			ow = 0;
			}
		else
			owi = 1;
		}
	if (x0kind & ASL_need_funnel)
		funnelhes(asl);
	s = asl->P.dOscratch;
	lap = r->lap;
	lape = lap + r->n;
	while(lap < lape) {
		la = *lap++;
		v = la->v;
		v->dO.r = *s++;
		v->adO = v->aO = 0.;
		}
	if (cei = cei0 = r->cei) {
		i = *cei0++;
		ceie = (cei = cei0) + i;
		do {
			i = *cei++;
			hv_fwd0(asl, cexps + i, asl->P.vp[i]);
			}
			while(cei < ceie);
		}
	for(b = r->refs; b; b = b->next) {
		if ((i = b->conno) < 0) {
			i = -2 - i;
			if (i == nobj)
				t = owi;
			else if (ow)
				t = ow[i];
			else
				continue;
			p = asl->P.ops;
			}
		else {
			if (!y || !(t = y[i]))
				continue;
			if (cscale)
				t *= cscale[i];
			p = asl->P.cps;
			}
		if (b->groupno) {
			p += i;
			g = p->g + (b->groupno - 1);
			if (a->P.pshv_g1)
				t *= g->g1;
			}
		if (e = b->D.ef) {
			hv_fwd(e);
			e = b->D.ee;
			e->aO = 0;
			e->adO = t;
			hv_back(e);
			}
		else if ((e = b->D.e)->op != f_OPNUM)
			e->adO += t;
		}
	while(cei > cei0) {
		i = *--cei;
		c = cexps + i;
		v = asl->P.vp[i];
		if ((t = v->aO) && (L = c->L))
		    for(Le = L + c->nlin; L < Le; L++)
			((expr_v*)L->v.vp)->aO += t * L->fac;
		if (c->hfun)
			funnel_back(asl, c, v, t);
		else if (e = c->ee) {
			e->aO = t;
			e->adO = v->adO;
			hv_back(e);
			}
		else if ((e = c->e)->op != f_OPNUM) {
			e->aO += t;
			e->adO += v->adO;
			}
		}
	}

 void
#ifdef KR_headers
funpset_ASL(asl, f) ASL_pfgh *asl; register funnel *f;
#else
funpset_ASL(ASL_pfgh *asl, register funnel *f)
#endif
{
	register derp	*d;
	register cplist	*cl;

	for(; f; f = f->next) {
		memset(adjoints_nv1, 0, f->fcde.zaplen);
		cl = f->cl;
		do *cl->ca.rp = 0;
			while(cl = cl->next);
		d = f->fcde.d;
		*d->b.rp = 1.;
		do *d->a.rp += *d->b.rp * *d->c.rp;
			while(d = d->next);
		cl = f->cl;
		do *cl->cfa = *cl->ca.rp;
			while(cl = cl->next);
		}
	}

#ifdef __cplusplus
}
#endif
