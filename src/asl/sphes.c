/****************************************************************
Copyright (C) 1997-1998, 2000-2001 Lucent Technologies
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

#include "asl_pfgh.h"

 static void
hv_fwd(ASL_pfgh *asl, register expr *e)
{
	expr *e1, **ep;
	real dO;
	de *d;
	derp *D, *D1;

	for(; e; e = e->fwd) {
	    e->aO = e->adO = 0;
	    switch(e->a) {

		case Hv_timesR:
		case Hv_binaryR:
			e->dO.r = e->R.e->dO.r;
			break;

		case Hv_timesLR:
		case Hv_binaryLR:
			e->dO.r = e->L.e->dO.r + e->R.e->dO.r;
			break;

		case Hv_timesL:
		case Hv_unary:
			e->dO.r = e->L.e->dO.r;
			break;

		case Hv_vararg:
			e->dO.r = 0;
			for(d = ((expr_va*)e)->L.d; d->e; d++) {
				if ((e1 = d->ef)) {
					hv_fwd(asl, e1);
					e->dO.r += d->ee->dO.r;
					}
				else if ((e1 = d->e))
					e->dO.r += e1->dO.r;
				}
			if (!((expr_va*)e)->val && (D = ((expr_va*)e)->R.D)) {
				D->a.rp = adjoints_nv1;
				D1 = ((expr_va*)e)->d0;
				d = ((expr_va*)e)->L.d;
				((expr_va*)e)->val = d->e;
				do {
					D->next = d->d;
					while(D->next != D1)
						D = D->next;
					d->dlast = D;
					++d;
					}
					while(d->e);
				((expr_va*)e)->next = asl->P.valist;
				asl->P.valist = (expr_va*)e;
				}
			break;

		case Hv_if:
			e->dO.r = 0;
			if ((e1 = ((expr_if*)e)->Tf)) {
				hv_fwd(asl, e1);
				e->dO.r = ((expr_if *)e)->Te->dO.r;
				}
			else if ((e1 = ((expr_if*)e)->T))
				e->dO.r = e1->dO.r;
			if ((e1 = ((expr_if*)e)->Ff)) {
				hv_fwd(asl, e1);
				e->dO.r += ((expr_if *)e)->Fe->dO.r;
				}
			else if ((e1 = ((expr_if*)e)->F))
				e->dO.r += e1->dO.r;
			if (!((expr_if*)e)->val && (D = ((expr_if*)e)->D)) {
				((expr_if*)e)->val = ((expr_if*)e)->T;
				D->a.rp = adjoints_nv1;
				D1 = ((expr_if*)e)->d0;
				D->next = ((expr_if*)e)->dT;
				while(D->next != D1)
					D = D->next;
				((expr_if*)e)->dTlast = D;
				D->next = ((expr_if*)e)->dF;
				((expr_if*)e)->next = asl->P.iflist;
				asl->P.iflist = (expr_if*)e;
				}
			break;

		case Hv_plterm:
			e->dO.r = e->R.e->dO.r;
			break;

		case Hv_sumlist:
			ep = e->R.ep;
			for(dO = 0; (e1 = *ep); ep++)
				dO += e1->dO.r;
			e->dO.r = dO;
			break;

		case Hv_func:
			e->dO.r = 1.;
			break;

		case Hv_negate:
			e->dO.r = e->L.e->dO.r;
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
			e->dO.r = e->R.e->dO.r;
			break;

		case Hv_minusLR:
			e->dO.r = e->L.e->dO.r + e->R.e->dO.r;
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
	real aO, adO, t;

	aO = f->aO;
	adO = f->adO;
	da = f->da;
	for(dae = f->dae; da < dae; da++) {
		e = da->e;
		e->aO += aO;
		e->adO += adO;
		t = adO*e->dO.r;
		for(da1 = f->da; da1 < dae; da1++) {
			e = da1->e;
			e->aO += t;
			}
		}
	}

 static void
hv_back(register expr *e)
{
	register expr *e1, **ep, *e2;
	real adO, t1, t2;
	de *d;

	if (!e || (!e->aO && !e->adO))
		return;
	for(; e; e = e->bak)
	    switch(e->a) {
		case Hv_binaryR:
			e1 = e->R.e;
			e1->adO += e->adO;
			e1->aO += e->aO  +  e->adO * e1->dO.r;
			break;

		case Hv_binaryLR:
			e1 = e->L.e;
			e2 = e->R.e;
			adO = e->adO;
			t1 = adO * e1->dO.r;
			t2 = adO * e2->dO.r;
			e1->aO  += e->aO + t1 + t2;
			e2->aO  += e->aO + t1 + t2;
			e1->adO += adO;
			e2->adO += adO;
			break;

		case Hv_unary:
			e1 = e->L.e;
			e1->adO += e->adO;
			e1->aO  += e->aO  +  e->adO * e1->dO.r;
			break;

		case Hv_vararg:
			for(d = ((expr_va*)e)->L.d; d->e; d++) {
				if ((e1 = d->ee)) {
					e1->aO = e->aO;
					e1->adO = e->adO;
					hv_back(e1);
					}
				else {
					e1 = d->e;
					if (e1->op != f_OPNUM) {
						e1->aO = e->aO;
						e1->adO = e->adO;
						}
					}
				}
			break;

		case Hv_if:
			if ((e1 = ((expr_if *)e)->Te)) {
				e1->aO = e->aO;
				e1->adO = e->adO;
				hv_back(e1);
				}
			else {
				e1 = ((expr_if *)e)->T;
				if (e1->op != f_OPNUM) {
					e1->aO = e->aO;
					e1->adO = e->adO;
					}
				}
			if ((e1 = ((expr_if *)e)->Fe)) {
				e1->aO = e->aO;
				e1->adO = e->adO;
				hv_back(e1);
				}
			else {
				e1 = ((expr_if *)e)->F;
				if (e1->op != f_OPNUM) {
					e1->aO = e->aO;
					e1->adO = e->adO;
					}
				}
			break;

		case Hv_plterm:
			e->R.e->aO += e->aO;
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
			e1->aO += e->aO;
			e1->adO += e->adO;
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
			e1->aO += t1;
			e1->adO += t2;
			break;

		case Hv_timesR:
			e1 = e->R.e;
			e1->adO += e->adO;
			e1->aO += e->aO;
			break;

		case Hv_timesL:
			e1 = e->L.e;
			e1->adO += e->adO;
			e1->aO  += e->aO;
			break;

		case Hv_timesLR:
			e1 = e->L.e;
			e2 = e->R.e;
			adO = e->adO;
			e1->aO  += e->aO  +  adO * e2->dO.r;
			e2->aO  += e->aO  +  adO * e1->dO.r;
			e1->adO += adO;
			e2->adO += adO;
			break;

		default:/*DEBUG*/
			fprintf(Stderr, "bad e->a = %d in hv_back\n", e->a);
			exit(1);
		}
	}

 static void
hv_fwd0(ASL_pfgh *asl, register cexp *c, register expr_v *v)
{
	register linpart *L, *Le;
	real x;

	v->aO = v->adO = 0;
	if (c->ef) {
		hv_fwd(asl, c->ef);
		x = c->ee->dO.r;
		}
	else if (c->e->op != f_OPNUM)
		x = c->e->dO.r;
	else
		x = 0;
	if ((L = c->L))
		for(Le = L + c->nlin; L < Le; L++)
			x += ((expr_v*)L->v.vp)->dO.r;
	v->dO.r = x;
	}

 static void
pshv_prod1(ASL_pfgh *asl, range *r, int nobj, int ow, int y)
{
	int *cei, *cei0, *ceie, i;
	linarg *la, **lap, **lap1, **lape;
	linpart *L, *Le;
	expr_v *v;
	cexp *c;
	expr *e;
	real *s;
	psb_elem *b;

	s = asl->P.dOscratch;
	lap = lap1 = r->lap;
	lape = lap + r->n;
	while(lap1 < lape) {
		la = *lap1++;
		v = la->v;
		v->dO.r = *s++;
		v->adO = v->aO = 0.;
		}
	if ((cei = cei0 = r->cei)) {
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
			if (!ow && i != nobj)
				continue;
			}
		else if (!y)
			continue;
		if ((e = b->D.ef)) {
			hv_fwd(asl, e);
			e = b->D.ee;
			e->aO = 0;
			e->adO = 1.;
			hv_back(e);
			}
		else if ((e = b->D.e)->op != f_OPNUM) {
			e->aO = 0;
			e->adO = 1.;
			}
		}
	while(cei > cei0) {
		i = *--cei;
		c = cexps + i;
		v = asl->P.vp[i];
		if (v->aO && (L = c->L))
		    for(Le = L + c->nlin; L < Le; L++)
			((expr_v*)L->v.vp)->aO++;
		if ((e = c->ee)) {
			e->aO = 1.;
			e->adO = v->adO;
			hv_back(e);
			}
		else if ((e = c->e)->op != f_OPNUM) {
			e->aO = 1.;
			e->adO = v->adO;
			}
		}
	}

#ifdef __cplusplus
extern "C" {
static int compar(const void*, const void*);
}
#endif

 static int
compar(const void *a, const void *b)
{ return *(int*)a - *(int*)b; }

#undef nzc
#undef asl
#undef del_mblk
#define del_mblk(b,c) Del_mblk_ASL(a,b,(Char*)(c))

 static void
new_Hesoprod(ASL_pfgh *asl, ograd *L, ograd *R, real coef)
{
	Hesoprod *h, **hp, *h1, *h2;
	int kh;
	Char **mblk_free;

	ACQUIRE_DTOA_LOCK(HESOPROD_LOCK);
	if (!(h = asl->P.hop_free)) {
		mblk_free = asl->mblk_free;
		kh = asl->P.khesoprod;
		while(kh < 8 && !mblk_free[kh])
			kh++;
		asl->P.khesoprod = kh;
		h = h1 = (Hesoprod *)new_mblk(kh);
		h2 = h + (sizeof(Char*) << kh)/sizeof(Hesoprod) - 1;
		while(h1 < h2)
			h1 = h1->next = h1 + 1;
		h1->next = 0;
		}
	asl->P.hop_free = h->next;
	FREE_DTOA_LOCK(HESOPROD_LOCK);
	h->left = L;
	h->right = R;
	h->coef = coef;
	hp = asl->P.otodo + R->varno;
	h->next = *hp;
	*hp = h;
	}

 static void
del_Hesoprod(ASL_pfgh *asl, Hesoprod *x)
{
	x->next = asl->P.hop_free;
	asl->P.hop_free = x;
	}

 static real *
saveog(ASL_pfgh *asl, int no, int noe, int y, int *kp)
{
	real *o, *ogsave;
	int i, k, n;
	ps_func *p, *pe;
	psg_elem *g, *ge;
	ograd *og;

	n = 0;
	if (asl->P.nobjgroups)
		for(i = no; i < noe; i++) {
			p = asl->P.ops + i;
			g = p->g;
			for(ge = g + p->ng; g < ge; g++)
				for(og = g->og; og; og = og->next)
					n++;
			}
	if (asl->P.ncongroups && y) {
		p = asl->P.cps;
		for(pe = p + nlc; p < pe; p++)
			for(g = p->g, ge = g + p->ng; g < ge; g++)
				for(og = g->og; og; og = og->next)
					n++;
		}
	if (!n)
		return 0;
	k = *kp = htcl(n*sizeof(real));
	o = ogsave = (real*)new_mblk(k);
	if (asl->P.nobjgroups)
		for(i = no; i < noe; i++) {
			p = asl->P.ops + i;
			g = p->g;
			for(ge = g + p->ng; g < ge; g++)
				for(og = g->og; og; og = og->next)
					*o++ = og->coef;
			}
	if (asl->P.ncongroups && y) {
		p = asl->P.cps;
		for(pe = p + nlc; p < pe; p++)
			for(g = p->g, ge = g + p->ng; g < ge; g++)
				for(og = g->og; og; og = og->next)
					*o++ = og->coef;
		}
	return ogsave;
	}

 static void
restog(ASL_pfgh *asl, real *ogsave, int no, int noe, int y, int k)
{
	real *o = ogsave;
	int i;
	ps_func *p, *pe;
	psg_elem *g, *ge;
	ograd *og;

	if (asl->P.nobjgroups)
		for(i = no; i < noe; i++) {
			p = asl->P.ops + i;
			g = p->g;
			for(ge = g + p->ng; g < ge; g++)
				for(og = g->og; og; og = og->next)
					og->coef = *o++;
			}
	if (asl->P.ncongroups && y) {
		p = asl->P.cps;
		for(pe = p + nlc; p < pe; p++)
			for(g = p->g, ge = g + p->ng; g < ge; g++)
				for(og = g->og; og; og = og->next)
					og->coef = *o++;
		}
	Del_mblk_ASL((ASL*)asl, k, ogsave);
	}

 static fint
bothadj(ASL_pfgh *asl, SputInfo *spi)
{
	/* Adjust to compute both triangles of Hessian */
	fint *hr, *hre, *hrn, *hrn0;
	int kz;
	size_t *hcs, *ucs;
	size_t i, i0, i1, j, k, k0, L, n, n1, nz;
	ssize_t nod, *ulc, *uli,  *z, *z0, *z1;

	n = n_var;
	if ((nod = spi->nod) >= 0) {
		if (!nod)
			return 0;
		goto done;
		}
	n1 = n + 1;
	hcs = spi->hcolstartsZ;
	nod = nz = hcs[n] - hcs[0];
	hr = spi->hrownos - 1;
	i = i0 = Fortran;
	for(j = i + n; i < j; i++, hcs++) {
		hr += k = hcs[1] - hcs[0];
		if (k && *hr == i)
			--nod;
		}
	/* nod = number of off-diagonal elements in upper triangle */
	if (!(spi->nod = nod))
		return 0;	/* diagonal Hessian */
	nz += nod;
	L = nz*sizeof(fint) + (2*(nod+n1))*sizeof(size_t);
	if (sizeof(size_t) != sizeof(fint))
		L += n1*sizeof(fint);
	spi->khinfob = kz = htcl(L);
	spi->ulinc0 = uli = (ssize_t*)new_mblk(kz);
	spi->ulcopy0 = ulc = uli + n1;
	spi->hcs[1] = hcs = (size_t*)(ulc + 2*nod);
	hrn0 = (fint*)(hcs + n1);
	if (sizeof(size_t) != sizeof(fint))
		hrn0 += n1;
	spi->hrn[1] = hrn0;
	z = z0 = (ssize_t*)new_mblk(kz = htcl(n*sizeof(ssize_t)));
	z1 = z - Fortran;
	ucs = spi->hcs[0];
	hre = spi->hrn[0];
	for(i = i0; i < j; i++, ucs++) {
		hr = hre;
		hre += *z++ = ucs[1] - ucs[0];
		while(hr < hre)
			if ((k = *hr++) != i)
				z1[k]++;
		}
	ucs = spi->hcs[0];
	hre = spi->hrn[0];
	*uli++ = 0;
	for(i = k = i0; i < j; i++, ucs++) {
		hr = hre;
		hre += L = ucs[1] - ucs[0];
		*hcs++ = k;
		k0 = k - i0;
		hrn = hrn0 + k0;
		*uli++ = z1[i] - L;
		k += z1[i];
		z1[i] = k0 + L;
		while(hr < hre)
			if ((i1 = *hrn++ = *hr++) != i) {
				*ulc++ = k0++;
				hrn0[*ulc++ = z1[i1]++] = i;
				}
		}
	*hcs = k;
	spi->ulcend = ulc;
	Del_mblk_ASL((ASL*)asl, kz, z0);
	spi->ulinc = spi->ulinc0;
	spi->ulcopy = spi->ulcopy0;
 done:
	spi->hrownos = spi->hrn[1];
	spi->hcolstartsZ = spi->hcs[1];
	return nod;
	}

 static void
upper_to_lower(ASL_pfgh *asl, SputInfo *spi, size_t nz)
{	/* convert upper to lower triangular */

	fint *hrownos, *rn;
	int k, k1, *u0, *utoL;
	size_t L, *cs, *hcolstarts;
	ssize_t f, i, j, j1, j2, n, n1, *rs, *z;

	f = Fortran;
	n = n_var;
	n1 = n + 1;
	hrownos = spi->hrownos;
	hcolstarts = spi->hcolstartsZ;
	L = nz*sizeof(fint) + n1*sizeof(size_t);
	if (sizeof(size_t) != sizeof(fint))
		L += n1*sizeof(fint);
	spi->khinfob = k = htcl(L);
	spi->hcolstartsZ = cs = (size_t*)new_mblk(k);
	spi->ulinc0 = (ssize_t*)cs;
	rn = (fint*)(cs + n1);
	if (sizeof(size_t) != sizeof(fint))
		rn += n1;
	spi->hrownos = rn;
	k = htcl((n+nz)*sizeof(ssize_t));
	rs = (ssize_t*)new_mblk(k);
	z = rs + n;
	memset(rs, 0, n*sizeof(size_t));
	for(i = 0; i < nz; i++)
		rs[hrownos[i]-f]++;
	for(i = j = 0; i < n; i++) {
		cs[i] = j + f;
		j1 = rs[i];
		rs[i] = j;
		j += j1;
		}
	cs[n] = nz + f;
	j1 = hcolstarts[1] - f;
	for(i = j = 0; i < nz; i++) {
		while(i >= j1)
			j1 = hcolstarts[++j + 1] - f;
		rn[z[i] = rs[hrownos[i]-f]++] = j + f;
		}
	for(i = j = 0; i < nz; i++) {
		if ((j1 = z[i]) <= i) {
			if (j1 < 0)
				z[i] = -(j1 + 1);
			continue;
			}
		j += 3;
		while((j2 = z[j1]) != i) {
			z[j1] = -(j2 + 1);
			j++;
			j1 = j2;
			}
		}
	if (j) {
		j += 2;
		k1 = htcl(j*sizeof(int));
		spi->uptolow = utoL = (int*)new_mblk(k1);
		*utoL++ = k1;
		for(i = 0; i < nz; i++) {
			if ((j = z[i]) <= i)
				continue;
			u0 = utoL++;
			*utoL++ = i;
			*utoL++ = j;
			while((j2 = z[j]) != i) {
				z[j] = -(j2 + 1);
				j = *utoL++ = j2;
				}
			*u0 = (utoL - u0) - 1;
			}
		*utoL = 0;
		}
	Del_mblk_ASL((ASL*)asl, k, rs);
	}

 fint
sphes_setup_ASL(ASL *a, SputInfo **pspi, int nobj, int ow, int y, int uptri)
{
	ASL_pfgh *asl;
	Hesoprod *hop, *hop1, **otodo, **otodoi, **otodoj;
	SputInfo *spi, *spi1;
	de *d;
	derp *D1;
	expr_if *iflist;
	expr_v *v;
	expr_va *valist;
	fint *hr, *hre, *hrownos, rv;
	int *ui, *zc, *zci;
	int i, j, k, khinfo, kog, kz, n, n1, nhinfo, no, noe, nqslim, nzc;
	int rfilter, robjno;
	linarg *la, **lap, **lap1, **lape;
	ograd *og, *og1, **ogp, **ogpe;
	ps_func *p, *pe;
	psb_elem *b;
	psg_elem *g, *ge;
	range *r, *r0, **rp, **rtodo;
	real *ogsave, *s, *si, t;
	size_t *hcolstarts, iz, jz, n1spi, *tf;
	uHeswork *uhw, *uhwi, **utodo, **utodoi, **utodoj;

	asl = pscheck_ASL(a, "sphes_setup");
	n1 = n_var + 1;
	if (!pspi)
		pspi = &asl->i.sputinfo_;
	if (nobj >= 0 && nobj < n_obj) {
		robjno = -2 - nobj;
		rfilter = n_obj > 1 || (!y && nlc > 0);
		ow = 0;
		no = nobj;
		noe = no + 1;
		}
	else {
		robjno = nobj = -1;
		rfilter = (!ow && n_obj > 0) || (!y && nlc > 0);
		no = noe = 0;
		if (ow) {
			noe = n_obj;
			ow = 1;
			}
		}
	if (y)
		y = 1;
	if ((n = nlvo) < nlvc)
		n = nlvc;
	if ((spi = *pspi)) {
		if (spi->ow == ow && spi->y == y && spi->nobj == nobj
		 && spi->uptri == uptri)
			goto done;
		del_mblk(spi->khinfo, spi);
		if (spi->ulinc0)
			del_mblk(spi->khinfob, spi->ulinc0);
		if ((ui = spi->uptolow))
			del_mblk(*ui, ui);
		*pspi = 0;
		}
	if (!asl->P.hes_setup_called)
		(*asl->p.Hesset)(a, 1, 0, nlo, 0, nlc);
	asl->P.hes_setup_called = 3;
	asl->P.iflist = 0;
	asl->P.valist = 0;
	otodo = otodoi = asl->P.otodo;
	rtodo = asl->P.rtodo;
	utodo = utodoi = asl->P.utodo;
	s = asl->P.dOscratch;
	nqslim = n >> 3;
	kz = htcl(2*sizeof(int)*n);
	zc = (int*)new_mblk_ASL(a, kz);
	zci = zc + n;
	memset(zc, 0, n*sizeof(int));
	n1spi = sizeof(SputInfo) + n1*sizeof(size_t);
	if (sizeof(size_t) != sizeof(fint))
		n1spi += n1*sizeof(fint);
	khinfo = htcl((n1 + 30)*sizeof(fint) + n1spi);
	spi = (SputInfo*)new_mblk_ASL(a, khinfo);
	hcolstarts = (size_t*)(spi+1);
	hr = hrownos = (fint*)((char*)spi + n1spi);
	nhinfo = ((sizeof(Char*)<<khinfo) - n1spi) / sizeof(fint);
	hre = hr + nhinfo;
	r0 = (range*)&asl->P.rlist;
	for(r = asl->P.rlist.next; r != r0; r = r->rlist.next) {
		if ((j = r->n) <= 0)
			continue;
		if (rfilter) {
			for(b = r->refs; b; b = b->next) {
				if (b->conno >= 0) {
					if (y)
						goto keep;
					}
				else if (b->conno == robjno)
					goto keep;
				}
			continue;
			}
 keep:
		i = r->lasttermno;
		rp = rtodo + i;
		r->hnext = *rp;
		*rp = r;
		}
	kog = 0; /* silence bogus "not-initialized" warning */
	ogsave = asl->P.npsgcomp ? saveog(asl, no, noe, y, &kog) : 0;
	if (asl->P.nobjgroups)
		for(i = no; i < noe; i++) {
			p = asl->P.ops + i;
			g = p->g;
			for(ge = g + p->ng; g < ge; g++)
			    if ((og = g->og)) {
				do og->coef = 1; while((og = og->next));
				og = g->og;
				new_Hesoprod(asl, og, og, 1.);
				}
			}
	if (asl->P.ncongroups && y) {
		p = asl->P.cps;
		for(pe = p + nlc; p < pe; p++)
			for(g = p->g, ge = g + p->ng; g < ge; g++)
			    if ((og = g->og)) {
				do og->coef = 1; while((og = og->next));
				og = g->og;
				new_Hesoprod(asl, og, og, 1.);
				}
		}
	for(i = 0; i < n; i++) {
		nzc = 0;
		rp = rtodo;
		uhwi = *utodoi;
		*utodoi++ = 0;
		while((r = *rp)) {
			rp = &r->hnext;
			lap = r->lap;
			lape = lap + r->n;
			if (r->n >= r->nv) {
				k = htcl(sizeof(uHeswork)
					+ (r->n - 1)*sizeof(ograd*));
				uhw = (uHeswork *)new_mblk_ASL(a, k);
				uhw->k = k;
				uhw->next = uhwi;
				uhwi = uhw;
				uhw->r = r;
				uhw->ui = ui = r->ui;
				uhw->uie = ui + r->nv;
				ogp = uhw->ogp;
				while(lap < lape)
					*ogp++ = (*lap++)->nz;
				}
			else {
				si = s;
				while(lap < lape) {
					*si = 1;
					pshv_prod1(asl, r, nobj, ow, y);
					*si++ = 0;
					lap1 = lap++;
					la = *lap1++;
					og = la->nz;
					v = la->v;
					if ((t = v->aO))
						new_Hesoprod(asl,og,og,t);
					while(lap1 < lape) {
					    la = *lap1++;
					    v = la->v;
					    if ((t = v->aO)) {
						og1 = la->nz;
						new_Hesoprod(asl,og,og1,t);
						new_Hesoprod(asl,og1,og,t);
						}
					    }
					}
				}
			}
		*rtodo++ = 0;	/* reset */
		while((uhw = uhwi)) {
			uhwi = uhwi->next;
			si = s;
			ogp = uhw->ogp;
			r = uhw->r;
			ogpe = ogp + r->n;
			si = s;
			do {
				if ((og = *ogp) && og->varno == i)
					*si = 1.; /* og->coef til 20080629 */
				si++;
				} while(++ogp < ogpe);
			pshv_prod1(asl, r, nobj, ow, y);

			lap = r->lap;
			lape = lap + r->n;
			do {
				la = *lap++;
				if (la->v->aO)
					for(og = la->nz; og; og = og->next)
						if ((j = og->varno) <= i
						 && !zc[j]++)
							zci[nzc++] = j;
				}
				while(lap < lape);

			ogp = uhw->ogp;
			si = s;
			do {
				if ((og = *ogp) && og->varno == i) {
					*si = 0;
					*ogp = og->next;
					}
				si++;
				} while(++ogp < ogpe);
			if ((ui = ++uhw->ui) >= uhw->uie)
				del_mblk(uhw->k, uhw);
			else {
				utodoj = utodo + *ui;
				uhw->next = *utodoj;
				*utodoj = uhw;
				}
			}

		hop1 = *otodoi;
		*otodoi++ = 0;
		while((hop = hop1)) {
			hop1 = hop->next;
			og = hop->left;
			og1 = hop->right;
			while((j = og->varno) <= i) {
				if (!zc[j]++)
					zci[nzc++] = j;
				if (!(og = og->next))
					break;
				}
			if ((og = og1->next)) {
				hop->right = og;
				otodoj = otodo + og->varno;
				hop->next = *otodoj;
				*otodoj = hop;
				}
			else
				del_Hesoprod(asl,hop);
			}
		hcolstarts[i] = hr - hrownos;
		if (nzc > hre - hr) {
			k = khinfo++;
			spi1 = (SputInfo*)new_mblk_ASL(a, khinfo);
			tf = (size_t*)(spi1+1);
			memcpy(tf, hcolstarts, (char*)hr - (char*)hcolstarts);
			del_mblk(k, spi);
			spi = spi1;
			hcolstarts = tf;
			hrownos = (fint*)((char*)spi1 + n1spi);
			hr = hrownos + hcolstarts[i];
			nhinfo = ((sizeof(Char*)<<khinfo) - n1spi) / sizeof(fint);
			hre = hrownos + nhinfo;
			}
		if (nzc > nqslim) {
			for(j = 0; j < n; j++)
				if (zc[j])
					zc[*hr++ = j] = 0;
			}
		else {
			if (nzc > 1)
				qsort(zci, nzc, sizeof(int), compar);
			for(j = 0; j < nzc; j++)
				zc[*hr++ = zci[j]] = 0;
			}
		}
	for(valist = asl->P.valist; valist; valist = valist->next) {
		D1 = valist->d0;
		for(d = valist->L.d; d->e; d++)
			d->dlast->next = D1;
		}
	for(iflist = asl->P.iflist; iflist; iflist = iflist->next)
		iflist->dTlast->next = iflist->d0;
	jz = hcolstarts[n] = hr - hrownos;
	for(i = n; ++i < n1; )
		hcolstarts[i] = jz;
	if ((j = Fortran)) {
		for(i = 0; i < n1; i++)
			hcolstarts[i] += j;
		iz = hcolstarts[n] - j;
		while(iz)
			hrownos[--iz] += j;
		}
	spi->hcs[0] = hcolstarts;
	spi->hrn[0] = hrownos;
	spi->nod = -1;
	spi->ulcend = 0;
	spi->khinfo = khinfo;
	spi->nobj = nobj;
	spi->ow = ow;
	spi->y = y;
	spi->uptri = uptri;
	*pspi = spi;
	if (ogsave)
		restog(asl, ogsave, no, noe, y, kog);
	spi->ulinc0 = spi->ulinc = 0;
	spi->ulcopy = 0;
	spi->uptolow = 0;
	del_mblk(kz, zc);
 done:
	spi->hrownos = spi->hrn[0];
	spi->hcolstartsZ = hcolstarts = spi->hcs[0];
	rv = hcolstarts[n] - hcolstarts[0];
	switch(uptri) {
	  case 0:
		rv += bothadj(asl, spi);
		break;
	  case 2:
		upper_to_lower(asl, spi, rv);
	  }
	hcolstarts = spi->hcolstartsZ;
	if (sizeof(size_t) == sizeof(fint))
		spi->hcolstarts = (fint*)hcolstarts;
	else {
		spi->hcolstarts = hr = (fint*)(hcolstarts + n1);
		for(i = 0; i < n1; ++i)
			hr[i] = hcolstarts[i];
		}
	return rv;
	}

 void
sphes_ASL(ASL *a, SputInfo **pspi, real *H, int nobj, real *ow, real *y)
{
	/* sparse upper triangle of Hessian */

	ASL_pfgh *asl;
	Hesoprod *hop, *hop1, **otodo, **otodoi, **otodoj;
	SputInfo *spi;
	expr_v *v;
	fint *hr;
	int i, j, k, kh, n, no, noe, *ui;
	linarg *la, **lap, **lap1, **lape;
	ograd *og, *og1, **ogp, **ogpe;
	ps_func *p, *pe;
	psg_elem *g, *ge;
	range *r, *r0, **rp, **rtodo;
	real *Hi, *H0, *H00;
	real *cscale, *owi, *s, *si, t, t1, *vsc0, *vsc1, *vsc, *y1;
	size_t *hcs;
	ssize_t *ulc, *uli;
	uHeswork *uhw, *uhwi, **utodo, **utodoi, **utodoj;

	asl = pscheck_ASL(a, "sputhes");
	xpsg_check_ASL(asl, nobj, ow, y);
	if (!pspi)
		pspi = &a->i.sputinfo_;
	i = j = 0;
	if (y)
		j = 1;
	if (nobj >= 0 && nobj < n_obj) {
		no = nobj;
		noe = no + 1;
		owi = ow ? ow + no : &edag_one_ASL;
		ow = 0;
		}
	else {
		nobj = -1;
		no = noe = 0;
		if ((owi = ow)) {
			noe = n_obj;
			i = 1;
			}
		}
	if (asl->P.hes_setup_called != 3)
		sphes_setup_ASL(a, pspi, nobj, ow != 0, y != 0, 0);
	spi = *pspi;
	if (spi->nobj != nobj || spi->ow < i || spi->y < j) {
		fprintf(Stderr,
		 "\nsphes() call inconsistent with previous sphsetup()\n");
		exit(1);
		}
	otodo = otodoi = asl->P.otodo;
	rtodo = asl->P.rtodo;
	utodo = utodoi = asl->P.utodo;
	s = asl->P.dOscratch;
	if ((n = nlvo) < nlvc)
		n = nlvc;
	Hi = H0 = (real*)new_mblk_ASL(a, kh = htcl(n*sizeof(real)));
	memset(Hi, 0, n * sizeof(real));
	H0 -= Fortran;
	r0 = (range*)&asl->P.rlist;
	for(r = asl->P.rlist.next; r != r0; r = r->rlist.next) {
		if ((j = r->n) <= 0)
			continue;
		i = r->lasttermno;
		rp = rtodo + i;
		r->hnext = *rp;
		*rp = r;
		}
	if (asl->P.nobjgroups)
	    for(; no < noe; no++)
		if ((t = *owi++)) {
		    p = asl->P.ops + no;
		    g = p->g;
		    for(ge = g + p->ng; g < ge; g++)
			if ((t1 = t*g->g2))
				for(og = g->og; og; og = og->next)
					if (og->coef) {
						new_Hesoprod(asl, og, og, t1);
						break;
						}
		}
	if (asl->P.ncongroups && y) {
		cscale = asl->i.lscale;
		p = asl->P.cps;
		y1 = y;
		for(pe = p + nlc; p < pe; p++, y1++)
			if ((t = cscale ? *cscale++ * *y1 : *y1))
				for(g = p->g, ge = g + p->ng; g < ge; g++)
				    if ((t1 = t*g->g2))
					for(og = g->og; og; og = og->next)
					    if (og->coef) {
						new_Hesoprod(asl, og, og, t1);
						break;
						}
		}
	hcs = spi->hcs[0];
	hr  = spi->hrn[0];
	uli = spi->ulinc;
	H00 = H;
	vsc = asl->i.vscale;
	vsc0 = vsc - Fortran;
	vsc1 = vsc;
	for(i = 0; i < n; i++) {
		rp = rtodo;
		uhwi = *utodoi;
		*utodoi++ = 0;
		while((r = *rp)) {
			rp = &r->hnext;
			lap = r->lap;
			lape = lap + r->n;
			if (r->n >= r->nv) {
				k = htcl(sizeof(uHeswork)
					+ (r->n - 1)*sizeof(ograd*));
				uhw = (uHeswork *)new_mblk_ASL(a, k);
				uhw->k = k;
				uhw->next = uhwi;
				uhwi = uhw;
				uhw->r = r;
				uhw->ui = ui = r->ui;
				uhw->uie = ui + r->nv;
				ogp = uhw->ogp;
				while(lap < lape)
					*ogp++ = (*lap++)->nz;
				}
			else {
				si = s;
				while(lap < lape) {
					*si = 1;
					pshv_prod_ASL(asl, r, nobj, ow, y);
					*si++ = 0;
					lap1 = lap++;
					la = *lap1++;
					og = la->nz;
					v = la->v;
					if ((t = v->aO))
						new_Hesoprod(asl,og,og,t);
					while(lap1 < lape) {
					    la = *lap1++;
					    v = la->v;
					    if ((t = v->aO)) {
						og1 = la->nz;
						new_Hesoprod(asl,og,og1,t);
						new_Hesoprod(asl,og1,og,t);
						}
					    }
					}
				}
			}
		*rtodo++ = 0;	/* reset */
		while((uhw = uhwi)) {
			uhwi = uhwi->next;
			si = s;
			ogp = uhw->ogp;
			r = uhw->r;
			ogpe = ogp + r->n;
			si = s;
			do {
				if ((og = *ogp) && og->varno == i)
					*si = og->coef;
				si++;
				} while(++ogp < ogpe);
			pshv_prod_ASL(asl, r, nobj, ow, y);

			lap = r->lap;
			lape = lap + r->n;
			do {
				la = *lap++;
				if ((t = la->v->aO))
					for(og = la->nz; og; og = og->next)
						if ((j = og->varno) <= i)
							Hi[j] += t*og->coef;
				}
				while(lap < lape);

			ogp = uhw->ogp;
			si = s;
			do {
				if ((og = *ogp) && og->varno == i) {
					*si = 0;
					*ogp = og->next;
					}
				si++;
				} while(++ogp < ogpe);
			if ((ui = ++uhw->ui) >= uhw->uie)
				del_mblk(uhw->k, uhw);
			else {
				utodoj = utodo + *ui;
				uhw->next = *utodoj;
				*utodoj = uhw;
				}
			}

		hop1 = *otodoi;
		*otodoi++ = 0;
		while((hop = hop1)) {
			hop1 = hop->next;
			og = hop->left;
			og1 = hop->right;
			t = hop->coef * og1->coef;
			while((j = og->varno) <= i) {
				Hi[j] += t*og->coef;
				if (!(og = og->next))
					break;
				}
			if ((og = og1->next)) {
				hop->right = og;
				otodoj = otodo + og->varno;
				hop->next = *otodoj;
				*otodoj = hop;
				}
			else
				del_Hesoprod(asl,hop);
			}
		k = (int)(hcs[1] - hcs[0]);
		hcs++;
		if (uli)
			H += *uli++;
		if (vsc) {
			t = *vsc1++;
			while(--k >= 0) {
				j = (int)*hr++;
				*H++ = t * vsc0[j] * H0[j];
				H0[j] = 0;
				}
			}
		else
			while(--k >= 0) {
				*H++ = H0[j = (int)*hr++];
				H0[j] = 0;
				}
		}
	del_mblk(kh, Hi);
	H = H00;
	if ((ulc = spi->ulcopy))
		for(uli = spi->ulcend; ulc < uli; ulc += 2)
			H[ulc[1]] = H[ulc[0]];
	else if ((ui = spi->uptolow))
		while((k = *++ui)) {
			t = H[j = *++ui];
			while(--k) {
				t1 = H[i = *++ui];
				H[i] = t;
				t = t1;
				}
			H[j] = t;
			}
	}
