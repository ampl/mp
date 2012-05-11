/****************************************************************
Copyright (C) 1997-1999, 2001 Lucent Technologies
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
#include "getstub.h"

 static ASL_pfgh *asl;

#ifdef __cplusplus
 extern "C" {
#endif

extern struct dcomsb_1_ {
    real acccg, ratio, radius, radmax, findmx, cma31;
    fint itercg, itcgmx, ngeval, iskip, ifixed, nsemib;
} dcomsb_;

extern struct dcomal_1_ {
    real rmu, rmutol, firstc, firstg;
    fint newsol;
} dcomal_;

 typedef void (*range_func) ANSI((fint*, fint*, real*, real*, fint*, fint*));

 extern void dauglg_ ANSI((
	fint *n, fint *ng, fint *nel, fint *ieling, fint *leling,
	 fint *istadg,
	fint *lstadg, fint *ielvar, fint *lelvar, fint *istaev,
	 fint *lstaev, fint *intvar, fint *lntvar,
	fint *istadh, fint *lstadh, fint *icna, fint *licna,
	 fint *istada, fint *lstada,
	real *a, fint *la, real *b, fint *lb, real *bl, fint *lbl,
	 real *bu, fint *lbu,
	real *gscale, fint *lgscal, real *escale, fint *lescal,
	 real *vscale, fint *lvscal,
	fint *gxeqx, fint *lgxeqx, fint *intrep, fint *lintre,
	 fint *kndofc, fint *lkndof,
	range_func, fint *inform, real *fobj, real *x, fint *lx,
	 real *u, fint *lu, real *gvals,
	fint *lgvals, real *ft, fint *lft, real *fuvals,
	 fint *lfuval, real *xt,
	fint *lxt, fint *icalcf, fint *lcalcf, fint *ncalcf,
	 fint *icalcg, fint *lcalcg,
	fint *ncalcg, fint *ivar, fint *livar, fint *nvar, real *q,
	 fint *lq, real *dgrad,
	fint *ldgrad, fint *ichose, fint *iter, fint *maxit,
	 fint *quadrt, char *vnames,
	fint *lvname, char *gnames, fint *lgname, real *stopg,
	 real *stopc, fint *iwk, fint *liwk, real *wk, fint *lwk,
	fint *iprint, fint *iout, fint vnames_len, fint gnames_len));

 static fint ICHOSE[6] = { 2, 8, 0, 0, 1, 2 };
 static fint LIWK, LWK, OUTLEV;
 static fint MAXIT = 1000;
 static real STOPG = 1e-5;
 static real STOPC = 1e-5;
 static fint lincmax =   0x800000;
 static fint nprob = 1;
 static fint pfreq = 1, pstart = 0, pstop = 1000000;
 static fint wantgr = 1;
 static fint solmsg = 1;
 static int cl_timing, nmax, timing;
 static int *Cmap;
 static psb_elem **Emap;
 static psg_elem **Gmap;
 static fint *IELVAR, *INTREP, *INTVAR, *ISTAEV, *Vused;
 static fint *Iscr, *Iscr1, **Rmap, *Rmapnext;
 static fint LICNA, NEL, NFCALL, NGCALL, NINTV, NINTV2,
	lwinc, maxsel, maxsin, nRmap, nvargp;
 static range Range0;
#ifdef WANT_IEEE0	/* for debugging */
 static fint call_ieee0 = 0;
#endif
 static real objscale;
 static real *Yscrch;
 typedef struct
Hdispatch {
	struct Hdispatch *next;
	ograd *og;
	int ivarno;
	} Hdispatch;
 static Hdispatch **hdisp, *hdfree;
 static real *Hscratch;

 static keyword
keywds[] = {	/* must be in alphabetical order */
 KW("cauchy", L_val, &ICHOSE[4], "1 = exact Cauchy point, 0 = approx."),
 KW("ctol", D_val, &STOPC, "tol. for maxnorm(constraint violation); default 1e-5"),
#if 0
 KW("fdiv", L_val, &ICHOSE[2], "gradients: 0 = analytic; 1 = finite diffs"),
#endif
 KW("firstc", D_val, &dcomal_.firstc, "constraint tol. for multiplier estimates; default .1"),
 KW("firstg", D_val, &dcomal_.firstg, "projected grad. tol. for multiplier estimates; default .1"),
 KW("groups", L_val, &wantgr, "1 = find group structure"),
 KW("gtol", D_val, &STOPG, "tol. for maxnorm(Lagrangian gradient); default 1e-5"),
 KW("imeth", L_val, &ICHOSE[1], "method for solving linear system (default 8):\n\
		1 = CG (no preconditioner)\n\
		2 = diagonal precond. CG\n\
		4 = expanding-band precond. CG\n\
		5 = Munksgaard's precond. CG\n\
		6 = Schnabel-Eskow modified Cholesky precond. CG\n\
		7 = Gill-Murray-Ponceleon-Saunders mod. Chol. precond. CG\n\
		8 = band matrix precond. CG (default)\n\
		11 = multifrontal direct\n\
		12 = modified multifront direct"),
#ifdef WANT_IEEE0
 KW("infnan", L_val, &call_ieee0, "1 = trap creation of NaN and Infinity"),
#endif
 KW("lincmax", L_val, &lincmax, "max increment to liwk and lwk"),
 KW("liwk", L_val, &LIWK, "length of IWK"),
 KW("lwinc", L_val, &lwinc, "increment to default work-array lengths"),
 KW("lwk", L_val, &LWK, "length of WK"),
 KW("maxfwd", IA_val, voffset_of(ASL,p.maxfwd_), "# of partials to forward recur; default 5"),
 KW("maxit", L_val, &MAXIT, "maximum iterations; default 1000"),
 KW("mu", D_val, &dcomal_.rmu, "initial penalty parameter; default .1"),
 KW("mutol", D_val, &dcomal_.rmutol, "max. mu value for updating Lagrange multipliers; default .1"),
 KW("objno", L_val, &nprob, "objective number (1 = first)"),
 KW("objscale", D_val, &objscale, "scale factor for objective"),
 KW("outlev", L_val, &OUTLEV, "printing: 0 = none; 1-6 = successively more"),
 KW("pfreq", L_val, &pfreq, "print every pfreq iterations (default 1)"),
 KW("pstart", L_val, &pstart, "iteration to start printing (default 0)"),
 KW("pstop", L_val, &pstop, "iteration to stop printing (default 1000000)"),
 KW("qmin", L_val, &ICHOSE[5], "2 = min. quad. model loosely; 1 = tightly"),
 KW("radius", D_val, &dcomsb_.radius, "initial trust radius; default -1"),
 KW("sbw", L_val, &dcomsb_.nsemib, "semi-bandwith for preconditioners; default 5"),
 KW("sdkind", L_val, &ICHOSE[3], "Hessians: 0 = anal.; 1 = BFGS, 2 = DFP, 3 = PSB, 4 = SR1"),
 KW("solmsg", L_val, &solmsg, "1 = echo solution msg, 0 = do not"),
 KW("timing",  I_val, &timing,  "report I/O and solution times: 1 = stdout, 2 = stderr, 3 = both"),
 KW("trtype", L_val, &ICHOSE[0], "trust-region type: 2 = box, 1 = spherical"),
 KW("version", Ver_val, 0, "report version"),
 KW("wantsol", WS_val, 0, WS_desc_ASL+5)
 };

 static keyword cl_options[] = {
	KW("t", IK1_val, &cl_timing, "report solve time, eval time, etc.")
	};

 static char xxxvers[] = "AMPL/LANCELOT-A\0\nAMPL/LANCELOT-A Driver Version 20020506\n";

 static Option_Info Oinfo = { "lancelot", "LANCELOT", "lancelot_options",
				keywds, nkeywds, 1, xxxvers, 0,0,0,
				cl_options, 1, 20020506 };

 static void
#ifdef KR_headers
ranges(IELEM, TRANSP, W1, W2, NELV, NINV)
	fint *IELEM, *TRANSP, *NELV, *NINV; real *W1, *W2;
#else
ranges(fint *IELEM, fint *TRANSP, real *W1, real *W2, fint *NELV, fint *NINV)
#endif
{
	fint i = *IELEM - 1;
	psb_elem *b = Emap[i];
	range *r = b->U;
	fint *rm = Rmap[i];
	linarg **L, **Le;
	ograd *og;
	real t;

	L = r->lap;
	Le = L + r->n;
	if (*TRANSP) {
		memset(W2, 0, r->nv*sizeof(real));
		while(L < Le) {
			t = *W1++;
			for(og = (*L++)->nz; og; og = og->next)
				W2[*rm++] += t*og->coef;
			}
		}
	else {
		while(L < Le) {
			t = 0;
			for(og = (*L++)->nz; og; og = og->next)
				t += W1[*rm++]*og->coef;
			*W2++ = t;
			}
		}
	}

 static fint
#ifdef KR_headers
nvsum(p, cgp) ps_func *p; cgrad **cgp;
#else
nvsum(ps_func *p, cgrad **cgp)
#endif
{
	cgrad *cg;
	fint *z, *zc, *zci;
	fint licna, nel, nintv, nintv2, nrmap, nrmap0, nv;
	int i;
	linarg **L, **Le;
	linpart *lp, *lpe;
	ograd *og;
	psb_elem *b, *be;
	psg_elem *g, *ge;
	range *r;

	z = Iscr;
	zci = zc = Iscr1;
	licna = nel = nintv = nintv2 = nrmap = nv = 0;
	if (cgp) {
		for(cg = *cgp; cg; cg = cg->next)
			z[*zci++ = cg->varno] = 1;
		licna += zci - zc;
		}
	for(b = p->b, be = b + p->nb; b < be; b++) {
		if (!(r = b->U)) {
			if (p->nb == 1 && p->b->D.e->op == f_OPNUM)
				continue;
			r = &Range0;
			}
		nel++;
		nv += r->nv;
		if (maxsel < r->nv)
			maxsel = r->nv;
		L = r->lap;
		Le = L + r->n;
		nrmap0 = nrmap;
		while(L < Le)
			for(og = (*L++)->nz; og; og = og->next) {
				nrmap++;
				if (!z[i = og->varno]++)
					*zci++ = i;
				}
		if (!cgp) {
			nvargp += zci - zc;
			while(zci > zc)
				z[*--zci] = 0;
			}
		if ((i = r->n) >= r->nv) {
			i = r->nv;
			nrmap = nrmap0;
			}
		if (maxsin < i)
			maxsin = i;
		nintv += i;
		nintv2 += i*(i+1) >> 1;
		}
	if (zci > zc) {
		nvargp += zci - zc;
		while(zci > zc)
			z[*--zci] = 0;
		}
	else
	     for(g = p->g, ge = g + p->ng; g < ge; g++) {
		nel += g->ns;
		if (i = g->nlin) {
			licna += i;
			lp = g->L;
			lpe = lp + i;
			do {
				if (!z[i = lp->v.i]++)
					*zci++ = i;
				}
				while(++lp < lpe);
			}
		for(b = g->E, be = b + g->ns; b < be; b++) {
			if (!(r = b->U))
				r = &Range0;
			nv += r->nv;
			if (maxsel < r->nv)
				maxsel = r->nv;
			nrmap0 = nrmap;
			L = r->lap;
			Le = L + r->n;
			while(L < Le)
				for(og = (*L++)->nz; og; og = og->next) {
					nrmap++;
					if (!z[i = og->varno]++)
						*zci++ = i;
					}
			if ((i = r->n) >= r->nv) {
				i = r->nv;
				nrmap = nrmap0;
				}
			if (maxsin < i)
				maxsin = i;
			nintv += i;
			nintv2 += i*(i+1) >> 1;
			}
		nvargp += zci - zc;
		while(zci > zc)
			z[*--zci] = 0;
		}
	LICNA += licna;
	NEL += nel;
	NINTV += nintv;
	NINTV2 += nintv2;
	nRmap += nrmap;
	return nv;
	}

 static int
#ifdef KR_headers
fintcomp(a, b) char *a; char *b;
#else
fintcomp(const void *a, const void *b)
#endif
{ return (int)(*(fint *)a - *(fint *)b); }

 static fint *
#ifdef KR_headers
elvars(Iev, r, j) fint *Iev; range *r; int j;
#else
elvars(fint *Iev, range *r, int j)
#endif
{
	linarg **lap, **lape;
	fint *ri, *vused, *z, *zc, *zci, *zc1;
	ograd *og;
	int i;

	if (!r)
		r = &Range0;
	z = Iscr;
	zc = Iscr1;
	zci = zc;
	lap = r->lap;
	lape = lap + r->n;
	while(lap < lape)
		for(og = (*lap++)->nz; og; og = og->next)
			if (!z[i = og->varno]++)
				*zci++ = i;
	if (r->n < r->nv) {
		ri = Rmapnext;
		Rmap[j] = ri;
		INTREP[j] = 1;
		INTVAR[j] = r->n;
		i = 0;
		zc1 = zc;
		while(zc1 < zci)
			z[*zc1++] = i++;
		lap = r->lap;
		while(lap < lape)
			for(og = (*lap++)->nz; og; og = og->next)
				*ri++ = z[og->varno];
		Rmapnext = ri;
		}
	else {
		Rmap[j] = 0;
		INTREP[j] = 0;
		INTVAR[j] = r->nv;
		if (r->nintv) {
			qsort((char*)zc, zci - zc, sizeof(fint), fintcomp);
			if (nmax < r->n)
				nmax = r->n;
			}
		}
	vused = Vused;
	while(zc < zci) {
		z[i = *zc++] = 0;
		vused[i]++;
		*Iev++ = i + 1;
		}
	return Iev;
	}

 static void
#ifdef KR_headers
namegen(s, n, c) char *s; fint n; int c;
#else
namegen(char *s, fint n, int c)
#endif
{
	fint i;
	int p;
	p = 1;
	for(i = n; i >= 10; i /= 10)
		p++;
	for(i = 1; i <= n; i++, s += 10)
		sprintf(s, "%c%-09.*ld", c, p, i);
	}

 static void
#ifdef KR_headers
derivs(j, g, h) fint j; real *g, *h;
#else
derivs(fint j, real *g, real *h)
#endif
{
	psb_elem *b;
	range *r;
	linarg *la, **L, **L1, **Le, **Le2;
	int i, k, nobj;
	fint *I, *I0, *I1, *Ie, *rm, *z, *zc, *zci;
	real *Adj0, *Adjoints, *Hs0, *s, t, *y, *yzap;
	ograd *og;
	Hdispatch *hd, *hdnext, **hdp, **hdp0, **hdp1;

	Adjoints = adjoints;
	Adj0 = Adjoints - 1;
	b = Emap[j];
	if (!(r = b->U))
		r = &Range0;
	L = L1 = r->lap;
	Le = L + r->n;
	I = I0 = IELVAR + ISTAEV[j] - 1;
	Ie = I + ISTAEV[j+1] - ISTAEV[j];
	while(I < Ie)
		Adj0[*I++] = 0;
	if (i = r->nintv) {
		Le2 = L + i;
		do Adjoints[(*L++)->v->a] = 0.;
			while(L < Le2);
		}
	if (k = b->D.zaplen) {
		memset(adjoints_nv1, 0, k);
		derprop(b->D.d);
		}
	if ((k = Cmap[j]) < 0) {
		nobj = nprob;
		y = yzap = 0;
		}
	else {
		nobj = -1;
		y = Yscrch;
		yzap = y + k;
		*yzap = 1.;
		}
	if (rm = Rmap[j]) {
		L = L1;
		do *g++ = Adjoints[(*L++)->v->a];
			while(L < Le);
		}
	else {
		if (r->nintv) {
			L = L1;
			Le2 = L + r->nintv;
			do {
				if (t = Adjoints[(la = *L++)->v->a]) {
					og = la->nz;
					do Adjoints[og->varno] += t*og->coef;
						while(og = og->next);
					}
				} while(L < Le2);
			}
		I = I0;
		while(I < Ie)
			*g++ = Adj0[*I++];
		}
	if (ICHOSE[3])
		goto done;
	r->refs = b;
	b->next = 0;
	s = asl->P.dOscratch;
	if (rm || !r->nintv) {
		for(; L1 < Le; L1++) {
			*s = 1.;
			pshv_prod(r, nobj, 0, y);
			*s++ = 0.;
			L = r->lap;
			do *h++ = (*L++)->v->aO;
				while(L <= L1);
			}
		goto done;
		}
	/* The complicated case: we have internal variables */
	/* that LANCELOT does not see. */
	z = Iscr;
	zc = Iscr1;
	i = 0;
	do {
		hd = hdfree;
		hd->og = og = L1[0]->nz;
		hdfree = hd->next;
		hdp = hdisp + og->varno;
		hd->next = *hdp;
		*hdp = hd;
		hd->ivarno = i++;
		}
		while(++L1 < Le);
	I = I0;
	Hs0 = Hscratch - 1;
	hdp0 = hdisp - 1;
	zci = zc;
	do {
		k = (int)*I++;
		hdp = hdp0 + k;
		for(hd = *hdp; hd; hd = hd->next)
			s[hd->ivarno] = hd->og->coef;
		pshv_prod(r, nobj, 0, y);
		for(hd = *hdp; hd; hd = hdnext) {
			hdnext = hd->next;
			s[hd->ivarno] = 0.;
			if (og = hd->og->next) {
				hd->og = og;
				hdp1 = hdisp + og->varno;
				hd->next = *hdp1;
				*hdp1 = hd;
				}
			else {
				hd->next = hdfree;
				hdfree = hd;
				}
			}
		*hdp = 0;
		L = r->lap;
		do {
			la = *L++;
			og = la->nz;
			if (og->varno >= k)
				break;
			if (t = la->v->aO) do {
				Hscratch[i = og->varno] += t*og->coef;
				if (!z[i = og->varno]++)
					*zci++ = i;
				}
				while(og = og->next);
			}
			while(L < Le);
		I1 = I0;
		do *h++ = Hs0[*I1++];
			while(I1 < I);
		while(zci > zc) {
			z[i = *--zci] = 0;
			Hscratch[i] = 0.;
			}
		}
		while(I < Ie);
 done:
	if (yzap)
		*yzap = 0;
	}

 static void
#ifdef KR_headers
funcvals(ICALCF, n, FUVALS, XT) fint *ICALCF; int n; real *FUVALS, *XT;
#else
funcvals(fint *ICALCF, int n, real *FUVALS, real *XT)
#endif
{
	int i, j;
	expr *e;
	real t;

	errno = 0;
	xp_check_ASL(asl, XT);
	NFCALL++;
	for(i = 0; i < n; i++) {
		j = ICALCF[i] - 1;
		e = Emap[j]->D.e;
		FUVALS[j] = e->op(e C_ASL);
		}
	}

 static void
#ifdef KR_headers
funcgrads(ICALCF, n, FUVALS0, INTVAR, ISTADH)
	fint *ICALCF, *INTVAR, *ISTADH; int n; real *FUVALS0;
#else
funcgrads(fint *ICALCF, int n, real *FUVALS0, fint *INTVAR, fint *ISTADH)
#endif
{
	int i, j;

	NGCALL++;
	for(i = 0; i < n; i++) {
		j = ICALCF[i] - 1;
		derivs(j, FUVALS0 + INTVAR[j], FUVALS0 + ISTADH[j]);
		}
	}

 static void
#ifdef KR_headers
groupvals(ICALCG, n, ft, gv) fint *ICALCG; int n; real *ft, *gv;
#else
groupvals(fint *ICALCG, int n, real *ft, real *gv)
#endif
{
	psg_elem *g, **gm1;
	expr *e;
	int i, j;

	errno = 0;
	gm1 = Gmap - 1;
	--ft;
	--gv;
	for(i = 0; i < n; i++)
		if (g = gm1[j = ICALCG[i]]) {
			g->esum.v = ft[j];
			e = g->g;
			gv[j] = e->op(e C_ASL);
			}
		else
			gv[j] = ft[j];
	}

 static void
#ifdef KR_headers
groupgrads(ICALCG, n, g1, g2) fint *ICALCG; int n; real *g1, *g2;
#else
groupgrads(fint *ICALCG, int n, real *g1, real *g2)
#endif
{
	psg_elem *g, **gm1;
	expr *e, *ee;
	real t, t1, t2;
	int i, j, k;

	gm1 = Gmap - 1;
	--g1;
	--g2;
	for(i = 0; i < n; i++) {
		if (!(g = gm1[j = ICALCG[i]])) {
			g1[j] = 1.;
			g2[j] = 0.;
			continue;
			}
		/* from conpval.c */
		e = g->g;
		ee = g->ge;
		if (e == ee) {
			t = e->dL;
			t2 = e->dL2;
			}
		else {
			t = e->dL;
			do {
				e = e->L.e;
				t *= e->dL;
				}
				while(e != ee);
			e = g->g;
			if (t != 0.) {
				t1 = ee->dL;
				t2 = ee->dL2 * (t / t1);
				for(;;) {
					ee = ee->R.e;
					t2 += (t / ee->dL) * t1 * ee->dL2;
					if (ee == e)
						break;
					t1 *= ee->dL;
					}
				
				}
			else
				for(k = 0, t2 = 1.;; e = e->L.e) {
					if (e->dL)
						t2 *= e->dL;
					else if (k++) {
						t2 = 0;
						break;
						}
					else
						t2 *= e->dL2;
					if (e == ee)
						break;
					}
			}
		g1[j] = g->g1 = t;
		g2[j] = g->g2 = t2;
		}
	}

 static real *
#ifdef KR_headers
shiftup(n, x, y) int n; real *x, *y;
#else
shiftup(int n, real *x, real *y)
#endif
{
	real *z = y;
	x += n;
	y += n;
	while(y > z)
		*--y = *--x;
	return z;
	}

 static void
#ifdef KR_headers
bad_inform(INFORM) fint INFORM;
#else
bad_inform(fint INFORM)
#endif
{
	fprintf(Stderr, "Unexpected INFORM value: %ld\n", INFORM);
	exit(1);
	}

 static real
#ifdef KR_headers
constterm(p, c) ps_func *p; cde *c;
#else
constterm(ps_func *p, cde *c)
#endif
{
	expr_n *e;

	switch(p->nb) {
	  case 0:
		if (p->ng == 0)
			return -((expr_n*)c->e)->v;
		break;
	  case 1:
		if (!(p->b->U) && (e = (expr_n*)p->b->D.e)->op == f_OPNUM_ASL)
			return -e->v;
	  }
	return 0;
	}

 static void
#ifdef KR_headers
time_out(f, T) FILE *f; real *T;
#else
time_out(FILE *f, real *T)
#endif
{
	fprintf(f, "\n LANCELOT seconds:\n read: %10.2f\n", T[1] - T[0]);
	fprintf(f, " solve: %9.2f", T[2] - T[1]);
	fprintf(f, "\tincluding %.2f for functions and derivatives", T[4]);
	fprintf(f, "\n write: %9.2f\n total: %9.2f\n\n", T[3] - T[2],
		T[3] - T[0]);
	fflush(f);
	}

 void
MAIN__(VOID)
{
	char *stub;
	FILE *nl;
	fint INFORM, ITER, LFUVAL, M, MXROW, MXCOL, N, N0, NCALCF, NCALCG,
		NEL1, NELd, NELVAR, NG, NG0, NG1, NG2, NO, NVAR, NZ, ne;
	fint *GXEQX, *ICALCF, *ICALCG, *ICNA, *IELING, *ISTADA,
		*ISTADG, *ISTADH, *IVAR, *IWK, *KNDOFC;
	fint *Iev, *vi;
	fint dlen, iielts, len, liw0, lw0, iprint, ptest, ti;
	int firsttime, i, ig, ig1, j, jg, k, n, nlj, nslacks, ntot, nvtot;
	int *cm;
	real FOBJ, Times[5], mfac, t;
	real *A, *B, *C, *DGRAD, *ESCALE, *FT, *FUVALS, *FUVALS0, *GGRD, *GHES,
		*GSCALE, *GVALS, *Q, *U, *U1, *VSCALE, *WK, *XT, *lu, *x, *y;
	char *GNAMES, *VNAMES;
	ps_func *p, *pe, *pp;
	psb_elem *b, *be, **em;
	psg_elem *g, *ge, **gm;
	linpart *L, *Le;
	range *r;
	cgrad *cg, **cgp;
	ograd *og;
	cde *c;
	Hdispatch *hd, *hde;
	Jmp_buf err_jmp0;
	char buf[400];
	static fint L0, L6 = 6;
	typedef struct { char *msg; int code, wantsol; } Sol_info;
	Sol_info *SI;
	static Sol_info solinfo[9] = {
		{ "problem solved", 0, 1 },
		{ "too many iterations", 400, 1 },
		{ "trust region got too small", 510, 1 },
		{ "step got too small", 511, 1 },
		{ "an integer array was too short", 500, 0 },
		{ "a floating-point array was too short", 501, 0 },
		{ "a logical array was too short", 502, 0 },
		{ "bad KNDOFC", 503, 0 },
		{ "could not find a feasible solution", 520, 1 }
		};
	extern char **xargv;

	Times[0] = xectim_();
	asl = (ASL_pfgh*)ASL_alloc(ASL_read_pfgh);

	/* set defaults */
	dcomsb_.nsemib = 5;
	dcomal_.firstc = dcomal_.firstg = dcomal_.rmu = dcomal_.rmutol = .1;

	if (!(stub = getstops(xargv, &Oinfo)))
		return;
	--nprob;

#ifdef WANT_IEEE0
	if (call_ieee0) {
		extern void ieee0(VOID);
		ieee0();	/* catch creation of NaNs, Infinities */
		}
#endif

	if (!(nl = jacdim(stub, &M, &N, &NO, &NZ, &MXROW, &MXCOL,
			(fint)strlen(stub))))
		return;

	X0 = (real *)Malloc(3*N*sizeof(real));
	LUv = X0 + N;
	Uvx = LUv + N;
	if (M) {
		pi0 = (real *)Malloc(3*M*sizeof(real));
		LUrhs = pi0 + M;
		}
	i = ASL_find_default_no_groups | ASL_GJ_zerodrop;
	pfgh_read(nl, wantgr ? i | ASL_findOgroups : i);
	asl->P.pshv_g1 = 0;

	NELVAR = NG = nRmap = nvargp = 0;
	pp = 0;
	mfac = -1.;
	Iscr = (fint*)Malloc(N*(2*sizeof(fint)));
	Iscr1 = Iscr + N;
	memset(Iscr, 0, N*sizeof(fint));
	if (nprob >= 0 && nprob < NO) {
		if (!objscale)
			objscale = objtype[nprob] ? -1. : 1.;
		if (objscale < 0.)
			mfac = 1.;
		pp = asl->P.ops + nprob;
		NG = pp->nb + pp->ng;
		if (og = Ograd[nprob])
			for(NG++; og; og = og->next)
				LICNA++;
		nvargp = LICNA;
		NELVAR = nvsum(pp, 0);
		if (!NG)
			NG = 1;
		}
	NG0 = NG;
	NG += M;
	cgp = Cgrad;
	lu = LUrhs;
	nlj = nslacks = 0;
	for(p = asl->P.cps, pe = p + M; p < pe; p++) {
		if (lu[0] < lu[1])
			nslacks++;
		lu += 2;
		NELVAR += nvsum(p, cgp++);
		}
	free(Iscr);
	want_deriv = 1;
	N0 = N;
	if (nslacks) {
		LICNA += nslacks;
		nvargp += nslacks;
		n = n_var;
		N += nslacks;
		X0 = (real *)Realloc(X0, 3*N*sizeof(real));
		Uvx = shiftup(n, X0 + 2*n, X0 + 2*N);
		LUv = shiftup(n, X0 + n, X0 + N);
		C = (real *)Malloc(M*sizeof(real));
		ne = 0;
		conval(X0, C, &ne);
		if (ne)
			memset(C, 0, M*sizeof(real));
		}
	LFUVAL = NEL + NINTV + NINTV2 + 2*N + nvargp;
	NG1 = NG + 1;
	NELd = NEL1 = NEL + 1;
	if (NEL)
		NELd = NEL;	/* Fortran 77 requires array dimensions >= 1 */
	if (!NELVAR)
		NELVAR = 1;
	len = 3*NELd + 2*NG1 + NELVAR + 3*NEL1 + LICNA + 4*NG + N;
	NG2 = 2*NG;
	IELING = (fint*)Malloc(len*sizeof(fint));
	ISTADG = IELING + NELd;
	IELVAR = ISTADG + NG1;
	ISTAEV = IELVAR + NELVAR;
	INTVAR = ISTAEV + NEL1;
	ISTADH = INTVAR + NEL1;
	ICNA   = ISTADH + NEL1;
	ISTADA = ICNA   + LICNA;
	INTREP = ISTADA + NG1;
	KNDOFC = INTREP + NELd;
	ICALCF = KNDOFC + NG;
	ICALCG = ICALCF + NELd;
	GXEQX  = ICALCG + NG;
	IVAR   = GXEQX  + NG2;

	len = LICNA + 7*NG + 4*N + NELd + M;
	A = (real *)Malloc(len*sizeof(real));
	B      = A      + LICNA;
	GSCALE = B      + NG;
	ESCALE = GSCALE + NG;
	VSCALE = ESCALE + NELd;
	GVALS  = VSCALE + N;
	GGRD   = GVALS  + NG;
	GHES   = GGRD   + NG;
	FT     = GHES   + NG;
	XT     = FT     + NG;
	U      = XT     + N;
	Q      = U      + NG;
	DGRAD  = Q      + N;
	Yscrch = DGRAD  + N;

	if (M) {
		memset(Yscrch, 0, M*sizeof(real));
		U1 = U + NG0;
		for(i = 0; i < M; i++)
			U1[i] = mfac * pi0[i];
		}

	VNAMES = Malloc((N+NG)*10 + 1);
	namegen(VNAMES, N, 'x');
	GNAMES = VNAMES + 10*N;
	namegen(GNAMES, NG0, 'g');
	namegen(GNAMES + 10*NG0, M, 'c');

	Emap = (psb_elem**)Malloc((2*NEL+NG)*sizeof(psb_elem*)
					+ NEL*sizeof(int)
					+ N0*sizeof(fint)
					+ (2*N+nRmap)*sizeof(fint));
	Gmap = (psg_elem**)(Emap + NEL);
	Rmap = (fint**)(Gmap + NG);
	Rmapnext = (fint*)(Rmap + NEL);
	Iscr = Rmapnext + nRmap;
	Iscr1 = Iscr + N;
	Vused = Iscr1 + N;
	Cmap = (int*)(Vused + N0);
	for(i = 0; i < N; i++) {
		VSCALE[i] = 1.;
		Iscr[i] = 0;
		}
	memset(Vused, 0, N0*sizeof(fint));

	len = N;
	if (ICHOSE[1] >= 11)
		len += N;
	if (len < NEL)
		len = NEL;
	if ((liw0 = N) < NG)
		liw0 = NG;
	liw0 += 4*N + len + maxsin*(maxsin + 1)
		+ 2*NG + NEL + 3*nvargp + 3;
	/* liw0 is still missing 2*iielts, which is computed below. */

	if ((len = NG) < maxsel)
		len = maxsel;
	len += 2*maxsin;
	if ((lw0 = maxsel) < NINTV)
		lw0 = NINTV;
	lw0 += N + NINTV;
	if (lw0 < len)
		lw0 = len;
	lw0 += 9*N;

	NFCALL = NGCALL = 0;
	firsttime = 1;
	--pfreq;
	FUVALS = (real*)Malloc(LFUVAL*sizeof(real));
	FUVALS0 = FUVALS - 1;
	hesset(1, (int)nprob, 1, 0, n_con);
	timing &= 3;
	if (cl_timing && !timing)
		timing = 1;
	if (timing) {
		Times[1] = xectim_();
		Times[4] = 0;
		}
 retry:
	ITER = iprint = 0;
	ptest = OUTLEV > 0 && pstop > pstart && pfreq >= 0 ? pstart : -1;
	n = n_var;
	cm = Cmap;
	em = Emap--;
	gm = Gmap;
	Rmapnext = (fint*)(Rmap + NEL);

	Iev = IELVAR;
	ISTAEV[0] = ISTADG[0] = ISTADA[0] = 1;
	i = ig1 = 1;
	j = ig = jg = 0;
	nmax = 0;
	if (pp) {
		t = constterm(pp, obj_de + nprob);
		if (og = Ograd[nprob]) {
			*gm++ = 0;
			for(; og; og = og->next) {
				A[jg] = og->coef;
				ICNA[jg++] = og->varno + 1;
				}
			B[0] = t;
			t = 0;
			GSCALE[0] = objscale;
			GXEQX[0] = 1;
			ISTADA[ig1] = jg + 1;
			ig = ig1++;
			ISTADG[i++] = em - Emap;
			}
		for(b = pp->b, be = b + pp->nb; b < be; b++) {
			*gm++ = 0;
			B[ig] = t;
			t = 0.;
			GSCALE[ig] = objscale;
			GXEQX[ig] = 1;
			ISTADA[ig1] = jg + 1;
			ig = ig1++;
			if ((r = b->U)
			 || pp->nb != 1
			 || pp->b->D.e->op != f_OPNUM) {
				*cm++ = -1;
				*em++ = b;
				Iev = elvars(Iev, r, j++);
				ISTAEV[j] = Iev - IELVAR + 1;
				}
			ISTADG[i++] = em - Emap;
			}
		for(g = pp->g, ge = g + pp->ng; g < ge; g++) {
			*gm++ = g;
			if (k = g->nlin) {
				L = g->L;
				Le = L + k;
				do {
					A[jg] = L->fac;
					ICNA[jg++] = L->v.i + 1;
					}
					while(++L < Le);
				}
			B[ig] = -g->g0;
			GSCALE[ig] = g->scale * objscale;
			GXEQX[ig] = 0;
			ISTADA[ig1] = jg + 1;
			ig = ig1++;
			for(b = g->E, be = b + g->ns; b < be; b++) {
				*cm++ = -1;
				*em++ = b;
				Iev = elvars(Iev, b->U, j++);
				ISTAEV[j] = Iev - IELVAR + 1;
				}
			ISTADG[i++] = em - Emap;
			}
		if (gm == Gmap) { /* minimize 1 */
			*gm++ = 0;
			B[ig] = t;
			GSCALE[ig] = objscale;
			GXEQX[ig] = 1;
			ISTADA[ig1] = jg + 1;
			ig = ig1++;
			ISTADG[i++] = em - Emap;
			}
		}

	cgp = Cgrad;
	lu = LUrhs;
	c = con_de;
	for(p = asl->P.cps, pe = p + M; p < pe; p++) {
		*gm++ = 0;
		cg = *cgp++;
		for(b = p->b, be = b + p->nb; b < be; b++)
			if ((r = b->U)
			 || p->nb != 1
			 || p->b->D.e->op != f_OPNUM) {
				*cm++ = p - asl->P.cps;
				*em++ = b;
				Iev = elvars(Iev, r, j++);
				ISTAEV[j] = Iev - IELVAR + 1;
				}
		for(; cg; cg = cg->next) {
			A[jg] = cg->coef;
			ICNA[jg++] = cg->varno + 1;
			}
		t = constterm(p, c++);
		if (lu[0] < lu[1]) {
			LUv[n] = lu[0] + t;
			Uvx[n] = lu[1] + t;
			A[jg] = -1.;
			X0[n] = C[p - asl->P.cps] + t;
			t = 0.;
			ICNA[jg++] = ++n;
			}
		else
			t += lu[0];
		B[ig] = t;
		lu += 2;
		GSCALE[ig] = 1.;
		GXEQX[ig] = 1;
		ISTADA[ig1] = jg + 1;
		ig = ig1++;
		ISTADG[i++] = em - Emap;
		}

	for(i = 0; i < NEL; i++) {
		ESCALE[i] = 1.;
		IELING[i] = i+1;
		}

	for(i = 0; i < NG0; i++)
		KNDOFC[i] = 1;
	while(i < NG)
		KNDOFC[i++] = 2;
	++Emap;

	if (firsttime) {
		firsttime = 0;
		dlen = len = 10*(N+M+10) + N*N + lwinc;
		if (LIWK)
			dlen = 0;
		else {
			iielts = N;
			vi = Vused + N0;
			while(vi > Vused)
				if ((ti = *--vi - 1) > 0)
					iielts += ti;
			liw0 += 2*iielts;
			LIWK = liw0 + len;
			}
		if (LWK)
			dlen = 0;
		else
			LWK = lw0 + len;

		liw0 = LIWK;
		lw0 = LWK;
		}

	IWK = (fint*)Malloc(LIWK*sizeof(fint));
	WK = (real*)Malloc(LWK*sizeof(real));
	if (nmax) {
		Hscratch = (real *)Malloc(N*(sizeof(real)+sizeof(Hdispatch*))
					+ nmax*sizeof(Hdispatch));
		hdisp = (Hdispatch**)(Hscratch + N);
		memset(Hscratch, 0, N*(sizeof(Hdispatch*)+sizeof(real)));
		hdfree = hd = (Hdispatch *)(hdisp+N);
		hde = hd + nmax - 1;
		while(hd < hde)
			hd = hd->next = hd + 1;
		hd->next = 0;
		}
	INFORM = 0;
	err_jmp = &err_jmp0;
	if (setjmp(err_jmp0.jb)) {
		INFORM = -11;	/* 1/0, sqrt'(0), etc. */
		nlj++;
		}
	goto for_start;
	for(;;) {
		if (timing)
			Times[4] += xectim_() - t;
 for_start:
		if (ITER == ptest) {
			if (ITER == pstop) {
				iprint = 0;
				ptest = -1;
				}
			else if (iprint) {
				iprint = 0;
				ptest += pfreq;
				if (ptest > pstop)
					ptest = -1;
				}
			else {
				iprint = OUTLEV;
				if (pfreq) {
					ptest += pfreq;
					if (ptest > pstop)
						ptest = pstop;
					}
				else
					ptest = pstop;
				}
			}
		dauglg_(&N, &NG, &NEL, IELING, &NELd, ISTADG,
			&NG1, IELVAR, &NELVAR, ISTAEV, &NEL1, INTVAR, &NEL1,
			ISTADH, &NEL1, ICNA, &LICNA, ISTADA, &NG1,
			A, &LICNA, B, &NG, LUv, &N, Uvx, &N,
			GSCALE, &NG, ESCALE, &NELd, VSCALE, &N,
			GXEQX, &NG2, INTREP, &NELd, KNDOFC, &NG,
			ranges, &INFORM, &FOBJ, X0, &N, U, &NG, GVALS,
			&NG, FT, &NG, FUVALS, &LFUVAL, XT,
			&N, ICALCF, &NELd, &NCALCF, ICALCG, &NG,
			&NCALCG, IVAR, &N, &NVAR, Q, &N, DGRAD,
			&N, ICHOSE, &ITER, &MAXIT, &L0, VNAMES,
			&N, GNAMES, &NG, &STOPG, &STOPC, IWK, &LIWK, WK, &LWK,
			&iprint, &L6, 10L, 10L);

		if (INFORM >= 0)
			break;
		if (timing)
			t = xectim_();
		switch(INFORM+7) {

		  case 0: /* INFORM == -7 */
		  case 4: /* INFORM == -3 */
			funcvals(ICALCF, (int)NCALCF, FUVALS, XT);
			continue;

		  case 2: /* INFORM == -5 */
			groupgrads(ICALCG, (int)NCALCG, GGRD, GHES);
			/* no break */

		  case 1: /* INFORM == -6 */
			funcgrads(ICALCF, (int)NCALCF, FUVALS0, INTVAR, ISTADH);
			continue;

		  case 3: /* INFORM == -4 */
			groupvals(ICALCG, (int)NCALCG, FT, GVALS);
			continue;

		  case 5: /* INFORM == -2 */
			groupvals(ICALCG, (int)NCALCG, FT, GVALS);
			groupgrads(ICALCG, (int)NCALCG, GGRD, GHES);
			continue;

		  case 6: /* INFORM == -1 */
			funcvals(ICALCF, (int)NCALCF, FUVALS, XT);
			funcgrads(ICALCF, (int)NCALCF, FUVALS0, INTVAR, ISTADH);
			continue;

		  default:
			bad_inform(INFORM);
		  }
		}
	if (INFORM > 8)
		bad_inform(INFORM);
	switch(INFORM) {
	 case 4:
	 case 5:
		if (dlen) {
			if (dlen < lincmax) {
				LIWK += dlen;
				LWK += dlen;
				dlen += dlen;
 retry1:
				if (nmax)
					free(Hscratch);
				free(IWK);
				free(WK);
				goto retry;
				}
			}
		else
			if (LIWK < lincmax && LWK < lincmax) {
				if ((LIWK <<= 1) > lincmax)
					LIWK = lincmax;
				if ((LWK <<= 1) > lincmax)
					LWK = lincmax;
				goto retry1;
				}
		fprintf(Stderr, "LIWK = %ld, LWK = %ld, LFUVAL = %ld\n",
			LIWK, LWK, LFUVAL);
	 }
	if (timing)
		Times[2] = xectim_();
	SI = solinfo + INFORM;
	solve_result_num = SI->code;
	i = sprintf(buf, "LANCELOT: %s.", SI->msg);
	x = y = 0;
	if (SI->wantsol) {
		x = X0;
		y = pi0;
		i += sprintf(buf+i, "\n%ld iterations; ", ITER);
		if (nprob >= 0 && nprob < NO) {
			i += sprintf(buf+i, "objective ");
			i += g_fmtop(buf+i, FOBJ/objscale);
			}
		else
			i += sprintf(buf+i, "no objective.");
		for(j = 0; j < M; j++)
			pi0[j] = mfac * U1[j];
		}
	i += sprintf(buf+i,
		"\nM = %ld, N = %ld, nslacks = %d, NEL = %ld, NG = %ld",
		M, N, nslacks, NEL, NG);
	for(j = k = ntot = nvtot = 0; j < NEL; j++)
		if (Rmap[j]) {
			k++;
			r = Emap[j]->U;
			ntot += r->n;
			nvtot += r->nv;
			}
	if (k) {
		j = k == 1;
		i += sprintf(buf+i,
		  "\n%d nontrivial linear map%s reduce%s %d variables to %d.",
			k, "s" + j, "s" + (1-j), nvtot, ntot);
		}
	i += sprintf(buf+i, "\nNFCALL = %ld, NGCALL = %ld", NFCALL, NGCALL);
	if (LIWK > liw0)
		i += sprintf(buf+i, "\nLIWK = %ld (started at %ld)",
			LIWK, liw0);
	if (LWK > lw0)
		i += sprintf(buf+i, "\nLWK = %ld (started at %ld)",
			LWK, lw0);
	if (nlj)
		sprintf(buf+i, "\n%d failed function evaluations.", nlj);
	write_sol(buf, x, y, &Oinfo);
	if (timing & 3) {
		Times[3] = xectim_();
		if (timing & 1)
			time_out(stdout, Times);
		if (timing & 2)
			time_out(Stderr, Times);
		}
	}

#ifdef __cplusplus
	}
#endif
