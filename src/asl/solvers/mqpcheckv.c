/****************************************************************
Copyright (C) 1997-1999, 2000 Lucent Technologies
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

#include "avltree.h"
#include "nlp.h"
#include "obj_adj.h"
#include "r_qp.hd"

 typedef struct
dyad {
	struct dyad *next;
	ograd *Lq, *Rq;
	} dyad;

 typedef struct
term {
	dyad	*Q, *Qe;
	ograd	*L, *Le;
	} term;

 typedef struct
dispatch {
	struct dispatch *next;
	int i, j, jend;
	} dispatch;

 enum {Memblock_gulp = 8190};

 typedef struct
Memblock {
	struct Memblock *next, *prev;
	void *x[Memblock_gulp];
	} Memblock;

 typedef struct
Static {
	ASL_fg *asl;
	AVL_Tree *AQ;
	Memblock *mb, *mb0, *mblast;
	dispatch *cd0, **cdisp;
	double *_s_x;
	dyad *_freedyad, **_s_q;
	int *_s_s, *_s_z, *w, *zct;
	ograd *_freeog, **oq;
	term **_cterms, *_freeterm;
	void **v, **ve;
	int _zerodiv, nvinc, nzct;
	} Static;

#define cterms		S->_cterms
#define freedyad	S->_freedyad
#define freeog		S->_freeog
#define freeterm	S->_freeterm
#define s_q		S->_s_q
#define s_s		S->_s_s
#define s_x		S->_s_x
#define s_z		S->_s_z
#define zerodiv		S->_zerodiv

 static void*
M2alloc(Static *S, size_t len, int dbl_align)
{
	Memblock *mb, *mb1;
	size_t n;
	void *rv;

	if (dbl_align)
		S->v = (void**)(((size_t)S->v + (sizeof(void*)-1)) & ~(sizeof(void*)-1));
	n = (len + (sizeof(void*)-1))/sizeof(void*);
	if (S->v + n >= S->ve) {
		mb = S->mb;
		if (!(mb1 = mb->next)) {
			mb1 = Malloc(sizeof(Memblock));
			S->mblast = mb->next = mb1;
			mb1->prev = mb;
			mb1->next = 0;
			}
		S->mb = mb1;
		S->v = mb1->x;
		S->ve = S->v + Memblock_gulp;
		}
	rv = (void*)S->v;
	S->v += n;
	return rv;
	}

 static void
free_term(Static *S, term *t)
{
	t->Q = (dyad *)freeterm;
	freeterm = t;
	}

 static term *
new_term(Static *S, ograd *o)
{
	term *rv;

	if ((rv = freeterm))
		freeterm = (term *)rv->Q;
	else
		rv = (term*)M2alloc(S, sizeof(term), 0);
	rv->L = rv->Le = o;
	rv->Q = rv->Qe = 0;
	return rv;
	}

 static void
free_og(Static *S, ograd *o)
{
	o->next = freeog;
	freeog = o;
	}

 static ograd *
new_og(Static *S, ograd *next, int i, real v)
{
	ograd *rv;

	if ((rv = freeog))
		freeog = rv->next;
	else
		rv = (ograd*)M2alloc(S, sizeof(ograd), 1);
	rv->next = next;
	rv->varno = i;
	rv->coef = v;
	return rv;
	}

 static ograd *
ogdup(Static *S, ograd *og, ograd **oge)
{
	ograd *og0, *og1;

	og0 = og1 = new_og(S, 0, og->varno, og->coef);
	while((og = og->next))
		og1 = og1->next = new_og(S, 0, og->varno, og->coef);
	if (oge)
		*oge = og1;
	return og0;
	}

 static int
count(Static *S, ograd **ogp)
{
	int i, rv, nz, *s, *z;
	double t, *x;
	ograd *og, *og1;

	s = s_s;
	x = s_x;
	z = s_z;

	t = 0;
	nz = rv = 0;
	for(og = *ogp; og; og = og1) {
		og1 = og->next;
		if ((i = og->varno) < 0)
			t += og->coef;
		else if (!s[i]++)
			x[z[nz++] = i] = og->coef;
		else
			x[i] += og->coef;
		free_og(S, og);
		}
	while(nz > 0) {
		s[i = z[--nz]] = 0;
		if (x[i]) {
			og = new_og(S, og, i, x[i]);
			rv++;
			}
		}
	if (t)
		og = new_og(S, og, -1, t);
	*ogp = og;
	return rv;
	}

 static void
free_dyad(Static *S, dyad *t)
{
	t->next = freedyad;
	freedyad = t;
	}

 static dyad *
new_dyad(Static *S, dyad *next, ograd *L, ograd *R, int permute)
{
	dyad *rv;
	ograd *t;

	if (permute) {
		if (L == R) {
			count(S, &L);
			R = L;
			}
		else if (count(S, &L) > count(S, &R)) {
			t = L;
			L = R;
			R = t;
			}
		if (!L) /* e.g., 0*x*x from <<0;0,0>>x*x */
			/* with AMPL version < 20000216. */
			return next;
		}
	if ((rv = freedyad))
		freedyad = rv->next;
	else
		rv = (dyad*)M2alloc(S, sizeof(dyad), 0);
	rv->next = next;
	rv->Lq = L;
	rv->Rq = R;
	return rv;
	}

 static term *
termsum(Static *S, term *L, term *R)
{
	if (!L || !R)
		return 0;
	if (L->Qe && (L->Qe->next = R->Q))
		L->Qe = R->Qe;
	else if (R->Q) {
		L->Q = R->Q;
		L->Qe = R->Qe;
		}
	if (L->Le && (L->Le->next = R->L))
		L->Le = R->Le;
	else if (R->L) {
		L->L = R->L;
		L->Le = R->Le;
		}
	free_term(S, R);
	return L;
	}

 static term *
scale(Static *S, term *T, register double t)
{
	register ograd *og;
	register dyad *d;

	if (T) {
		for(d = T->Q; d; d = d->next) {
			if (d->Lq == d->Rq)
				d->Rq = ogdup(S, d->Lq, 0);
			for(og = d->Lq; og; og = og->next)
				og->coef *= t;
			}
		for(og = T->L; og; og = og->next)
			og->coef *= t;
		}
	return T;
	}

 static term *ewalk ANSI((Static*, expr*));

 static term *
comterm(Static *S, int i)
{
	int nlin;
	cexp *c;
	cexp1 *c1;
	linpart *L, *Le;
	expr_v ev, *vp;
	term *T;
	ASL_fg* asl = S->asl;

	if (i < ncom0) {
		c = cexps + i;
		T = ewalk(S, c->e);
		L = c->L;
		nlin = c->nlin;
		}
	else {
		c1 = cexps1 + (i - ncom0);
		T = ewalk(S, c1->e);
		L = c1->L;
		nlin = c1->nlin;
		}
	if (L && T)
		for(Le = L + nlin; L < Le; L++) {
			vp = (expr_v *)((char *)L->v.rp -
				((char *)&ev.v - (char *)&ev));
			T = termsum(S, T, new_term(S,
				new_og(S, 0, (int)(vp - var_e), L->fac)));
			}
	return T;
	}

 static term *
termdup(Static *S, term *T)
{
	term *rv;
	ograd *og, *oge;
	dyad *Q, *Q1;

	if ((og = oge = T->L))
		og = ogdup(S, og, &oge);
	rv = new_term(S, og);
	rv->Le = oge;
	if (!(Q = T->Q))
		return rv;
	Q1 = rv->Qe = new_dyad(S, 0, ogdup(S, Q->Lq,0), ogdup(S, Q->Rq,0), 1);
	while((Q = Q->next))
		Q1 = new_dyad(S, Q1, ogdup(S, Q->Lq,0), ogdup(S, Q->Rq,0), 1);
	rv->Q = Q1;
	return rv;
	}

 static term *
ewalk(Static *S, expr *e)
{
	term *L, *R, *T;
	ograd *o, *oR;
	expr **ep, **epe;
	int i;
	ASL_fg *asl;

	switch(Intcast e->op) {

	  case OPNUM:
		return new_term(S, new_og(S, 0, -1 , ((expr_n *)e)->v));

	  case OPPLUS:
		return termsum(S, ewalk(S, e->L.e), ewalk(S, e->R.e));

	  case OPMINUS:
		return termsum(S, ewalk(S, e->L.e), scale(S, ewalk(S, e->R.e), -1.));

	  case OPUMINUS:
		return scale(S, ewalk(S, e->L.e), -1.);

	  case OPMULT:
		if (!(L = ewalk(S, e->L.e))
		 || !(R = ewalk(S, e->R.e)))
			break;
		if (L->Q) {
			if (R->Q)
				break;
 qscale:
			o = R->L;
			if (o->next || o->varno >= 0)
				break;
			scale(S, L, o->coef);
			free_og(S, o);
			free_term(S, R);
			return L;
			}
		if (R->Q) {
			T = L;
			L = R;
			R = T;
			goto qscale;
			}
		o = L->L;
		oR = R->L;
		if (o->next || o->varno >= 0) {
			if (oR->next || oR->varno >= 0) {
				L->Q = L->Qe = new_dyad(S, 0,o,oR,1);
				L->L = L->Le = 0;
				}
			else {
				scale(S, L, oR->coef);
				free_og(S, oR);
				}
			free_term(S, R);
			return L;
			}
		scale(S, R, o->coef);
		free_og(S, o);
		free_term(S, L);
		return R;

	  case OPDIV:
		/* only allow division by a constant */
		if (!(R = ewalk(S, e->R.e)))
			break;
		o = R->L;
		if (R->Q || o->next || o->varno >= 0)
			break;
		if (!(L = ewalk(S, e->L.e)))
			break;
		if (!o->coef) {
			zerodiv++;
			L = 0;
			}
		else
			scale(S, L, 1./o->coef);
		free_og(S, o);
		free_term(S, R);
		return L;

	  case OPSUMLIST:
		ep = e->L.ep;
		epe = e->R.ep;
		L = ewalk(S, *ep);
		while(L && ++ep < epe)
			L = termsum(S, L, ewalk(S, *ep));
		return L;

	  case OP2POW:
		L = ewalk(S, e->L.e);
		if (!L || L->Q)
			break;
		o = L->L;
		if (!o->next && o->varno < 0) {
			o->coef *= o->coef;
			return L;
			}
		L->Q = L->Qe = new_dyad(S, 0,o,o,1);
		L->L = L->Le = 0;
		return L;

	  case OPVARVAL:
		asl = S->asl;
		if ((i = (expr_v *)e - var_e) < n_var)
			return new_term(S, new_og(S, 0, i, 1.));
		i -= S->nvinc;
		if (!(L = cterms[i -= n_var])) {
			if (!(L = comterm(S, i)))
				return 0;
			cterms[i] = L;
			S->zct[S->nzct++] = i;
			}
		return termdup(S, L);
		}
	return 0; /* nonlinear */
	}

#ifdef __cplusplus
extern "C" {
 static int comp(const void*, const void*, void*),
	    lcmp(const void*, const void*, void*);
}
#endif

 static int
comp(const void *a, const void *b, void *v)
{
	Not_Used(v);
	return (*(ograd **)a)->varno - (*(ograd **)b)->varno;
	}

 static ograd *
sortq(ograd *og0, ograd **q)
{
	ograd *og, **q1;
	int n;

	for(q1 = q, og = og0; og; og = og->next)
		*q1++ = og;
	if ((n = q1 - q) > 1) {
		qsortv(q, n, sizeof(ograd *), comp, NULL);
		og0 = 0;
		do {
			og = *--q1;
			og->next = og0;
			og0 = og;
			} while(q1 > q);
		}
	return og0;
	}

 static double
dsort(Static *S, term *T, ograd **q, cgrad **cgp, ograd **ogp, int arrays)
{
	cgrad *cg;
	ograd *og, *og1;
	double t, t1, rv, *x;
	dyad *Q;

	x = s_x;

	rv = 0;
	count(S, &T->L);	/* accumulate */
	if (arrays) {
		if (ogp)
			for(og = *ogp; og; og = og->next)
				x[og->varno] = og->coef;
		else
			for(cg = *cgp; cg; cg = cg->next)
				x[cg->varno] = cg->coef;
		}
	if ((og = T->L) && og->varno < 0) {
		rv = og->coef;
		og = og->next;
		}
	for(; og; og = og->next)
		x[og->varno] += og->coef;

	for(Q = T->Q; Q; Q = Q->next) {
		og = Q->Lq;
		og1 = Q->Rq;
		t = t1 = 0;
		if (og->varno < 0) {
			t = og->coef;
			og = og->next;
			}
		if (og1->varno < 0) {
			t1 = og1->coef;
			Q->Rq = og1 = og1->next;
			rv += t*t1;
			}
		if (t)
			for(; og1; og1 = og1->next)
				x[og1->varno] += t*og1->coef;
		if (t1)
			for(og1 = og; og1; og1 = og1->next)
				x[og1->varno] += t1*og1->coef;
		Q->Lq = sortq(og, q);
		Q->Rq = og == Q->Rq ? Q->Lq : sortq(Q->Rq, q);
		}
	if (arrays) {
		if (ogp)
			for(og = *ogp; og; og = og->next)
				og->coef = x[og->varno];
		else
			for(cg = *cgp; cg; cg = cg->next)
				cg->coef = x[cg->varno];
		}
	return rv;
	}

 static void
free_oglist(Static *S, ograd *og)
{
	ograd *og1;

	for(; og; og = og1) {
		og1 = og->next;
		free_og(S, og);
		}
	}

 static void
cterm_free(Static *S)
{
	dyad *d, *d1;
	term **ct, *t;
	int i, n, *zct;

	if (!(n = S->nzct))
		return;
	zct = S->zct;
	S->nzct = 0;
	ct = cterms;
	while(n > 0) {
		if ((t = ct[i = zct[--n]])) {
			ct[i] = 0;
			free_oglist(S, t->L);
			d1 = t->Q;
			while((d = d1)) {
				d1 = d->next;
				free_oglist(S, d->Lq);
				if (d->Rq != d->Lq)
					free_oglist(S, d->Rq);
				free_dyad(S, d);
				}
			}
		}
	}

 static int
lcmp(const void *a, const void *b, void *v)
{
	Not_Used(v);
	return (int)(*(int *)a - *(int *)b);
	}

 static int
vcomp(void *v, const Element *a, const Element *b)
{
	if (a == b)
		return 0;
	return a < b ? -1 : 1;
	}

 ssize_t
mqpcheckv_ASL(ASL *a, int co, QPinfo **QPIp, void **vp)
{
	ASL_fg *asl;
	AVL_Node *NQ, *NQ0;
	AVL_Tree *AQ;
	Memblock *mb;
	QPinfo *qpi;
	Objrep *od, **pod;
	Static *S;
	cde *c;
	cgrad *cg, **cgp, **cgq, *cq;
	dispatch *cd, *cd0, **cdisp, **cdisp0, *cdnext, **cdp;
	dyad *d, *d1, **q, **q1, **q2;
	expr *e;
	expr_n *en;
	int *cm, *colno, *qm, *rowq, *rowq0, *rowq1, *s, *vmi, *w, *z;
	int arrays, co0, ftn, i, icol, j, ncol, ncom, nv, nva, nz, pass;
	ograd *og, *og1, *og2, **ogp;
	real *L, *U, *delsq, *delsq0, *delsq1, objadj, t, *x;
	size_t  *colq, *colq1, nelq;
	term *T;

	ASL_CHECK(a, ASL_read_fg, "nqpcheck");
	asl = (ASL_fg*)a;
	if (co >= n_obj || co < -n_con)
		return -3L;
	colno = 0;
	od = 0;
	co0 = co;
	if (co >= 0) {
		if ((pod = asl->i.Or) && (od = pod[co])) {
			co = od->ico;
			goto use_Cgrad;
			}
		else {
			c = obj_de + co;
			ogp = Ograd + co;
			cgp = 0;
			}
		}
	else {
		co = -1 - co;
		if ((cm = asl->i.cmap))
			co = cm[co];
 use_Cgrad:
		c = con_de + co;
		if (!(cgp = asl->i.Cgrad0))
			cgp = Cgrad;
		cgp += co;
		ogp = 0;
		}

	e = c->e;
	if (e->op == f_OPNUM)
		return 0;

	if (asl->i.vmap && !asl->i.vminv)
		get_vminv_ASL(a);
	nv = n_var;
	ncom = ncom0 + ncom1;
	if (!(S = *(Static**)vp)) {
		i = asl->i.n_var0 + asl->i.nsufext[0];
		if ((nva = nv) < i)
			nva = i;
		x = (double *)Malloc(nva*(sizeof(double)
					+sizeof(dyad*)
					+sizeof(ograd*)
					+sizeof(dispatch*)
					+sizeof(dispatch)
					+3*sizeof(int))
					+ sizeof(Memblock)
					+ sizeof(Static));
		mb = (Memblock*)(x + nva);
		mb->prev = mb->next = 0;
		S = (Static*)(mb + 1);
		*vp = (void*)S;
		memset(S, 0, sizeof(Static));
		S->mb0 = S->mblast = mb;
		s_x = x;
		S->asl = asl;
		s_q = q = (dyad**)(S+1);
		S->oq = (ograd**)(q + nva);
		S->cdisp = cdisp = (dispatch**)(S->oq + nva);
		S->cd0 = cd0 = (dispatch*)(cdisp + nva);
		s_z = z = (int*)(cd0 + nva);
		s_s = s = z + nva;
		S->w = (int*)(s + nva);
		memset(s, 0, nva*sizeof(int));
		memset(cdisp, 0, nva*sizeof(dispatch*));
		memset(q, 0, nva*sizeof(dyad *));
		memset(S->w, 0, nva*sizeof(int));
		if (ncom) {
			cterms = (term **)Malloc(ncom*(sizeof(term*)+sizeof(int)));
			memset(cterms, 0, ncom*sizeof(term*));
			S->zct = (int*)(cterms + ncom);
			}
		S->AQ = AVL_Tree_alloc2(0, vcomp, mymalloc, 0);
		}
	else {
		q = s_q;
		x = s_x;
		z = s_z;
		s = s_s;
		cdisp = S->cdisp;
		cd0 = S->cd0;
		}
	S->mb = mb = S->mb0;
	S->v  = &mb->x[0];
	S->ve = &mb->x[Memblock_gulp];
	w = S->w;
	freedyad = 0;
	freeog = 0;
	freeterm = 0;
	AQ = S->AQ;
	ftn = Fortran;
	cdisp0 = cdisp - ftn;
	S->nvinc = nv - asl->i.n_var0 + asl->i.nsufext[ASL_Sufkind_var];

	delsq = delsq0 = delsq1 = 0; /* silence buggy "not-initialized" warning */
	colq = colq1 = 0;				/* ditto */
	rowq = rowq0 = rowq1 = 0;			/* ditto */

	arrays = 0;
	if (QPIp) {
		*QPIp = 0;
		arrays = 1;
		}
	zerodiv = 0;
	if (!(T = ewalk(S, e)) || zerodiv)
		return T ? -2L : -1L;

	if (S->nzct)
		cterm_free(S);
	if (od) {
		cgq = &od->cg;
		for(i = 0, cg = *cgp; cg; cg = cg->next) {
			if (cg->coef != 0.)
				++i;
			}
		if (i) {
			cq = M1alloc(i*sizeof(cgrad));
			for(cg = *cgp; cg; cg = cg->next) {
				*cgq = cq;
				cgq = &cq->next;
				*cq = *cg;
				++cq;
				}
			}
		*cgq = 0;
		}

	q = s_q;
	objadj = dsort(S, T, S->oq, cgp, ogp, arrays);

	nelq = ncol = nz = 0;
	qpi = 0;
	for(pass = 0; pass < 2; pass++) {
		if (pass) {
			if (!nelq)
				break;
			if (!arrays) {
				for(qm = (int*)AVL_first(AQ, &NQ); qm; ) {
					*qm = 0;
					NQ0 = NQ;
					qm = (int*) AVL_next(&NQ);
					AVL_delnode(AQ, &NQ0);
					}
				break;
				}
			qpi = *QPIp = (QPinfo*)Malloc(sizeof(QPinfo)
						+ nelq*(sizeof(real) + sizeof(int))
						+ ncol*sizeof(int)
						+ (ncol + 1)*sizeof(size_t));
			qpi->nz = nelq;
			qpi->nc = ncol;
			qpi->delsq = delsq = delsq1 = (double *)(qpi+1);
			qpi->colbeg = colq = (size_t *)(delsq + nelq);
			qpi->rowno = rowq = (int *)(colq + ncol + 1);
			qpi->colno = colno = rowq + nelq;
			nelq = ftn;
			delsq0 = delsq - ftn;
			rowq0 = rowq - ftn;
			for(d = T->Q; d; d = d->next) {
				og = d->Rq;
				og1 = d->Lq;
				i = og->varno;
				while(og1 && og1->varno < i)
					og1 = og1->next;
				if (og1) {
					q1 = q + i;
					if (!w[i]++)
						AVL_vinsert(AQ, 0, (Element*)&w[i], 0);
					*q1 = new_dyad(S, *q1, og, og1, 0);
					}
				og1 = d->Lq;
				i = og1->varno;
				while(og && og->varno < i)
					og = og->next;
				if (og) {
					q1 = q + i;
					if (!w[i]++)
						AVL_vinsert(AQ, 0, (Element*)&w[i], 0);
					*q1 = new_dyad(S, *q1, og1, og, 0);
					}
				}
			}
		else {
			for(d = T->Q; d; d = d->next) {
				og = d->Rq;
				og1 = d->Lq;
				q1 = q + (i = og->varno);
				if (!w[i]++) {
					++ncol;
					AVL_vinsert(AQ, 0, (Element*)&w[i], 0);
					}
				*q1 = new_dyad(S, *q1, og, og1, 0);
				q1 = q + (i = og1->varno);
				if (!w[i]++) {
					++ncol;
					AVL_vinsert(AQ, 0, (Element*)&w[i], 0);
					}
				*q1 = new_dyad(S, *q1, og1, og, 0);
				}
			}
		vmi = asl->i.vmap ? get_vminv_ASL((ASL*)asl) : 0;
		for(qm = (int*)AVL_first(AQ, &NQ); qm; ) {
			NQ0 = NQ;
			icol = qm - w;
			if (pass) {
				*qm = 0;
				*colno++ = icol + ftn;
				*colq++ = nelq;
				if ((cd = cdisp[icol])) {
				    cdisp[icol] = 0;
				    do {
					cdnext = cd->next;
					s[i = cd->i]++;
					x[z[nz++] = i] = delsq0[cd->j++];
					if (cd->j < cd->jend) {
						cdp = cdisp0 + rowq0[cd->j];
						cd->next = *cdp;
						*cdp = cd;
						}
					} while((cd = cdnext));
				    }
				}
			if ((d = q[icol])) {
			    q[icol] = 0;
			    do {
				og = d->Lq;
				og1 = d->Rq;
				switch(pass) {
				  case 0:
					for(; og1; og1 = og1->next)
						if (!s[i = og1->varno]++)
							z[nz++] = i;
					break;
				  case 1:
					t = og->coef;
					for(; og1; og1 = og1->next) {
						if (!s[i = og1->varno]++)
							x[z[nz++] = i] =
								t*og1->coef;
						else
							x[i] += t*og1->coef;
						}
					if ((og1 = og->next)) {
					  og2 = d->Rq;
					  while (og2->varno < og1->varno)
					    if (!(og2 = og2->next)) {
						while((og1 = og->next))
							og = og1;
						break;
						}
					  d->Rq = og2;
					  }
					}
				d1 = d->next;
				if ((og = og->next)) {
					i = og->varno;
					if (pass) {
						og1 = d->Rq;
						while(og1->varno < i)
							if (!(og1 = og1->next))
								goto d_del;
						d->Rq = og1;
						}
					d->Lq = og;
					q2 = q + i;
					if (!w[i]++) {
						++ncol;
						AVL_vinsert(AQ, 0, (Element*)&w[i], 0);
						}
					d->next = *q2;
					*q2 = d;
					}
				else {
 d_del:
					free_dyad(S, d);
					}
				}
				while((d = d1));
			    }
			if (nz) {
				if (pass) {
					if (nz > 1)
						qsortv(z, nz, sizeof(int), lcmp, NULL);
					for(i = 0; i < nz; i++) {
						if ((t = x[j = z[i]])) {
							*delsq++ = t;
							if (vmi)
								j = vmi[j];
							*rowq++ = j + ftn;
							nelq++;
							}
						s[j] = 0;
						}
					for(i = 0; i < nz; i++)
					    if ((j = z[i]) > icol && x[j]) {
						cd0->i = icol;
						cd0->j = colq[-1] + i;
						cd0->jend = nelq;
						cdp = cdisp + j;
						cd0->next = *cdp;
						*cdp = cd0++;
						if (!w[j]++)
							AVL_vinsert(AQ, 0, (Element*)&w[j], 0);
						break;
						}
					nz = 0;
					}
				else {
					nelq += nz;
					while(nz > 0)
						s[z[--nz]] = 0;
					}
				}
			qm = (int*) AVL_next(&NQ);
			if (pass)
				AVL_delnode(AQ, &NQ0);
			}
		}
	if (colq)
		*colq = nelq;
	if (arrays) {
		if (nelq)
			nelq -= ftn;
		en = (expr_n *)mem(sizeof(expr_n));
		en->op = f_OPNUM_ASL;
		if (od) {
			od->opify = qp_opify_ASL;
			if ((t = od->c12) != 1.)
				for(i = 0; i < nelq; ++i)
					delsq1[i] *= t;
			objadj = t*objadj + od->c0a;
			for(i = 0, cg = *cgp; cg; cg = cg->next)
				++i;
			ogp = Ograd + co0;
			og2 = i ? (ograd*)M1alloc(i*sizeof(ograd)) : 0;
			for(cg = *cgp; cg; cg = cg->next) {
				*ogp = og = og2++;
				ogp = &og->next;
				og->varno = cg->varno;
				og->coef = t*cg->coef;
				}
			*ogp = 0;
			c = obj_de + co0;
			}
		else if (cgp && objadj != 0.) {
			if (Urhsx) {
				L = LUrhs + co;
				U = Urhsx + co;
				}
			else {
				L = LUrhs + 2*co;
				U = L + 1;
				}
			if (*L > negInfinity)
				*L -= objadj;
			if (*U < Infinity)
				*U -= objadj;
			objadj = 0.;
			}
		en->v = objadj;
		c->e = (expr *)en;
		}
	return nelq;
	}

 void
mqpcheckv_free_ASL(ASL *asl, void **vp)
{
	Memblock *m, *p;
	Static *S;

	if (vp && (S = *(Static**)vp)) {
		for(m = S->mblast; (p = m->prev); m = p)
			free(m);
		if (cterms)
			free(cterms);
		free(S->_s_x);
		*vp = 0;
		}
	}
