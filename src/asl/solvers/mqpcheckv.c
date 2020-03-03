/*******************************************************************
Copyright (C) 2018 AMPL Optimization, Inc.; written by David M. Gay.

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
	Memblock *mb[2], *mb0, *mblast[2];
	dispatch *cd0, **cdisp;
	double *_s_x;
	dyad *_freedyad[2], **_s_q;
	int *_s_s, *_s_z, *w;
	ograd *_freeog[2], **oq;
	term **_cterms, *_freeterm[2];
	void **v[2], **ve[2];
	int _zerodiv, comterms, flags, nvinc;
	} Static;

#define cterms		S->_cterms
#define freedyad	S->_freedyad[S->comterms]
#define freeog		S->_freeog[S->comterms]
#define freeterm	S->_freeterm[S->comterms]
#define s_q		S->_s_q
#define s_s		S->_s_s
#define s_x		S->_s_x
#define s_z		S->_s_z
#define zerodiv		S->_zerodiv

 static void*
M2alloc(Static *S, size_t len, int dbl_align)
{
	Memblock *mb, *mb1;
	int ct;
	size_t n;
	void *rv, ***vp;

	vp = &S->v[ct = S->comterms];
	if (dbl_align)
		*vp = (void**)(((size_t)*vp + (sizeof(void*)-1)) & ~(sizeof(void*)-1));
	n = (len + (sizeof(void*)-1))/sizeof(void*);
	if (*vp + n >= S->ve[ct]) {
		mb = S->mb[ct];
		if (!mb || !(mb1 = mb->next)) {
			S->mblast[ct] = mb1 = (Memblock*)Malloc(sizeof(Memblock));
			if (mb)
				mb->next = mb1;
			mb1->prev = mb;
			mb1->next = 0;
			}
		S->mb[ct] = mb1;
		S->v[ct] = mb1->x;
		S->ve[ct] = (S->v[ct] = mb1->x) + Memblock_gulp;
		vp = &S->v[ct];
		}
	rv = (void*)*vp;
	*vp += n;
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
	ASL_fg* asl;
	cexp *c;
	cexp1 *c1;
	expr_v ev, *vp;
	int nlin;
	linpart *L, *Le;
	ograd *og;
	term *T;

	asl = S->asl;
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
	if (L && T) {
		for(Le = L + nlin; L < Le; L++) {
			vp = (expr_v *)((char *)L->v.rp -
				((char *)&ev.v - (char *)&ev));
			T = termsum(S, T, new_term(S,
				new_og(S, 0, (int)(vp - var_e), L->fac)));
			}
		}
	return T;
	}

 static void
Comeval(Static *S, int i, int ie)
{
	term **Cterms;

	S->comterms = 1;
	for(Cterms = cterms; i < ie; ++i)
		Cterms[i] = comterm(S, i);
	S->comterms = 0;
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
	ASL_fg *asl;
	expr **ep, **epe;
	int i;
	ograd *o, *oR;
	term *L, *R, *T;

	switch(Intcast e->op) {

	  case OPNUM:
		return new_term(S, new_og(S, 0, -1, ((expr_n *)e)->v));

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
			/* c_cexp1st and o_cexp1st may not have been allocated */
			S->comterms = 1;
			if (!(L = comterm(S, i)))
				return 0;
			cterms[i] = L;
			S->comterms = 0;
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
	int *C1, *cm, *colno, *qm, *rowq, *rowq0, *rowq1, *s, *vmi, *w, *z;
	int C10, arrays, co0, dv0, dv1, dvbit, ftn, i, icol, icolf, j;
	int ncol, ncom, nv, nva, nz, nz1, pass;
	ograd *og, *og1, *og2, **ogp;
	real *L, *U, *delsq, *delsq0, *delsq1, objadj, t, *x;
	size_t  *colq, *colq1, nelq, nelq0;
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
			if (!(cgp = asl->i.Cgrad0))
				cgp = Cgrad;
			goto use_Cgrad;
			}
		dv0 = combc;
		dv1 = dv0 + como;
		dvbit = 2;
		C1 = o_cexp1st;
		c = obj_de + co;
		ogp = Ograd + co;
		cgp = 0;
		}
	else {
		co = -1 - co;
		if ((cm = asl->i.cmap))
			co = cm[co];
		cgp = Cgrad;
 use_Cgrad:
		dv0 = comb;
		dv1 = combc;
		dvbit = 3;
		C1 = c_cexp1st;
		c = con_de + co;
		cgp += co;
		ogp = 0;
		}

	e = c->e;
	if (e->op == f_OPNUM)
		return 0;

	if (asl->i.vmap && !asl->i.vminv)
		get_vminv_ASL(a);
	nv = n_var;
	C10 = ncom0;
	ncom = C10 + ncom1;
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
		S->mb0 = S->mblast[0] = mb;
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
			cterms = (term **)Malloc(ncom*sizeof(term*));
			memset(cterms, 0, ncom*sizeof(term*));
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
	S->mb[0] = mb = S->mb0;
	S->v[0]  = &mb->x[0];
	S->ve[0] = &mb->x[Memblock_gulp];
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
	if (comb && !(S->flags & 1)) {
		Comeval(S, 0, comb);
		S->flags |= 1;
		}
	if (!(S->flags & dvbit)) {
		S->flags |= dvbit;
		if (dv1 > dv0)
			Comeval(S, dv0, dv1);
		}
	if (C1 && C1[co] < C1[co+1])
		Comeval(S, C10 + C1[co], C10 + C1[co+1]);
	if (!(T = ewalk(S, e)) || zerodiv)
		return T ? -2L : -1L;

	if (od) {
		cgq = &od->cg;
		for(i = 0, cg = *cgp; cg; cg = cg->next) {
			if (cg->coef != 0.)
				++i;
			}
		if (i) {
			cq = (cgrad*)M1alloc(i*sizeof(cgrad));
			for(cg = *cgp; cg; cg = cg->next) {
				if (cg->coef != 0.) {
					*cgq = cq;
					cgq = &cq->next;
					*cq = *cg;
					++cq;
					}
				}
			}
		*cgq = 0;
		}

	objadj = dsort(S, T, S->oq, cgp, ogp, arrays);

	icolf = nelq = ncol = nz = nz1 = 0;
	qpi = 0;
	/* In pass 0, we the count nonzeros in the lower triangle. */
	/* In pass 1, we compute the lower triangle and use column dispatch */
	/* (via the cdisp array) to copy the strict lower triangle to the */
	/* strict upper triangle.  This ensures symmetry. */
	for(pass = 0; pass < 2; pass++) {
		if (pass) {
			if (!nelq)
				break;
			nelq += nelq - nz1; /* nz1 = number of diagonal elements */
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
			qpi->delsq = delsq = delsq1 = (double *)(qpi+1);
			qpi->colbeg = colq = (size_t *)(delsq + nelq);
			qpi->rowno = rowq = (int *)(colq + ncol + 1);
			qpi->colno = colno = rowq + nelq;
			qpi->nc = ncol;
			qpi->nz = nelq;
			nelq = ftn;
			delsq0 = delsq - ftn;
			rowq0 = rowq - ftn;
			}
		for(d = T->Q; d; d = d->next) {
			og = d->Rq;
			og1 = d->Lq;
			i = og->varno;
			while(og1 && og1->varno < i)
				og1 = og1->next;
			if (og1) {
				q1 = q + i;
				if (!w[i]) {
					w[i] = 1;
					AVL_vinsert(AQ, 0, (Element*)&w[i], 0);
					}
				*q1 = new_dyad(S, *q1, og, og1, 0);
				}
			og1 = d->Lq;
			i = og1->varno;
			while(og && og->varno < i)
				og = og->next;
			if (og) {
				q1 = q + i;
				if (!w[i]) {
					w[i] = 1;
					AVL_vinsert(AQ, 0, (Element*)&w[i], 0);
					}
				*q1 = new_dyad(S, *q1, og1, og, 0);
				}
			}
		vmi = asl->i.vmap ? get_vminv_ASL((ASL*)asl) : 0;
		for(qm = (int*)AVL_first(AQ, &NQ); qm; ) {
			NQ0 = NQ;
			icol = qm - w;
			nelq0 = nelq;
			if (pass) {
				*qm = 0;
				icolf = icol + ftn;
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
					goto get_d1;
					}
				  d->Rq = og2;
				  }
 get_d1:
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
					if (!w[i]) {
						w[i] = 1;
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
					for(i = nz1 = 0; i < nz; i++) {
						if ((t = x[j = z[i]])) {
							*delsq++ = t;
							if (vmi)
								j = vmi[j];
							*rowq++ = j + ftn;
							nelq++;
							z[nz1++] = j;
							}
						s[j] = 0;
						}
					if (nelq > nelq0) {
						*colq++ = nelq0;
						*colno++ = icolf;
						}
					for(i = 0; i < nz1; i++)
					    if ((j = z[i]) > icol) {
						cd0->i = icol;
						cd0->j = nelq0 + i;
						cd0->jend = nelq;
						cdp = cdisp + j;
						cd0->next = *cdp;
						*cdp = cd0++;
						break;
						}
					nz = 0;
					}
				else {
					while(nz > 0) {
						s[i = z[--nz]] = 0;
						if (x[i]) {
							++nelq;
							if (i == icol)
								++nz1;
							else {
								if (!w[i])
						AVL_vinsert(AQ, 0, (Element*)&w[i], 0);
								w[i] = 2;
								}
							}
						}
					if (nelq > nelq0 || w[icol] == 2)
						++ncol;
					}
				}
			else if (!pass && w[icol] == 2)
				++ncol;
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
	Memblock *m, *m0, *p;
	Static *S;
	int i;

	if (vp && (S = *(Static**)vp)) {
		m0 = S->mb0;
		for(i = 0; i < 2; ++i, m0 = 0) {
			for(m = S->mblast[i]; m != m0; m = p) {
				p = m->prev;
				free(m);
				}
			}
		if (cterms)
			free(cterms);
		AVL_Tree_free(&S->AQ);
		free(S->_s_x);
		*vp = 0;
		}
	}
