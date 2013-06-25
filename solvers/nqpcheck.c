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

#include "nlp.h"
#include "obj_adj.h"
#include "r_qp.hd"

#define GULP		200

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
Static {
	ASL_fg *asl;
	fint *_s_s, *_s_z;
	double *_s_x;
	term *_freeterm, *_term_block;
	ograd *_freeog, *_ograd_block;
	dyad *_freedyad, *_dyad_block;
	int _zerodiv;
	term **_cterms;
	int _dyad_ntogo, nvinc, _ograd_ntogo, _term_ntogo;
	Char **_M1state1, **_M1state2;
	} Static;

#define M1state1	S->_M1state1
#define M1state2	S->_M1state2
#define cterms		S->_cterms
#define dyad_block	S->_dyad_block
#define dyad_ntogo	S->_dyad_ntogo
#define freedyad	S->_freedyad
#define freeog		S->_freeog
#define freeterm	S->_freeterm
#define ograd_block	S->_ograd_block
#define ograd_ntogo	S->_ograd_ntogo
#define s_s		S->_s_s
#define s_x		S->_s_x
#define s_z		S->_s_z
#define term_block	S->_term_block
#define term_ntogo	S->_term_ntogo
#define zerodiv		S->_zerodiv

 static void
free_blocks(Static *S)
{
	M1free_ASL(&S->asl->i, M1state1, M1state2);
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
	else {
		if (!term_ntogo) {
			term_block = (term *)M1alloc_ASL(&S->asl->i,
				GULP*sizeof(term));
			term_ntogo = GULP;
			}
		rv = term_block++;
		--term_ntogo;
		}
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
	else {
		if (!ograd_ntogo) {
			ograd_block = (ograd *)M1alloc_ASL(&S->asl->i,
				GULP*sizeof(ograd));
			ograd_ntogo = GULP;
			}
		rv = ograd_block++;
		--ograd_ntogo;
		}
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
	int i, rv, nz;
	fint *s, *z;
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
	else {
		if (!dyad_ntogo) {
			dyad_block = (dyad *)M1alloc_ASL(&S->asl->i,
				GULP*sizeof(dyad));
			dyad_ntogo = GULP;
			}
		rv = dyad_block++;
		--dyad_ntogo;
		}
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
		if (!(L = cterms[i -= n_var])
		 && !(L = cterms[i] = comterm(S, i)))
			return 0;
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
cterm_free(Static *S, term **cte)
{
	term **ct, *t;
	dyad *d, *d1;

	for(ct = cterms; ct < cte; ct++)
		if ((t = *ct)) {
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
	free(cterms);
	}

 static int
lcmp(const void *a, const void *b, void *v)
{
	Not_Used(v);
	return (int)(*(fint *)a - *(fint *)b);
	}

 fint
mqpcheck_ASL(ASL *a, int co, fint **rowqp, fint **colqp, real **delsqp)
{
	typedef struct dispatch {
		struct dispatch *next;
		fint i, j, jend;
		} dispatch;
	ASL_fg *asl;
	Objrep *od, **pod;
	Static SS, *S;
	cde *c;
	cgrad *cg, **cgp, **cgq, *cq;
	dispatch *cd, *cd0, **cdisp, **cdisp0, *cdnext, **cdp;
	dyad *d, *d1, **q, **q1, **q2, **qe;
	expr *e;
	expr_n *en;
	fint *colq, *colq1, *rowq, *rowq0, *rowq1, *s, *z;
	fint ftn, i, icol, j, ncom, nelq, nv, nz;
	int arrays, *cm, co0, pass, *vmi;
	ograd *og, *og1, *og2, **ogp;
	real *L, *U, *delsq, *delsq0, *delsq1, objadj, t, *x;
	term *T;

	ASL_CHECK(a, ASL_read_fg, "nqpcheck");
	asl = (ASL_fg*)a;
	if (co >= n_obj || co < -n_con)
		return -3L;
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

	memset(S = &SS, 0, sizeof(Static));
	SS.asl = asl;
	if (asl->i.vmap && !asl->i.vminv)
		/* keep vminv from being lost in free_blocks(S) below */
		get_vminv_ASL(a);
	M1state1 = asl->i.Mbnext;
	M1state2 = asl->i.Mblast;
	nv = n_var;
	s_x = x = (double *)Malloc(nv*(sizeof(double)+2*sizeof(fint)));
	s_z = z = (fint *)(x + nv);
	s_s = s = z + nv;
	memset(s, 0, nv*sizeof(fint));
	ftn = Fortran;
	SS.nvinc = nv - asl->i.n_var0 + asl->i.nsufext[ASL_Sufkind_var];

	delsq = delsq0 = delsq1 = 0; /* silence buggy "not-initialized" warning */
	colq = colq1 = rowq = rowq0 = rowq1 = 0;	/* ditto */
	cd0 = 0;					/* ditto */
	cdisp = cdisp0 = 0;				/* ditto */

	if ((ncom = ncom0 + ncom1)) {
		cterms = (term **)Malloc(ncom*sizeof(term*));
		memset(cterms, 0, ncom*sizeof(term*));
		}

	arrays = 1;
	if (rowqp)
		*rowqp = 0;
	else
		arrays = 0;
	if (colqp)
		*colqp = 0;
	else
		arrays = 0;
	if (delsqp)
		*delsqp = 0;
	else
		arrays = 0;

	zerodiv = 0;
	if (!(T = ewalk(S, e)) || zerodiv) {
		free_blocks(S);
		free(x);
		return T ? -2L : -1L;
		}

	if (cterms)
		cterm_free(S, cterms + ncom);
	if (od) {
		od->cg = 0;
		for(i = 0, cg = *cgp; cg; cg = cg->next) {
			if (cg->coef != 0.)
				++i;
			}
		if (i) {
			cgq = &od->cg;
			cq = Malloc(i*sizeof(cgrad));
			for(cg = *cgp; cg; cg = cg->next) {
				*cgq = cq;
				cgq = &cq->next;
				*cq = *cg;
				++cq;
				}
			*cgq = 0;
			}
		}
				
	q = (dyad **)Malloc(nv*sizeof(dyad *));
	qe = q + nv;
	objadj = dsort(S, T, (ograd **)q, cgp, ogp, arrays);

	nelq = nz = 0;
	for(pass = 0; pass < 2; pass++) {
		if (pass) {
			if (!nelq)
				break;
			free(q);
			delsq1 = delsq = (double *)Malloc(nelq*sizeof(real));
			rowq1 = rowq = (fint *)Malloc(nelq*sizeof(fint));
			colq1 = colq = (fint *)Malloc((nv+2)*sizeof(fint));
			nelq = ftn;
			delsq0 = delsq - ftn;
			rowq0 = rowq - ftn;
			q = (dyad **)Malloc(nv*(sizeof(dyad*)
						+ sizeof(dispatch *)
						+ sizeof(dispatch)));
			qe = q + nv;
			cdisp = (dispatch**) qe;
			cdisp0 = cdisp - ftn;
			memset(cdisp, 0, nv*sizeof(dispatch*));
			cd0 = (dispatch *)(cdisp + nv);
			}
		memset(q, 0, nv*sizeof(dyad *));

		if (pass)
			for(d = T->Q; d; d = d->next) {
				og = d->Rq;
				og1 = d->Lq;
				i = og->varno;
				while(og1 && og1->varno < i)
					og1 = og1->next;
				if (og1) {
					q1 = q + i;
					*q1 = new_dyad(S, *q1, og, og1, 0);
					}
				og1 = d->Lq;
				i = og1->varno;
				while(og && og->varno < i)
					og = og->next;
				if (og) {
					q1 = q + i;
					*q1 = new_dyad(S, *q1, og1, og, 0);
					}
				}
		else
			for(d = T->Q; d; d = d->next) {
				og = d->Rq;
				og1 = d->Lq;
				q1 = q + og->varno;
				*q1 = new_dyad(S, *q1, og, og1, 0);
				q1 = q + og1->varno;
				*q1 = new_dyad(S, *q1, og1, og, 0);
				}
		icol = -1;
		vmi = asl->i.vmap ? get_vminv_ASL((ASL*)asl) : 0;
		for(q1 = q; q1 < qe; q1++) {
		    if (pass) {
			*colq++ = nelq;
			for(cd = cdisp[++icol]; cd; cd = cdnext) {
				cdnext = cd->next;
				s[i = cd->i]++;
				x[z[nz++] = i] = delsq0[cd->j++];
				if (cd->j < cd->jend) {
					cdp = cdisp0 + rowq0[cd->j];
					cd->next = *cdp;
					*cdp = cd;
					}
				}
			}
		    if ((d = *q1))
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
					d->next = *q2;
					*q2 = d;
					}
				else {
 d_del:
					free_dyad(S, d);
					}
				}
				while((d = d1));
		if (nz) {
			if (pass) {
				if (nz > 1)
					qsortv(z, nz, sizeof(fint), lcmp, NULL);
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
		    }
		}
	free(q);
	free_blocks(S);
	free(x);
	if (od && od->cg)
		M1record(od->cg);
	if (nelq) {
		if (arrays) {
			/* allow one more for obj. adjustment */
			*colq = colq[1] = nelq;
			*rowqp = rowq1;
			*colqp = colq1;
			*delsqp = delsq1;
			}
		else {
			free(colq1);
			free(rowq1);
			free(delsq1);
			}
		nelq -= ftn;
		}
	if (arrays) {
		en = (expr_n *)mem(sizeof(expr_n));
		en->op = f_OPNUM_ASL;
		if (od) {
			od->opify = 1;
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
			objadj = 0;
			}
		en->v = objadj;
		c->e = (expr *)en;
		}
	return nelq;
	}


 fint
nqpcheck_ASL(ASL *asl, int co, fint **rowqp, fint **colqp, real **delsqp)
{
	fint rv = mqpcheck_ASL(asl, co, rowqp, colqp, delsqp);
	if (rowqp && *rowqp) {
		M1record(*delsqp);
		M1record(*rowqp);
		M1record(*colqp);
		}
	return rv;
	}
