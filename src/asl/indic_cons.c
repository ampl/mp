/****************************************************************
Copyright (C) 2011 AMPL Optimization LLC; written by David M. Gay.

Permission to use, copy, modify, and distribute this software and its
documentation for any purpose and without fee is hereby granted,
provided that the above copyright notice appear in all copies and that
both that the copyright notice and this permission notice and warranty
disclaimer appear in supporting documentation.

The author and AMPL Optimization LLC disclaim all warranties with
regard to this software, including all implied warranties of
merchantability and fitness.  In no event shall the author be liable
for any special, indirect or consequential damages or any damages
whatsoever resulting from loss of use, data or profits, whether in an
action of contract, negligence or other tortious action, arising out
of or in connection with the use or performance of this software.
****************************************************************/

#include "nlp.h"
#include "opcode.hd"

 static int Neg[6] = {2,3,0,1,5,4};

#ifndef LCADJ_GULP
#define LCADJ_GULP 1023
#endif
 enum { Gulp = LCADJ_GULP };

 typedef struct
Lconstr {
	real LU[2];
	struct Lconstr *next;
	real *x;
	int *z;
	int nx;
	} Lconstr;

 typedef union Tchunk { union Tchunk *prev; real r; } Tchunk;

 typedef struct
LCADJ_Info {
	ASL *asl;
	Lconstr *lc[2];	/* "true" (==>) constraints and "false" (else) constraints */
	Lconstr **plc[2];
	ograd *freeog, **s;
	int *z;
	real rhs;
	real *tfree, *tfree0;
	Tchunk *tchunks;	/* temp. memory (unlikly to be needed; freed after use) */
	int b[2];
	int nlv[6];
	int ntfree;	/* elements remaining in tfree */
	int n1lc;	/* numbers of reals in one Lconstr */
	int n1og;	/* numbers of reals in one ograd */
	int zerodiv;
	int vk, vnb;
	void *v;
	int *errinfo;
	Add_Indicator add_indic;
	} LCADJ_Info;

 static int
bincheck(LCADJ_Info *lci, expr *e)
{
	expr_n *en;
	expr_v *v;
	int j, k;

	switch(Intcast e->op) {
	  case GT: k = 1; break;
	  case GE: k = 2; break;
	  case LE: k = 3; break;
	  case LT: k = 4; break;
	  case EQ: k = 5; break;
	  case NE: k = 6; break;
	  default: return 0;
	  }
	v = (expr_v*)e->L.e;
	switch (Intcast v->op) {
	  case OPVARVAL:
		en = (expr_n*) e->R.e;
		if ((Intcast en->op) == OPNUM) {
 hit:
			if ((j = v->a) < lci->nlv[0] || (j >= lci->nlv[1] && j < lci->nlv[2])
			 || (j >= lci->nlv[3] && j < lci->nlv[4]) || j >= lci->nlv[5])
				return 0;
			lci->vk = j;
			lci->rhs = en->v;
			return k;
			}
		break;
	  case OPNUM:
		en = (expr_n*)v;
		v = (expr_v*)e->R.e;
		if (k < 5)
			k = 5 - k;
		if ((Intcast v->op) == OPVARVAL)
			goto hit;
	  }
	return 0;
	}

 static void
new_tchunk(LCADJ_Info *lci, int tneed)
{
	Tchunk *tc;

	if (tneed < Gulp)
		tneed = Gulp;
	tc = (Tchunk*)Malloc((tneed+1)*sizeof(real));
	tc->prev = lci->tchunks;
	lci->tchunks = tc;
	lci->tfree = &tc[1].r;
	lci->ntfree = tneed;
	}

 static real*
tmem(LCADJ_Info *lci, size_t L)
{
	int n = (L + sizeof(real) - 1) / sizeof(real);
	real *r;

	if (n > lci->ntfree)
		new_tchunk(lci, n);
	r = lci->tfree;
	lci->tfree = r + n;
	lci->ntfree -= n;
	return r;
	}

 static ograd *
new_og(LCADJ_Info *lci, int varno, real coef)
{
	int n1;
	ograd *og;

	if ((og = lci->freeog))
		lci->freeog = og->next;
	else {
		if (lci->ntfree < (n1 = lci->n1og))
			new_tchunk(lci, n1);
		og = (ograd*)lci->tfree;
		lci->tfree += n1;
		lci->ntfree -= n1;
		}
	og->next = 0;
	og->varno = varno;
	og->coef = coef;
	return og;
	}

 static void
free_og(LCADJ_Info *lci, ograd *og1, ograd *og1e)
{
	og1e->next = lci->freeog;
	lci->freeog = og1;
	}

#define voffset_of(t,c) ((void *)&((t*)0)->c)

 static ograd *
finish_plus(LCADJ_Info *lci, ograd *og1, ograd *og2, ograd **oglp)
{
	ograd *og, *og0, *oge, **ogp;
	ssize_t d;

	og0 = oge = 0;
	ogp = &og0;
	for(;;) {
		d = og1->varno - og2->varno;
		if (d < 0) {
			*ogp = og1;
			oge = og1;
			if ((og1 = og1->next)) {
				ogp = &oge->next;
				continue;
				}
			oge->next = og2;
 finish_og2:
			while((og = og2->next))
				og2 = og;
			*oglp = og2;
			break;
			}
		if (d > 0) {
			*ogp = og2;
			oge = og2;
			if ((og2 = og2->next)) {
				ogp = &oge->next;
				continue;
				}
			oge->next = og1;
 finish_og1:
			while((og = og1->next))
				og1 = og;
			*oglp = og1;
			break;
			}
		og1->coef += og2->coef;
		og = og2->next;
		og2->next = lci->freeog;
		lci->freeog = og2;
		og2 = og;
		if (og1->coef != 0.) {
			*ogp = oge = og1;
			ogp = &og1->next;
			og1 = *ogp;
			}
		else {
			og = og1->next;
			og1->next = lci->freeog;
			lci->freeog = og1;
			og1 = og;
			}
		if (!og1) {
			if ((*ogp = og2))
				goto finish_og2;
			*oglp = oge;
			break;
			}
		if (!og2) {
			*ogp = og1;
			goto finish_og1;
			}
		}
	if (!og0)
		og0 = *oglp = new_og(lci, -1, 0.);
	return og0;
	}

 static ograd *
linform(LCADJ_Info *lci, expr *e, ograd **oglp)
{
	ASL_fg *asl;
	cexp *ce;
	cexp1 *ce1;
	expr **ep, **epe;
	expr_n *en;
	int i;
	linpart *L, *Le;
	ograd *og, *og1, *og1e, *og2, *og2e, *og2x, **ogp;
	real t;

	switch(Intcast e->op) {

	  case OPNUM:
		return *oglp = new_og(lci, -1, ((expr_n *)e)->v);

	  case OPPLUS:
		if (!(og1 = linform(lci, e->L.e, &og1e)))
			return og1;
		if (!(og2 = linform(lci, e->R.e, &og2e))) {
			free_og(lci, og1, og1e);
			return og2;
			}
		return finish_plus(lci, og1, og2, oglp);

	  case OPMINUS:
		if (!(og1 = linform(lci, e->L.e, &og1e)))
			return og1;
		if (!(og2 = linform(lci, e->R.e, &og2e))) {
			free_og(lci, og1, og1e);
			return og2;
			}
		for(og = og2; og; og = og->next)
			og->coef = -og->coef;
		return finish_plus(lci, og1, og2, oglp);

	  case OPUMINUS:
		if ((og1 = linform(lci, e->L.e, oglp))) {
			og2 = og1;
			do og2->coef = -og2->coef;
				while((og2 = og2->next));
			}
		return og1;

	  case OPMULT:
		if (!(og1 = linform(lci, e->L.e, &og1e)))
			return og1;
		if (!(og2 = linform(lci, e->R.e, &og2e))) {
			free_og(lci, og1, og1e);
			return og2;
			}
		if (og1->varno < 0 && !og1->next) {
			t = og1->coef;
			free_og(lci, og1, og1e);
			}
		else if (og2->varno < 0 && !og2->next) {
			t = og2->coef;
			free_og(lci, og2, og2e);
			og2 = og1;
			og2e = og1e;
			}
		else {
			free_og(lci, og1, og1e);
			free_og(lci, og2, og2e);
			return 0;
			}
		for(og = og2; og; og = og->next)
			og->coef *= t;
		*oglp = og2e;
		return og2;

	  case OPDIV:
		/* only allow division by a constant */
		if (!(og1 = linform(lci, e->L.e, &og1e)))
			return og1;
		if (!(og2 = linform(lci, e->L.e, &og2e))) {
			free_og(lci, og1, og1e);
			return og2;
			}
		if (og2->varno < 0 && !og2->next) {
			t = og2->coef;
			free_og(lci, og2, og2e);
			}
		else {
			free_og(lci, og1, og1e);
			free_og(lci, og2, og2e);
			return 0;
			}
		for(og = og1; og; og = og->next)
			og->coef /= t;
		*oglp = og1e;
		return og1;

	  case OPSUMLIST:
		ep = e->L.ep;
		epe = e->R.ep;
		if (!(og1 = linform(lci, *ep, &og1e)))
			return og1;
		while(++ep < epe) {
			if (!(og2 = linform(lci, *ep, &og2e))) {
				free_og(lci, og1, og1e);
				return og2;
				}
			if (og1->varno > og2->varno) {
				og = og1;
				og1 = og2;
				og2 = og;
				og = og1e;
				og1e = og2e;
				og2e = og;
				}
			else for(og = og1; og; og = og->next) {
				if (!og2) {
					og2e = og1e;
					break;
					}
				if (og2->varno != og->varno)
					break;
				og->coef += og2->coef;
				og2x = og2->next;
				og2->next = lci->freeog;
				lci->freeog = og2;
				og2 = og2x;
				}
			og1e->next = og2;
			og1e = og2e;
			}
		*oglp = og1e;
		return og1;

	  case OPVARVAL:
		asl = (ASL_fg*)lci->asl;
		if ((i = (expr_v *)e - var_e) < n_var)
			return *oglp = new_og(lci, i, 1.);
		if ((i -= n_var) < ncom0) {
			ce = cexps + i;
			en = (expr_n*)ce->e;
			L =  ce->L;
			Le = L + ce->nlin;
			}
		else {
			ce1 = cexps1 + (i - ncom0);
			en = (expr_n*)ce1->e;
			L = ce1->L;
			Le = L + ce1->nlin;
			}
		if ((Intcast en->op) != OPNUM)
			return 0;
		ogp = &og2;
		if (en->v != 0.) {
			og2 = new_og(lci, -1, en->v);
			ogp = &og2->next;
			}
		for(og = 0; L < Le; L++) {
			i = (expr_v*)((char*)L->v.rp - (char*)voffset_of(expr_v,v)) - var_e;
			og = *ogp = new_og(lci, i, L->fac);
			ogp = &og->next;
			}
		*ogp = 0;
		*oglp = og;
		return og2;
	  }
	return 0;
	}

 static int
intcomp(const void *a, const void *b, void *v)
{
	return *(int*)a - *(int*)b;
	}

 static int
lincheck(LCADJ_Info *lci, expr *e, Lconstr ***pplc)
{
	Lconstr *lc;
	expr **ep1, **ep2;
	int i, j, k, m, *y, *z;
	ograd *freeog, *og, *og1, *og1e, *og2, *og2e, **s;
	real *LU, rhs, *x;

 top:
	switch(Intcast e->op) {
	  case OPAND:
		if ((k = lincheck(lci, e->L.e, pplc)))
			return k;
		e = e->R.e;
		goto top;
	  case ANDLIST:
		ep1 = e->L.ep;
		for(ep2 = e->R.ep; ep1 < ep2; ep1++)
			if ((k = lincheck(lci, *ep1, pplc)))
				return k;
		return 0;
	  case GE: k = 2; break;
	  case LE: k = 3; break;
	  case EQ: k = 5; break;
	  default: return 1;
	  }
	if (!(og1 = linform(lci, e->L.e, &og1e)))
		return 1;
	if (!(og2 = linform(lci, e->R.e, &og2e))) {
		free_og(lci, og1, og1e);
		return 1;
		}
	for(og = og2; og; og = og->next)
		og->coef = -og->coef;
	rhs = 0.;
	og1 = finish_plus(lci, og1, og2, &og1e);
	if (og1 && og1->varno < 0) {
		rhs = -og1->coef;
		og2 = og1->next;
		og1->next = lci->freeog;
		lci->freeog = og1;
		og1 = og2;
		}
	s = lci->s;
	z = lci->z;
	m = 0;
	freeog = lci->freeog;
	while((og = og1)) {
		og1 = og->next;
		if ((i = og->varno) < 0)
			rhs -= og->coef;
		else if (!(og2 = s[i])) {
			s[i] = og;
			z[m++] = i;
			continue;
			}
		else
			og2->coef += og->coef;
		og->next = freeog;
		freeog = og;
		}
	if (m > 1)
		qsortv(z, m, sizeof(int), intcomp, 0);
	x = tmem(lci, m*(sizeof(real) + sizeof(int)) + sizeof(Lconstr));
	**pplc = lc = (Lconstr*)(x + m);
	y = (int*)(lc + 1);
	lc->next = 0;
	*pplc = &lc->next;
	lc->x = x;
	lc->z = y;
	lc->nx = m;
	LU = lc->LU;
	for(i = 0; i < m; i++) {
		y[i] = j = z[i];
		og = s[j];
		s[j] = 0;
		x[i] = og->coef;
		og->next = freeog;
		freeog = og;
		}
	lci->freeog = freeog;
	LU[0] = LU[1] = rhs;
	switch(k) {
	  case 2:
		LU[1] = Infinity;
		break;
	  case 3:
		LU[0] = negInfinity;
	  }
	return 0;
	}

 static int
lcadj(LCADJ_Info *lci, expr *e)
{
	expr *e1, *e2, *e3;
	expr_if *eif;
	int b1, b2, i, k, vk;
	real rhs;

	e3 = 0;	/* silence buggy warning about possibly not being initialized */
	switch(Intcast e->op) {
	  case OPOR:
		e1 = e->L.e;
		e2 = e->R.e;
		rhs = 0.; /* silence buggy warning about possibly not being initialized */
		vk = 0; /* ditto */
		if ((b1 = bincheck(lci,e1))) {
			vk = lci->vk;
			rhs = lci->rhs;
			}
		if ((b2 = bincheck(lci,e2))) {
			if (b1 && !lincheck(lci, e2, &lci->plc[0])) {
				lci->b[0] = Neg[b1 - 1];
				lci->vk = vk;
				lci->rhs = rhs;
				return 1;
				}
			b1 = b2;
			e1 = e2;
			e2 = e->L.e;
			}
		else if (!b1)
			return 0;
		k = 1;
		break;
	  case OPIMPELSE:
		eif = (expr_if*)e;
		e1 = eif->e;
		if (!(b1 = bincheck(lci,eif->e)))
			return 0;
		k = 2;
		e2 = eif->F;
		e3 = eif->T;
		break;
	  default:
		return 0;
	  }
	b1 = Neg[b1 - 1];
	for(i = 0; i < k; i++, e2 = e3) {
		/* must split e2 into const op elin, pass elin to xqpcheck */
		if (lincheck(lci, e2, &lci->plc[i]))
			return 0;
		lci->b[i] = b1;
		b1 = Neg[b1];
		}
	return k;
	}

 static void
chunkfree(LCADJ_Info *lci)
{
	Tchunk *tc, *tc1;
	for(tc1 = lci->tchunks; (tc = tc1); ) {
		tc1 = tc->prev;
		free(tc);
		}
	lci->tchunks = 0;
	}

 static int
add_indicator(int ci, LCADJ_Info *lci, expr *e)
{
	Lconstr *lc;
	int compl, i, j, k, sense, vk;
	real *LU, rhs;

	if (lci->tchunks)
		chunkfree(lci);
	lci->freeog = 0;
	lci->tfree = lci->tfree0;
	lci->ntfree = Gulp;
	for(i = 0; i < 2; ++i)
		*(lci->plc[i] = &lci->lc[i]) = 0;
	k = lcadj(lci, e);
	if (!k) {
		lci->errinfo[0] = ci;
		return 1;
		}
	vk = lci->vk;
	rhs = lci->rhs;
	compl = 0;
	switch(lci->b[0]) {
	  case 0: /* > */
		if (rhs < 0. || rhs >= 1.) {
 bad_cmpop:
			lci->errinfo[0] = ci;
			lci->errinfo[1] = vk;
			return 2;
			}
		break;
	  case 1: /* >= */
		if (rhs <= 0. || rhs > 1.)
			goto bad_cmpop;
		break;
	  case 2: /* <= */
		if (rhs < 0. || rhs >= 1.)
			goto bad_cmpop;
		compl = 1;
		break;
	  case 3: /* < */
		if (rhs <= 0. || rhs > 1.)
			goto bad_cmpop;
		compl = 1;
		break;
	  case 4: /* == */
		if (rhs == 1.)
			break;
		if (rhs == 0.) {
			compl = 1;
			break;
			}
		goto bad_cmpop;
	  case 5: /* != */
		if (rhs == 0.)
			break;
		if (rhs == 1.) {
			compl = 1;
			break;
			}
		goto bad_cmpop;
	  }
	for(j = 0; j < k; j++, compl = 1 - compl) {
		for(lc = lci->lc[j]; lc; lc = lc->next) {
			LU = lc->LU;
			if (LU[0] > negInfinity) {
				rhs = LU[0];
				sense = LU[1] == rhs ? 2: 1;
				}
			else {
				rhs = LU[1];
				sense = 0;
				}
			if ((i = lci->add_indic(lci->v, vk, compl, sense, lc->nx,
					lc->z, lc->x, rhs))) {
 badret:
				lci->errinfo[0] = i;
				return 3;
				}
			if (sense == 1 && LU[1] < Infinity
			 && (i = lci->add_indic(lci->v, vk, compl, sense, lc->nx,
					lc->z, lc->x, rhs)))
				goto badret;
			}
		}
	return 0;
	}

 int
indicator_constrs_ASL(ASL *asl, void *v, Add_Indicator add_indic, int errinfo[2])
{
	LCADJ_Info lci;
	cde *logc;
	int i, n, nlogc, rc;
	real chunk1[Gulp];
	static char who[] = "indicator_constrs_ASL";

	ASL_CHECK(asl, ASL_read_fg, who);
	if (!(nlogc = n_lcon))
		return 0;
	memset(&lci, 0, sizeof(lci));
	lci.v = v;
	lci.tfree0 = chunk1;
	n = n_var;
	lci.s = (ograd**)Malloc(n*(sizeof(int) + sizeof(ograd*)));
	lci.z = (int*)(lci.s + n);
	memset(lci.s, 0, n*sizeof(ograd*));
	lci.n1lc = (sizeof(Lconstr) + sizeof(real) - 1) / sizeof(real);
	lci.n1og = (sizeof(ograd) + sizeof(real) - 1) / sizeof(real);
	lci.asl = asl;
	lci.nlv[1] = i = nlvb;
	lci.nlv[0] = i - nlvbi;
	lci.nlv[3] = i = nlvc;
	lci.nlv[2] = i - nlvci;
	if (nlvo > nlvc)
		i = nlvo;
	lci.nlv[5] = i;
	lci.nlv[4] = i - nlvoi;
	lci.errinfo = errinfo;
	lci.add_indic = add_indic;
	logc = ((ASL_fg*)asl)->I.lcon_de_;
	for(i = rc = 0; i < nlogc; ++i)
		if ((rc = add_indicator(i, &lci, logc[i].e)))
			break;
	if (lci.tchunks)
		chunkfree(&lci);
	free(lci.s);
	return rc;
	}
