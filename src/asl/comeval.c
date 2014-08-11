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

#include "nlp.h"

 void
#ifdef KR_headers
comeval_ASL(asl, i, ie) ASL_fg *asl; int i, ie;
#else
comeval_ASL(ASL_fg *asl, int i, int ie)
#endif
{
	register cexp *c, *ce;
	register expr *e;
	register expr_v *V = var_ex + i;
	register linpart *L, *Le;
	real t;
	c = cexps + i;
	ce = cexps + ie;
	do {
		cv_index = ++i;	/* identify var in case of error */
		e = c->e;
		t = (*e->op)(e C_ASL);
		if (L = c->L)
			for(Le = L + c->nlin; L < Le; L++)
				t += L->fac * *L->v.rp;
		(*V++).v = t;
		}
		while(++c < ce);
	cv_index = 0;
	}

 void
#ifdef KR_headers
com1eval_ASL(asl, i, ie)ASL_fg *asl;  int i, ie;
#else
com1eval_ASL(ASL_fg *asl, int i, int ie)
#endif

{
	register cexp1 *c, *ce;
	register expr *e;
	register expr_v *V = var_ex1 + i;
	register linpart *L, *Le;
	real t;
	c = cexps1 + i;
	ce = cexps1 + ie;
	i += ncom0;
	do {
		cv_index = ncom0 + ++i;	/* identify var in case of error */
		e = c->e;
		t = (*e->op)(e C_ASL);
		if (L = c->L)
			for(Le = L + c->nlin; L < Le; L++)
				t += L->fac * *L->v.rp;
		(*V++).v = t;
		}
		while(++c < ce);
	cv_index = 0;
	}

 void
#ifdef KR_headers
funnelset_ASL(asl, f) ASL_fg *asl; register funnel *f;
#else
funnelset_ASL(ASL_fg *asl, register funnel *f)
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
