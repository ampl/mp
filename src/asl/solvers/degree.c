/*******************************************************************
Copyright (C) 2016 AMPL Optimization, Inc.; written by David M. Gay.

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

/* After qp_read() or fg_wread(), degree(asl, co, pv) determines whether
   the indicated objective (co >= 0, 0 <= co < n_obj) or constraint
   (co < 0 indicates constraint 1 - co, 0 <= 1-co < n_con) is linear,
   quadratic, or higher-order.  In contrast to nqpcheck() and variants,
   this routine does not attempt to compute the number of nonzeros in
   the Hessian, so it can run much faster.  Possible return values:

	-1 ==> invalid co value
	 0 ==> constant
	 1 ==> linear
	 2 ==> quadratic
	 3 ==> more general nonlinear

   The pv argument (of type void **pv) should be NULL if degree_ASL(...)
   is to be called just once (e.g., to see if the objective is quadratic).
   If multiple calls are expected (e.g., for objectives and constraints),
   it may save time to use the pattern

	void *v = 0;
	for(...) { ... degree(asl, co, &v); ... }
	if (v) free(v);
*/

#include "nlp.h"
#include "obj_adj.h"
#include "r_qp.hd"

 typedef struct
Dhelp {
	cexp *c0;
	cexp1 *c1;
	expr_v *vare;
	int *vk, nc0, nv;
	} Dhelp;

 static int
kind(expr *e, Dhelp *h)
{
	cexp *c0;
	cexp1 *c1;
	expr **ep, **epe;
	int i, j;

 top:
	switch(Intcast e->op) {

	  case OPNUM:
		return 0;

	  case OPPLUS:
	  case OPMINUS:
		i = kind(e->L.e, h);
		if (i < 3) {
			j = kind(e->R.e, h);
			if (i < j)
				i = j;
			}
		return i;

	  case OPUMINUS:
		e = e->L.e;
		goto top;

	  case OPMULT:
		i = kind(e->L.e, h);
		if (i < 3)
			i += kind(e->R.e, h);
		return i;

	  case OPDIV:
		i = kind(e->L.e, h);
		if (i <= 2) {
			j = kind(e->R.e, h);
			if (j > 0)
				i = 3;
			}
		return i;

	  case OPSUMLIST:
		ep = e->L.ep;
		epe = e->R.ep;
		i = kind(*ep, h);
		while(i < 3 && ++ep < epe) {
			j = kind(*ep, h);
			if (i < j)
				i = j;
			}
		return i;

	  case OP2POW:
		i = kind(e->L.e, h);
		if (i > 1)
			i = 3;
		else
			i <<= 1;
		return i;

	  case OPVARVAL:
		i = (expr_v *)e - h->vare;
		if ((i -= h->nv) < 0)
			return 1;
		if ((j = h->vk[i]) == -2) {
			if (i >= h->nc0) {
				c1 = &h->c1[i-h->nc0];
				if (!(j = kind(c1->e, h)) && c1->nlin)
					j = 1;
				}
			else {
				c0 = &h->c0[i];
				if (!(j = kind(c0->e, h)) && c0->nlin)
					j = 1;
				}
			h->vk[i] = j;
			}
		return j;
		}
	return 3; /* nonlinear */
	}

 int
degree_ASL(ASL *a, int co, void **pv)
{
	ASL_fg *asl;
	Dhelp h;
	Objrep *od, **pod;
	cde *c;
	cgrad *cg;
	int *cm, i, ncom, rv;
	ograd *og;

	ASL_CHECK(a, ASL_read_fg, "degree");
	asl = (ASL_fg*)a;
	if (co >= n_obj || co < -n_con)
		return -1;
	h.nv = n_var;
	h.vare = var_e;
	h.vk = 0;
	h.nc0 = ncom0;
	if ((ncom = h.nc0 + ncom1)) {
		h.c0 = cexps;
		h.c1 = cexps1;
		if (!pv || !(h.vk = *(int**)pv)) {
			h.vk = (int*)Malloc(ncom*sizeof(int));
			for(i = 0; i < ncom; ++i)
				h.vk[i] = -2;
			if (pv)
				*pv = h.vk;
			}
		}
	od = 0;
	if (co >= 0) {
		if ((pod = asl->i.Or) && (od = pod[co])) {
			co = od->ico;
			goto use_Cgrad;
			}
		c = obj_de + co;
		og = Ograd[co];
		cg = 0;
		}
	else {
		co = -1 - co;
		if ((cm = asl->i.cmap))
			co = cm[co];
 use_Cgrad:
		c = con_de + co;
		cg = Cgrad[co];
		og = 0;
		}
	rv = kind(c->e, &h);
	if (h.vk && !pv)
		free(h.vk);
	if (rv > 3)
		rv = 3;
	else if (rv == 0) {
		while(og) {
			if (og->coef) {
				rv = 1;
				goto ret;
				}
			og = og->next;
			}
		while(cg) {
			if (cg->coef) {
				rv = 1;
				goto ret;
				}
			cg = cg->next;
			}
		}
 ret:
	return rv;
	}
