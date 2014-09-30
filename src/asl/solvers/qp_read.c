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

#include "nlp.h"
#include "r_opn.hd"	/* for N_OPS */

extern efunc *r_ops[];

 int
#ifdef KR_headers
qp_read_ASL(asl, nl, flags) ASL *asl; FILE *nl; int flags;
#else
qp_read_ASL(ASL *asl, FILE *nl, int flags)
#endif
{
	efunc *rops[N_OPS];
	int i;

	ASL_CHECK(asl, ASL_read_fg, "edqpread");
	for(i = 0; i < N_OPS; i++)
		rops[i] = (efunc*)(size_t)i;
	((ASL_fg*)asl)->I.r_ops_ = rops;
	i = fg_read_ASL(asl, nl, flags);
	((ASL_fg*)asl)->I.r_ops_ = 0;
	return i;
	}

 static void
#ifdef KR_headers
ed1opwalk(e, f, varval) expr *e; efunc **f, *varval;
#else
ed1opwalk(expr *e, efunc **f, efunc *varval)
#endif
{
	int n;
	argpair *ap, *ape;
	de *d;
	efunc *op;
	expr **ep, **epe;
	expr_f *ef;
	expr_if *eif;

 top:
	if ((op = e->op) == varval || op == (efunc*)f_OPNUM_ASL)
		return;
	n = Intcast op;
	if (n < 0 || n >= N_OPS) {
		fprintf(Stderr, "qp_opify: bad op field\n");
		exit(1);
		}
	e->op = f[n];
	switch(optypeb[n]) {

	  case 2: /* binary */
		ed1opwalk(e->R.e, f, varval);
		/* no break */

	  case 1: /* unary */
		e = e->L.e;
		goto top;

	  case 3: /* vararg (min, max) */
		for(d = ((expr_va*)e)->L.d; d->e; d++)
			ed1opwalk(d->e, f, varval);

	  case 4: /* piece-wise linear */
	  case 8:
	  case 9:
	  case 10:
	  case 11:
		break;

	  case 5: /* if */
		eif = (expr_if *)e;
		ed1opwalk(eif->T, f, varval);
		ed1opwalk(eif->F, f, varval);
		e = eif->e;
		goto top;

	  case 6: /* sumlist */
		ep = e->L.ep;
		epe = e->R.ep;
		do ed1opwalk(*ep++, f, varval);
			while(ep < epe);
		break;

	  case 7: /* funcall */
		ef = (expr_f*)e;
		ap = ef->ap;
		for(ape = ef->sape; ap < ape; ap++)
			ed1opwalk(ap->e, f, varval);
		break;

	  default:
		fprintf(Stderr, "ed1opwalk bug! optype[%d] = %d\n",
			n, optypeb[n]);
		exit(1);
	  }
	}

 static void
#ifdef KR_headers
ed1oploop(c, n, f, varval) cde *c; int n; efunc **f, *varval;
#else
ed1oploop(cde *c, int n, efunc **f, efunc *varval)
#endif
{
	cde *ce = c + n;
	while(c < ce)
		ed1opwalk((c++)->e, f, varval);
	}

 void
#ifdef KR_headers
qp_opify_ASL(a) ASL *a;
#else
qp_opify_ASL(ASL *a)
#endif
{
	cexp *c, *ce;
	cexp1 *c1, *c1e;
	expr_v *v, *ve;
	ASL_fg *asl;
	efunc **f, *varval;
	int nv;

	ASL_CHECK(a, ASL_read_fg, "qp_opify");
	asl = (ASL_fg*)a;
	if (asl->i.rflags & ASL_opified)
		return;
	asl->i.rflags |= ASL_opified;
	f = r_ops_ASL;
	varval = r_ops_ASL[&f_OPVARVAL - &r_ops[0]];
	if ((nv = c_vars) < o_vars)
		nv = o_vars;
	nv += comb + comc + como + comc1 + como1;
	v = var_e;
	ve = v + nv;
	while(v < ve)
		(v++)->op = f_OPVARVAL;
	ed1oploop(obj_de, n_obj, f, varval);
	ed1oploop(con_de, asl->i.n_con0, f, varval);
	c1 = cexps1;
	c1e = c1 + comc1 + como1;
	while(c1 < c1e)
		ed1opwalk((c1++)->e, f, varval);
	c = cexps;
	ce = c + comb + comc + como;
	while(c < ce)
		ed1opwalk((c++)->e, f, varval);
	}
