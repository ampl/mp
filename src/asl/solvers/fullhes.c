/****************************************************************
Copyright (C) 1997, 2001 Lucent Technologies
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
#ifdef KR_headers
add_op(H, og0, t, n) real *H; ograd *og0; real t; fint n;
#else
add_op(real *H, ograd *og0, real t, fint n)
#endif
{
	ograd *og, *og1;
	real *Hj, t1;
	int j;

	for(og = og0; og; og = og->next) {
		if (t1 = t * og->coef) {
			j = og->varno;
			Hj = H + n*j;
			for(og1 = og0;; og1 = og1->next) {
				Hj[og1->varno] += t1*og1->coef;
				if (og == og1)
					break;
				}
			}
		}
	}

 void
#ifdef KR_headers
fullhes_ASL(a, H, LH, nobj, ow, y)
	ASL *a; real *H, *ow, *y; fint LH; int nobj;
#else
fullhes_ASL(ASL *a, real *H, fint LH, int nobj, real *ow, real *y)
#endif
{
	/* full Hessian */
	int i, j, n, no, noe;
	linarg *la, **lap, **lap1, **lape;
	expr_v *v;
	range *r, *r0;
	real *Hi, *Hj, *Hje, *cscale, *owi, *s, *si, t, t1;
	ograd *og, *og1;
	ps_func *p, *pe;
	psg_elem *g, *ge;
	ASL_pfgh *asl;

	asl = pscheck_ASL(a, "fullhes");
	xpsg_check_ASL(asl, nobj, ow, y);
	if (nobj >= 0 && nobj < n_obj) {
		no = nobj;
		noe = no + 1;
		owi = ow ? ow + no : &edag_one_ASL;
		}
	else {
		nobj = -1;
		no = noe = 0;
		if (owi = ow)
			noe = n_obj;
		}

	if (!asl->P.hes_setup_called)
		(*asl->p.Hesset)(a, 1, 0, nlo, 0, nlc);
	s = asl->P.dOscratch;
	n = c_vars >= o_vars ? c_vars : o_vars;
	if (n <= 0)
		return;

	/* First compute upper triangle */

	for(Hj = H, i = 1; i <= n; Hj += LH - i++)
		for(Hje = Hj + i; Hj < Hje; )
			*Hj++ = 0;

	r0 = (range*)&asl->P.rlist;
	for(r = asl->P.rlist.next; r != r0; r = r->rlist.next) {
		if ((j = r->n) <= 0)
			continue;
		lap = r->lap;
		lape = lap + j;
		si = s;
		while(lap < lape) {
			*si = 1;
			pshv_prod_ASL(asl, r, nobj, ow, y);
			*si++ = 0;
			la = *lap++;
			for(og = la->nz; og; og = og->next) {
				t = og->coef;
				j = og->varno;
				Hj = H + LH*j;
				for(lap1 = r->lap; lap1 < lape; ) {
					la = *lap1++;
					v = la->v;
					if (!(t1 = t * v->aO))
						continue;
					for(og1 = la->nz;
					    og1 && (i = og1->varno) <= j;
					    og1 = og1->next)
						Hj[i] += t1*og1->coef;
					}
				}
			}
		}
	if (asl->P.nobjgroups)
		for(; no < noe; no++)
			if (t = *owi++) {
				p = asl->P.ops + no;
				g = p->g;
				for(ge = g + p->ng; g < ge; g++)
					if (g->g2)
						add_op(H, g->og, t*g->g2, LH);
				}
	if (asl->P.ncongroups && y) {
		cscale = asl->i.lscale;
		p = asl->P.cps;
		for(pe = p + n_con; p < pe; p++, y++)
			if (t = cscale ? *cscale++ * *y : *y)
				for(g = p->g, ge = g + p->ng; g < ge; g++)
					if (g->g2)
						add_op(H, g->og, t*g->g2, LH);
		}

	if (s = asl->i.vscale) {
		Hi = H;
		for(i = 0; i < n; i++, Hi += LH) {
			t = s[i];
			for(j = 0; j <= i; j++)
				Hi[j] *= t * s[j];
			}
		}

	/* Now copy upper triangle to lower triangle. */

	for(Hj = H, i = 1; i < n; i++) {
		Hi = H + i;
		Hj = H + i*LH;
		Hje = Hj + i;
		while(Hj < Hje) {
			*Hi = *Hj++;
			Hi += LH;
			}
		}
	}
