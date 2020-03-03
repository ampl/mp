/*******************************************************************
Copyright (C) 2017 AMPL Optimization, Inc.; written by David M. Gay.

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

#include "asl_pfg.h"
#ifdef __cplusplus
extern "C" {
#endif

#undef funnelset
#define funnelset(x) funnelset_ASL((ASL_fg*)asl, x)

 int
xp1_check_ASL(ASL_pfg *asl, real *x)
{
	cexp *c, *c1, *ce;
	expr *e;
	expr_v *v, *v0;
	int *dvsp0, i0, i1, i, *vm;
	linarg *la;
	linpart *L, *Le;
	ograd *og;
	real t, *vscale, *xe;

	if (x0len == 0) {
		x0kind = 0;
		return 0;
		}
	if (x0kind == ASL_first_x)
		x0kind = 0;
	else if (!memcmp(Lastx, x, x0len))
		return 0;

	if (asl->i.Derrs)
		deriv_errclear_ASL(&asl->i);
	want_deriv = want_derivs;
	memcpy(Lastx, x, x0len);
	asl->i.nxval++;

	xe = (real*)((char*)x + x0len);
	v = v0 = var_e;
	vscale = asl->i.vscale;
	if ((vm = asl->i.vmap)) {
		if (vscale)
			while(x < xe)
				v0[*vm++].v = *vscale++ * *x++;
		else
			while(x < xe)
				v0[*vm++].v = *x++;
		}
	else {
		if (vscale)
			while(x < xe)
				(v++)->v = *vscale++ * *x++;
		else
			while(x < xe)
				(v++)->v = *x++;
		}

	for(la = asl->P.lalist; la; la = la->lnext) {
		og = la->nz;
		t = var_e[og->varno].v*og->coef;
		while((og = og->next))
			t += var_e[og->varno].v*og->coef;
		la->v->v = t;
		}
	if (asl->P.ncom) {
		c = cexps;
		dvsp0 = asl->P.dvsp0;
		i0 = *dvsp0;
		i = 0;
		/* Normally v == var_ex here, but v < var_ex is possible if sos_add */
		/* overestimated the number of new variables that sos_finish added. */
		v = var_ex;
		for(ce = c1 = c + asl->P.ncom; c < ce; c++) {
			for(i1 = *++dvsp0; i0 < i1; i0++, c1++) {
				cv_index = i0 + 1;
				e = c1->e;
				asl->P.vp[i0]->v = (*e->op)(e C_ASL);
				if (c1->funneled)
					funnelset(c1->funneled);
				}
			e = c->e;
			cv_index = ++i;
			t = (*e->op)(e C_ASL);
			if ((L = c->L))
				for(Le = L + c->nlin; L < Le; L++)
					t += L->fac * ((expr_v*)L->v.vp)->v;
			else if (!c->d && (og = asl->P.dv[c-cexps].ll)) {
				if (og->varno < 0) {
					t += og->coef;
					og = og->next;
					}
				while(og) {
					t += og->coef*v0[og->varno].v;
					og = og->next;
					}
				}
			(v++)->v = t;
			if (c->funneled)
				funnelset(c->funneled);
			}
		cv_index = 0;
		}
	return 1;
	}

 int
xp1known_ASL(ASL *asl, real *x, fint *nerror)
{
	Jmp_buf err_jmp0;
	int ij, rc;

	ASL_CHECK(asl, ASL_read_pfg, "xp1known");
	rc = 1;
	if (asl->i.xknown_ignore)
		goto ret;
	if (nerror && *nerror >= 0) {
		err_jmp = &err_jmp0;
		ij = setjmp(err_jmp0.jb);
		if ((*nerror = ij))
			goto done;
		}
	errno = 0;	/* in case f77 set errno opening files */
	co_index = nlo ? -1 : 0;
	rc = xp1_check_ASL((ASL_pfg*)asl, x);
	asl->i.x_known = 1;
 done:
	err_jmp = 0;
 ret:
	return rc;
	}

#ifdef __cplusplus
	}
#endif
