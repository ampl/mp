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
#define SKIP_NL2_DEFINES
#undef f_OPNUM
#include "psinfo.h"
#include "jacpdim.h"
#undef cde
#undef ps_func

 struct
MPEC_Adjust {
	int *cc;		/* indices of constraints originally involved in complementarities */
	int *cce;		/* end of cc array */
	int *ck;		/* kind of complementarity modification */
	real *rhs1;		/* start of added constraint right-hand sides */
	cgrad **Cgrda;		/* linear terms added to existing constraints */
	int incc;		/* for incrementing constraint lhs and rhs */
	int incv;		/* for incrementing variable   lhs and rhs */
	int m0;			/* original number of constraints */
	int n0;			/* original number of variables */
	};

 static expr_n ZeroExpr = { f_OPNUM_ASL, 0. };

 static void
reverse(int *a, int *b)
{
	int t;

	while(--b > a) {
		t = *b;
		*b = *a;
		*a++ = t;
		}
	}

 void
mpec_adjust_ASL(ASL *asl)
{
	MPEC_Adjust *mpa;
	cde *cd;
	cde2 *cd2;
	cgrad **Cgrd, **Cgrd1, **Cgrda, *cg, *cg1, *ncg, **pcg;
	char *hx0;
	int *cc, *ck, *cv, *ind1, *ind2, *map, *mapinv;
	int i, incc, incv, j, k, m, m0, n, n0, n1, nb, ncc, ncc0, nib, nib0;
	int nnv, nz, nz0, nznew, v1, v2, v3, v4;
	real *Lc, *Lc0, *Lc1, *Lv, *Lv0, *Lv1, *Uc, *Uc0, *Uc1, *Uv, *Uv0, *Uv1;
	real a, b, *x;
	extern void f_OPVARVAL_ASL(), f2_VARVAL_ASL();

	n = n1 = n0 = n_var;
	nib = niv + nbv;
	nib0 = n - nib;	/* offset of first linear integer or binary variable */
	m = m0 = n_con;
	nz = nz0 = nzc;
	cv = cvar;
	Cgrd = Cgrad;
	Cgrd1 = Cgrd + m;
	incc = incv = 1;
	Lc0 = LUrhs;
	if (!(Uc0 = Urhsx)) {
		Uc0 = Lc0 + 1;
		incc = 2;
		}
	Lv0 = LUv;
	if (!(Uv0 = Uvx)) {
		Uv0 = Lv0 + 1;
		incv = 2;
		}
	ncc = ncc0 = n_cc;
	Lc1 = Lc0 + m*incc;
	Uc1 = Uc0 + m*incc;
	Lv1 = Lv0 + n*incv;
	Uv1 = Uv0 + n*incv;

	for(i = k = 0; i < m0; ++i)
		if ((j = cv[i])) {
			++k;
			Lc = Lc0 + incc*i;
			Uc = Uc0 + incc*i;
			nb = (*Lc > negInfinity) + (*Uc < Infinity);
			/* nb == 0 or 1 */
			if (!nb) {
				m += 2;
				n += 4;
				nz += 6;
				++ncc;
				}
			else {
				Lv = Lv0 + incv*--j;
				if (*Lv != 0.) {
					++m;
					++n;
					nz += 2;
					}
				/* Even if constraint i has the form v >= 0, */
				/* add a new variable v1 >= 0 and change the */
				/* constraint to v1 = v - rhs, in case v is  */
				/* involved in more than one complementarity */
				++n;
				++nz;
				}
			}
	if (k != ncc0) {
		fprintf(Stderr,
			"\nERROR: mpec_adjust saw %d rather than %d incoming complementarities.\n",
			k, ncc0);
		exit(1);
		}
	n_var = n;
	n_con = m;
	nnv = n - n0;
	if (n_obj)
		adjust_zerograds_ASL(asl, nnv);
	if (n_conjac[1] >= m0)
		n_conjac[1] = m;
	nzc = nz;
	n_cc = ncc;
	nznew = nz - nz0;
	ncg = (cgrad*)M1alloc(2*(ncc + ncc0)*sizeof(int) + nznew*sizeof(cgrad)
			+ ncc0*sizeof(cgrad*) + sizeof(MPEC_Adjust));
	asl->i.mpa = mpa = (MPEC_Adjust*)(ncg + nznew);
	Cgrda = mpa->Cgrda = (cgrad**)(mpa + 1);
	asl->i.ccind1 = ind1 = (int*)(Cgrda + ncc0);
	asl->i.ccind2 = ind2 = ind1 + ncc;
	mpa->cc = cc = ind2 + ncc;
	mpa->ck = ck = mpa->cce = cc + ncc0;
	mpa->m0 = m0;
	mpa->n0 = n0 - nib;
	mpa->rhs1 = Lc1;
	mpa->incc = incc;
	mpa->incv = incv;
	if (nib) {
		map = get_vcmap_ASL(asl, ASL_Sufkind_var);
		/* Three reverse calls move nib values of map up nnv places. */
		j = n0 - nib;
		reverse(map+j, map + n0 + nnv);
		reverse(map+j, map + j + nnv);
		reverse(map + j + nnv, map + n0 + nnv);
		i = n0 + nnv;
		while(--i >= n0) {
			j = i - nnv;
			Lv0[incv*i] = Lv0[incv*j];
			Uv0[incv*i] = Uv0[incv*j];
			}
		if ((x = X0)) {
			i = n0 + nnv;
			while(--i >= n0)
				x[i] = x[i-nnv];
			for(i = n0 - nnv; i < n0; ++i)
				x[i] = 0.;
			if ((hx0 = havex0)) {
				for(i = n0 + nnv; --i >= n0; )
					hx0[i] = hx0[i-nnv];
				for(i = n0 - nnv;i < n0; ++i)
					hx0[i] = 0;
				}
			}
		Lv1 -= j = incv*nib;
		Uv1 -= j;
		}
	else {
		if ((map = asl->i.vmap)) {
			j = asl->i.n_var0;
			for(i = n0; i < n; ++i)
				map[i] = -1;
			}
		if ((x = X0)) {
			memset(x + n0, 0, nnv*sizeof(real));
			if ((hx0 = havex0))
				memset(hx0 + n0, 0, nnv);
			}
		}
#define vset(x,y) *x = y; x += incv;
	for(i = 0; i < m0; ++i)
		if ((j = cv[i])) {
			if (j > nib0)
				j += nnv;
			*cc++ = i;
			pcg = &Cgrd[i];
			cg = 0;
			while((cg1 = *pcg))
				pcg = &(cg = cg1)->next;
			*Cgrda++ = cg;
			Lc = Lc0 + incv*i;
			Uc = Uc0 + incc*i;
			Lv = Lv0 + incv*--j;
			Uv = Uv0 + incv*j;
			a = *Lc;
			b = *Uc;
			*ck++ = nb = (a > negInfinity) + (b < Infinity);
			if (nb == 0) {
				/* change L <= v = _svar[j] <= U */
				/* and -Infinity <= body <=  Infinity into */
				/* v1 = v - L >= 0, v2 = U - v >= 0, */
				/* v3 - v4 = body, v3 >= 0, v4 >= 0, */
				/* v1 complements v3, v2 complements v4 */

				*Lc = *Uc = 0.;
				v1 = n1++;
				v2 = n1++;
				v3 = n1++;
				v4 = n1++;
				for(k = 0; k < 4; ++k) {
					vset(Lv1, 0.);
					vset(Uv1, Infinity);
					}
				ncg[1].varno = v4;
				ncg[1].coef = 1.;
				ncg[1].next = 0;
				ncg[0].varno = v3;
				ncg[0].coef = -1.;
				ncg[0].next = &ncg[1];
				*pcg = ncg;
				ncg += 2;
				ncg[1].varno = v1;
				ncg[1].coef = -1.;
				ncg[1].next = 0;
				ncg[0].varno = j;
				ncg[0].coef = 1.;
				ncg[0].next = &ncg[1];
				*Lc1 = *Uc1 = *Lv;
				Lc1 += incc;
				Uc1 += incc;
				*Cgrd1++ = ncg;
				ncg += 2;
				ncg[1].varno = v2;
				ncg[1].coef = 1.;
				ncg[1].next = 0;
				ncg[0].varno = j;
				ncg[0].coef = 1.;
				ncg[0].next = &ncg[1];
				*Lc1 = *Uc1 = *Uv;
				Lc1 += incc;
				Uc1 += incc;
				*Cgrd1++ = ncg;
				ncg += 2;
				*ind1++ = v1;
				*ind2++ = v3;
				*ind1++ = v2;
				*ind2++ = v4;
				}
			else {
				/*nb == 1*/
				v1 = j;
				if (*Lv != 0.) {
					/* For v = _svar[j], replace */
					/* v >= a with v1 = v - a, v1 >= 0, or */
					/* v <= b with v1 = b - v, v1 >= 0 */
					v1 = n1++;
					vset(Lv1, 0.);
					vset(Uv1, Infinity);
					ncg[1].varno = v1;
					ncg[1].next = 0;
					ncg[0].varno = j;
					ncg[0].coef = 1.;
					ncg[0].next = &ncg[1];
					if (*Lv > negInfinity) {
						ncg[1].coef = -1.;
						*Lc1 = *Uc1 = *Lv;
						}
					else {
						ncg[1].coef = 1.;
						*Lc1 = *Uc1 = *Uv;
						}
					Lc1 += incc;
					Uc1 += incc;
					*Cgrd1++ = ncg;
					ncg += 2;
					}
				ncg->varno = v2 = n1++;
				ncg->next = 0;
				vset(Lv1, 0.);
				vset(Uv1, Infinity);
				if (*Lv > negInfinity) {
					ncg->coef = -1.;
					*Uc = *Lc;
					}
				else {
					ncg->coef = 1.;
					*Lc = *Uc;
					}
				*pcg = ncg++;
				*ind1++ = v1;
				*ind2++ = v2;
				}
			}
#undef vset
	if (map) {
		ind1 -= ncc;
		ind2 -= ncc;
		mapinv = get_vminv_ASL(asl);
		for(i = 0; i < ncc; ++i) {
			ind1[i] = mapinv[ind1[i]];
			ind2[i] = mapinv[ind2[i]];
			}
		}
	if ((map = asl->i.cmap)) {
		j = asl->i.n_con0;
		Cgrd1 = asl->i.Cgrad0;
		for(i = m0; i < m; ++i) {
			map[i] = -1;
			Cgrd1[j++] = Cgrd[i];
			}
		}
	i = m0;
	k = m - m0;
	switch(asl->i.ASLtype) {
	  case ASL_read_pfg:
		memset(((ASL_pfg*)asl)->P.cps + m0, 0, k*sizeof(ps_func));
		cd = ((ASL_pfg*)asl)->I.con_de_;
		goto have_cd;
	  case ASL_read_f:
	  case ASL_read_fg:
		cd = ((ASL_fg*)asl)->I.con_de_;
 have_cd:
		while(i < m)
			cd[i++].e = (expr*)&ZeroExpr;
		break;
	  case ASL_read_fgh:
		cd2 = ((ASL_fgh*)asl)->I.con2_de_;
		goto have_cd2;
	  case ASL_read_pfgh:
		memset(((ASL_pfgh*)asl)->P.cps + m0, 0, k*sizeof(ps_func2));
		cd2 = ((ASL_pfgh*)asl)->I.con2_de_;
 have_cd2:
		while(i < m)
			cd2[i++].e = (expr2*)&ZeroExpr;
	  }

	}

 void
mpec_auxvars_ASL(ASL *asl, real *c, real *x)
{
	/* Adjust variables added by mpec_adjust_ASL() so the constraints */
	/* added by mpec_adjust_ASL() are satisfied. */

	MPEC_Adjust *mpa;
	cgrad **Cg, **Cga, *cg;
	int *cc, *cce, *ck, *cv, i, incc, incv, j, m0, n0;
	real *Lc, *Lc0, *Lc1, *Lv0, *ca, t;

	mpa = asl->i.mpa;
	cv = cvar;
	cc = mpa->cc;
	cce = mpa->cce;
	ck = mpa->ck;
	Cga = mpa->Cgrda;
	m0 = mpa->m0;
	n0 = mpa->n0;
	Cg = Cgrad + m0;
	ca = c + m0;
	Lc0 = LUrhs;
	Lc1 = mpa->rhs1;
	Lv0 = LUv;
	incc = mpa->incc;
	incv = mpa->incv;
	do {
		t = c[i = *cc++];
		c[i] = 0.;
		j = cv[i] - 1;
		for(cg = *Cga++; cg->varno < n0; cg = cg->next);
		Lc = Lc0 + i*incc;
		if (!*ck++) {
			if (t >= 0.)
				x[cg->varno] = t;
			else {
				cg = cg->next;
				x[cg->varno] = -t;
				}
			cg = (*Cg++)->next;
			x[cg->varno] = x[j] - *Lc1;
			*ca++ = *Lc1;
			Lc1 += incc;
			cg = (*Cg++)->next;
			x[cg->varno] = *Lc1 - x[j];
			*ca++ = *Lc1;
			Lc1 += incc;
			}
		else {
			x[cg->varno] = cg->coef*(*Lc - t);
			c[i] = *Lc;
			if (Lv0[incv*j] != 0.) {
				cg = (*Cg++)->next;
				x[cg->varno] = cg->coef*(*Lc1 - x[j]);
				*ca++ = *Lc1;
				Lc1 += incc;
				}
			}
		} while(cc < cce);
	}
