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

/* For replacing objectives of the form v or c*v, where c is a constant	  */
/* and variable v has no bounds and appears linearly in one constraint    */
/* (of the right sense if an inequality) with the body of the constraint, */
/* then removing v and the constraint from the problem seen by the solver.*/

#include "nlp.h"
#define SKIP_NL2_DEFINES
#undef f_OPNUM
#include "psinfo.h"
#include "jacpdim.h"
#undef ps_func
#undef psb_elem
#include "r_qp.hd" /* for OPNUM */

 struct
Objrep {
	/* For minimization (or maximization, with suitable inequality */
	/* sense adjustments) of Objective = c0 + c1*v */
	/* with v = _var[ivo] defined by constraint ico of the form */
	/* (_con[ico].body = c2*v + f(x)) >= rhs  with c2 > 0   or */
	/* (_con[ico].body = c2*v + f(x)) <= rhs  with c2 < 0, */
	/* We change the Objective to */
	/* Objective = c0 + c1*(rhs - f(x))/c2 = c0a + c12*f(x) */
	/* with c12 = -c1 / c2  and  c0a = c0 - c12*rhs. */
	/* Objective gradient = c12 * grad(f(x)). */
	/* We remove constraint _con[ico] and variable v; writesol */
	/* calls obj_adj_xy_ASL to compute v = (Objectve - c0)/c1 */
	/* and dual variable value -c12 for the removed constraint. */

	int ico;	/* index of constraint that gives the objective */
	int ivo;	/* index of objective variable */
	int nxval;	/* X generation for f */
	real c0, c0a, c1, c12, f;
	};

 static void
obj_adj1(ASL *asl, int no)
{
	Objrep *od, **pod;
	cgrad **Cgrd, **Cgrd0, *cg, *cgo, **pcg, **pcge;
	efunc_n *op;
	expr_n *e;
	int co, cv, flags, i, incc, incv, j, k, m, n;
	int *Cvar, *cm, *vm, *zg, **zgp;
	ograd *og;
        ps_func *P;
        ps_func2 *P2;
	real *Lc, *Lv, *Uc, *Uv, c1, c2, rhs, t;

	op = f_OPNUM_ASL;
	switch (asl->i.ASLtype) {
	 case ASL_read_fg:
		e = (expr_n*)((ASL_fg*)asl)->I.obj_de_[no].e;
		break;
	 case ASL_read_fgh:
		e = (expr_n*)((ASL_fgh*)asl)->I.obj2_de_[no].e;
		break;
         case ASL_read_pfg:
		e = (expr_n*)((ASL_pfg*)asl)->I.obj_de_[no].e;
		P = &((ASL_pfg*)asl)->P.ops[no];
		if (P->nb || P->ng)
			return;
		op = (efunc_n*)OPNUM;
		break;
	 case ASL_read_pfgh:
		e = (expr_n*)((ASL_pfgh*)asl)->I.obj2_de_[no].e;
		P2 = &((ASL_pfgh*)asl)->P.ops[no];
		if (P2->nb || P2->ng)
			return;
		op = (efunc_n*)OPNUM;
		break;
	 default:
		fprintf(Stderr, "Bug: surprise ASLtype = %d in obj_adj\n", asl->i.ASLtype);
		exit(1);
	 }
	if (e->op != op)
		return;
	og = Ograd[no];
	if (!og || og->next)
		return;
	cv = og->varno;
	if (cv < nlvc)
		return;
	if (!(c1 = og->coef))
		return;
	if (Uvx) {
		Lv = LUv + cv;
		Uv = Uvx + cv;
		incv = 1;
		}
	else {
		Lv = LUv + 2*cv;
		Uv = Lv + 1;
		incv = 2;
		}
	if (*Lv > negInfinity || *Uv < Infinity)
		return;
	pcg = Cgrd = Cgrad;
	cgo = 0;
	for(pcge = pcg + n_con; pcg < pcge; ++pcg) {
		for(cg = *pcg; cg; cg = cg->next)
			if (cg->varno == cv) {
				if (cgo)
					return;
				cgo = cg;
				co = pcg - Cgrd;
				for(k = 0, cg = *pcg; cg; cg = cg->next)
					++k;
				break;
				}
		}
	if (!cgo)
		return;
	if (n_cc && cvar[co])
		return;
	if ((c2 = cgo->coef) == 0.)
		return;
	t = c1 / c2;
	j = t < 0.;
	if (objtype[no])
		j = 1 - j;
	if (Urhsx) {
		Lc = LUrhs + co;
		Uc = Urhsx + co;
		incc = 1;
		}
	else {
		Lc = LUrhs + 2*co;
		Uc = Lc + 1;
		incc = 2;
		}
	if (j) {
		if ((rhs = *Uc) >= Infinity)
			return;
		}
	else {
		if ((rhs = *Lc) <= negInfinity)
			return;
		}
	flags = asl->i.rflags;
	if (*Lc < *Uc) {
		if (!(flags & ASL_obj_replace_ineq))
			return;
		}
	else {
		if (!(flags & ASL_obj_replace_eq))
			return;
		}

	if (co < nlc) {
		--nlc;
		++nlo;
		}
	nzc -= k;
	nzo += k - 1;
	pod = asl->i.Or;
	od = (Objrep*)M1alloc(sizeof(Objrep)
		+ (pod ? 0 : n_obj*sizeof(Objrep*)));
	if (!pod) {
		pod = asl->i.Or = (Objrep**)(od+1);
		for(k = n_obj; --k >= 0; )
			pod[k] = 0;
		}
	pod[no] = od;
	od->ico = co;
	od->ivo = cv;
	od->c0 = e->v;
	od->c0a = e->v + t*rhs;
	od->c1 = c1;
	od->c12 = -t;
	od->nxval = -1;
	od->f = 0.;

	pcg = &Cgrd[co];
	while((cg = *pcg) != cgo)
		pcg = &cg->next;
	*pcg = cgo->next;
	if ((cm = asl->i.cmap) && (Cgrd0 = asl->i.Cgrad0))
		Cgrd0[cm[co]] = Cgrd[co];

	m = --n_con;
	if (n_conjac[1] > m)
		n_conjac[1] = m;
	if (co != m) {
		cm = get_vcmap_ASL(asl, ASL_Sufkind_con);
		pcg = Cgrd;
		Cvar = cvar;
		for(i = co; i < m; i = j) {
			cm[i] = cm[j = i + 1];
			pcg[i] = pcg[j];
			if (Cvar)
				Cvar[i] = Cvar[j];
			}
		cm[m] = -1;
		}
	for(i = co; i < m; ++i) {
		*Lc = Lc[incc];
		*Uc = Uc[incc];
		Lc += incc;
		Uc += incc;
		}
	n = --n_var;
	if (cv != n) {
		vm = get_vcmap_ASL(asl, ASL_Sufkind_var);
		for(i = cv; i < n; ++i)
			vm[i] = vm[i+1];
		vm[n] = -1;
		}
	for(i = cv; i < n; ++i) {
		*Lv = Lv[incv];
		*Uv = Uv[incv];
		Lv += incv;
		Uv += incv;
		}
	if ((zgp = zerograds))
		for(zg = zgp[no]; *zg >= 0; ++zg)
			if (*zg >= cv)
				--*zg;
	}

 static real
objval_adj(ASL *asl, int no, real *X, fint *nerror)
{
	Objrep *od;
	real c;

	if (!(od = asl->i.Or[no]))
		return asl->p.Objval_nomap(asl, no, X, nerror);
	c = asl->p.Conival_nomap(asl, od->ico, X, nerror);
	if (nerror && *nerror)
		return 0.;
	od->nxval = asl->i.nxval;
	return od->f = od->c0a + od->c12*c;
	}

 static void
objgrd_adj(ASL *asl, int no, real *X, real *G, fint *nerror)
{
	Objrep *od;
	cgrad *gr;
	int k, *vmi;
	real c;

	if (!(od = asl->i.Or[no])) {
		asl->p.Objgrd_nomap(asl, no, X, G, nerror);
		return;
		}
	if (od->nxval != asl->i.nxval)
		objval_adj(asl, no, X, nerror);
	if ((k = asl->i.congrd_mode))
		asl->i.congrd_mode = 0;
	asl->p.Congrd_nomap(asl, od->ico, X, G, nerror);
	asl->i.congrd_mode = k;
	if ((c = od->c12) != 1. && (!nerror || !*nerror)) {
		vmi = get_vminv_ASL(asl);
		for(gr = asl->i.Cgrad0[od->ico]; gr; gr = gr->next)
			G[vmi[gr->varno]] *= c;
		}
	}

 static void
jac_adj(ASL *asl)
{
	cgrad **Cgrd, *cg;
	int i, k, k0, nc, nv, *x;

	nc = n_con;
	Cgrd = Cgrad;
	if (asl->i.rflags & ASL_rowwise_jac) {
		k = 0;
		for(i = 0; i < nc; ++i)
			for(cg = Cgrd[i]; cg; cg = cg->next)
				cg->goff = k++;
		return;
		}
	nv = asl->i.n_var0;
	x = (int*)Malloc(nv*sizeof(int));
	memset(x, 0, nv*sizeof(int));
	for(i = 0; i < nc; ++i)
		for(cg= Cgrd[i]; cg; cg = cg->next)
			++x[cg->varno];
	for(i = k = 0; i < nv; ++i) {
		k0 = k;
		k += x[i];
		x[i] = k0;
		}
	for(i = 0; i < nc; ++i)
		for(cg = Cgrd[i]; cg; cg = cg->next)
			cg->goff = x[cg->varno]++;
	free(x);
	}

 static void
paradj(ASL *asl, int *pno, real **pow, real **py)
{
	Objrep *Od, **Or;
	int *cm, i, j, k, nc, nc0, needow, no, nobj;
	real *ow, wo, *ws, *y, *ys;

	nobj = n_obj;
	no = *pno;
	Od = 0;
	ow = *pow;
	wo = 1.;
	Or = asl->i.Or;
	ws = asl->i.orscratch;
	ys = ws + nobj;
	memset(ys, 0, asl->i.n_con0*sizeof(real));
	cm = asl->i.cmap;
	if ((y = *py)) {
		*py = ys;
		nc = n_con;
		nc0 = asl->i.n_con0;
		for(i = 0; i < nc; ++i)
			if ((j = cm[i]) < nc0)
				ys[j] = y[i];
		}
	needow = 1;
	if (no >= 0 && no < nobj) {
		k = nobj;
		if (ow && (wo = ow[no]) == 0.)
			no = *pno = -1;
		else
			Od = Or[no];
		ow = 0;
		}
	else if (ow) {
		needow = 0;
		for(i = k = 0; i < nobj; ++i) {
			if (ow[i]) {
				if (Or[i]) {
					if (!Od) {
						k = i;
						Od = Or[i];
						wo = ow[i];
						continue;
						}
					}
				else
					++needow;
				}
			}
		if (!needow)
			ow = 0;
		}
	if (Od) {
		*pno = -1;
		if (!y)
			*py = ys;
		if (ow) {
			memcpy(*pow = ws, ow, nobj*sizeof(real));
			ow = ws;
			wo = ow[k];
			}
 loop:
		ys[Od->ico] = wo * Od->c12;
		if (ow) {
			ow[k] = 0.;
			while(++k < nobj)
				if ((wo = ow[k]) && (Od = Or[k]))
					goto loop;
			}
		}
	*pow = ow;
	}

 static void
hvcomp_adj(ASL *asl, real *hv, real *p, int no, real *ow, real *y)
{
	paradj(asl, &no, &ow, &y);
	asl->p.Hvcomp_nomap(asl, hv, p, no, ow, y);
	}

 static void
hvinit_adj(ASL *asl, int hid_limit, int no, real *ow, real *y)
{
	paradj(asl, &no, &ow, &y);
	asl->p.Hvinit_nomap(asl, hid_limit, no, ow, y);
	}

 static void
duthes_adj(ASL *asl, real *H, int no, real *ow, real *y)
{
	paradj(asl, &no, &ow, &y);
	asl->p.Duthes_nomap(asl, H, no, ow, y);
	}

 static void
fulhes_adj(ASL *asl, real *H, fint LH, int no, real *ow, real *y)
{
	paradj(asl, &no, &ow, &y);
	asl->p.Fulhes_nomap(asl, H, LH, no, ow, y);
	}

 static void
sphes_adj(ASL *asl, SputInfo **spi, real *H, int no, real *ow, real *y)
{
	paradj(asl, &no, &ow, &y);
	asl->p.Sphes_nomap(asl, spi, H, no, ow, y);
	}

 static fint
sphes_setup_adj(ASL *asl, SputInfo **spi, int no, int ow, int y, int uptri)
{
	Objrep *od, **odp;

	if (no >= 0 && no < n_obj && (odp = asl->i.Or) && (od = odp[no])) {
		no = -1;
		ow = 0;
		y = 1;
		}
	return asl->p.Sphset_nomap(asl, spi, no, ow, y, uptri);
	}

 void
obj_adj_ASL(ASL *asl)
{
	int nc0, no, nobj;

	if (A_vals)
		return; /* for now, no adjustment when A_vals is used */
	nobj = n_obj;
	nc0 = n_con;
	for(no = 0; no < nobj; ++no)
		obj_adj1(asl, no);
	if (asl->i.Or) {
		if (asl->i.ASLtype != ASL_read_f) {
			asl->p.Objval = objval_adj;
			asl->p.Objgrd = objgrd_adj;
			if (asl->i.ASLtype == ASL_read_pfgh) {
				asl->p.Hvcomp = hvcomp_adj;
				asl->p.Hvinit = hvinit_adj;
				asl->p.Duthes = duthes_adj;
				asl->p.Fulhes = fulhes_adj;
				asl->p.Sphes  = sphes_adj;
				asl->p.Sphset = sphes_setup_adj;
				}
			}
		jac_adj(asl);
		asl->i.orscratch = (real*)M1zapalloc((nc0 + nobj)*sizeof(real));
		}
	}

 void
obj_adj_xy_ASL(ASL *asl, real *x, real *x0, real *y)
{
	Objrep *od, **odp;
	fint nerror;
	int no, nobj;

	odp = asl->i.Or;
	nobj = n_obj;
	for(no = 0; no < nobj; ++no)
		if ((od = odp[no])) {
			if (od->nxval != asl->i.nxval) {
				nerror = 0;
				objval_adj(asl, no, x0, &nerror);
				if (nerror)
					continue;
				}
			x[od->ivo] = (od->f - od->c0) / od->c1;
			if (y)
				y[od->ico] = -od->c12;
			}
	}
