/****************************************************************
Copyright (C) 1997-1998 Lucent Technologies
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

#ifndef NO_MAXCHECK
#include "nlp.h"
#include "r_opn.hd"
#ifndef f_MINLIST	/* TEMPORARY: now in updated r_opn.hd */
#define f_MINLIST	r_ops[11]
#define f_MAXLIST	r_ops[12]
#define f_ABS		r_ops[15]
#endif
#endif
#include "getstub.h"

extern void cfsqp ANSI((int nparam, int nf, int nfsr, int nineqn,
	int nineq, int neqn,
	int neq, int ncsrl, int ncsrn, int *mesh_pts,
	int mode, int iprint, int miter, int *inform, double bigbnd,
	double eps, double epseqn, double udelta, double *bl, double *bu,
	double *x, double *f, double *g, double *lambda,
	void (*obj)(int,  int,  double *,  double *, void *),
	void (*constr)(int, int, double *, double *, void *),
	void (*gradob)(int, int, double *, double *,
		void (*)(int, int, double *, double *, void *), void *),
	void (*gradcn)(int, int, double *, double *,
		void (*)(int, int, double *, double *, void *), void *),
	void *cd));

 extern char fsqp_version[];

 typedef void (*Dummy) ANSI((int, int, real*, real*, void *));

 int x_is_new = 1;
 real gLgeps = -1., objeps = -1., objrep = -1.;
 static real objsign;
 static int abs_mode, *cperm, nf, nx;
#ifdef NO_MAXCHECK
#define maxchk(b) 1
#else
 static int maxcheck = 1;
 static expr *obj0, **objp;
 de *mde;
 derp *mderp;

 static int
#ifdef KR_headers
maxchk(d) cde *d;
#else
maxchk(cde *d)
#endif
{
	de *c;
	expr *e;
	expr_va *v;
	int a, nf;

	nf = 1;
	if (!maxcheck)
		goto ret;
	e = d->e;
	if (e->op == (objsign > 0. ? f_MAXLIST : f_MINLIST)) {
		v = (expr_va *)e;
		mde = c = v->L.d;
		obj0 = e;
		objp = &d->e;
		mderp = v->R.D;
		nf = a = 0;
		while(e = (*c++).e) {
			nf++;
			if (e->op == f_ABS)
				a++;
			}
		if (a < nf || objsign < 0.)
			goto ret;
		abs_mode = 1;
		for(c = mde; e = c->e; c++) {
			c->e = e->L.e;
			e->dL = 1.;
			}
		}
	else if (e->op == f_ABS && objsign > 0.) {
		abs_mode = 1;
		d->e = e->L.e;
		e->dL = 1.;
		}
 ret:
	return nf;
	}
#endif /* NO_MAXCHECK */

 static void
#ifdef KR_headers
new_x(asl, x) ASL *asl; real *x;
#else
new_x(ASL *asl, real *x)
#endif
{
	x_is_new = 0;
	nx++;
	xknown(x);
	}

 static void
#ifdef KR_headers
obj(np, j, x, fj, cd) int np; int j; real *x; real *fj; char *cd;
#else
obj(int np, int j, real *x, real *fj, void *cd)
#endif
{
	ASL *asl = (ASL*)cd;
	int i = obj_no;
	Not_Used(np);
	if (x_is_new)
		new_x(asl, x);
#ifdef NO_MAXCHECK
	Not_Used(j);
	nf++;
#else
	if (!--j)
		nf++;
	if (objp)
		*objp = mde[j].e;
#endif
	*fj = objsign * objval(i, x, 0);
	}

 static void
#ifdef KR_headers
gradobj(np, j, x, g, d, cd) int np; int j; real *x; real *g; Dummy d; char *cd;
#else
gradobj(int np, int j, real *x, real *g, Dummy d, void *cd)
#endif
{
	ASL *asl = (ASL*)cd;
	int i = obj_no;
	Not_Used(d);

	if (x_is_new) {
		new_x(asl, x);
#ifndef NO_MAXCHECK
		if (objp) {
			*objp = obj0;
			objval(i, x, 0);
			}
#endif
		}
#ifdef NO_MAXCHECK
	Not_Used(j);
#else
	if (objp) {
		derp *D;
		de *d1 = mde + j - 1;
		*objp = obj0;
		if (D = mderp) {
			D->a.rp = d1->dv.rp;
			D->next = d1->d;
			}
		}
#endif
	objgrd(i, x, g, 0);
	if (objsign < 0.)
		for(i = 0; i < np; i++)
			g[i] = -g[i];
	}

 static void
#ifdef KR_headers
constr(np, j, x, c, cd) int np; int j; real *x; real *c; char *cd;
#else
constr(int np, int j, real *x, real *c, void *cd)
#endif
{
	ASL *asl = (ASL*)cd;
	Not_Used(np);

	if (x_is_new)
		new_x(asl, x);
	j = cperm[j-1];
	if (j < 0) {
		j = -(j+1);
		*c = LUrhs[j] - conival(j, x, 0);
		}
	else
		*c = conival(j, x, 0) - Urhsx[j];
	}

 static void
#ifdef KR_headers
gradcon(np, j, x, J, d, cd)
	int np; int j; real *x; real *J; Dummy d; char *cd;
#else
gradcon(int np, int j, real *x, real *J, Dummy d, void *cd)
#endif
{
	ASL *asl = (ASL*)cd;
	Not_Used(d);

	if (x_is_new)
		new_x(asl, x);
	j = cperm[j-1];
	if (j < 0) {
		j = -(j+1);
		congrd(j, x, J, 0);
		for(j = 0; j < np; j++)
			J[j] = -J[j];
		}
	else
		congrd(j, x, J, 0);
	}

static int maxiter = 200, objno = 1;
static int always, iprint, monotone;
static real eps = 1e-8, epseqn = 1e-8;

 static keyword
keywds[] = {	/* must be alphabetical */
 KW("always",	I_val,	&always,	"always check feasibility in linesearch  (0 = no) [mode: C]"),
 KW("eps",	D_val,	&eps,		"final Newton dir. norm (1e-8) [eps]"),
 KW("epseqn",	D_val,	&epseqn,	"max. nonlin. equality constraint violation (1e-8) [epseqn]"),
 KW("glgeps",	D_val,	&gLgeps,	"tol. for norm of Lagrangian (-1) [gLeps]"),
 KW("iprint",	I_val,	&iprint,	"0, 1, 2, 3, or 10*M+{2,3} (0) [iprint]"),
#ifndef NO_MAXCHECK
 KW("maxcheck",	I_val,	&maxcheck,	"check for minimax objective (1 = yes)"),
#endif
 KW("maxiter",	I_val,	&maxiter,	"maximum iterations (default 200) [miter]"),
 KW("monotone",	I_val,	&monotone,	"monotone iterates (0 = no) [mode: B]"),
 KW("objeps",	D_val,	&objeps,	"abs func change tol for probs with no nonlin eqns (-1) [objeps]"),
 KW("objno",	I_val,	&objno, 	"objective (1 = first)"),
 KW("objrep",	D_val,	&objrep,	"rel func change tol for probs with no nonlin eqns (-1) [objrep]"),
 KW("outlev",	I_val,	&iprint,	"synonym for iprint [iprint]"),
 KW("wantsol",	WS_val,	0,		WS_desc_ASL+5)
 };

 static char *usage_msg[] = {
  "\n	See Sections 4 and 5 of the CFSQP manual for details on",
  "	many of the  name=value  assignment possibilities; names",
  "	in square [brackets] when you say \"fsqp -=\" are used in",
  "	the CFSQP manual.",
  0};

 static Option_Info Oinfo = {
	"fsqp", "CFSQP 2.5", "fsqp_options", keywds, nkeywds, 1,
	fsqp_version, usage_msg };

 static void
#ifdef KR_headers
bcount(n, L, U, pnin, pneq) int n; real *L; real *U; int *pnin; int *pneq;
#else
bcount(int n, real *L, real *U, int *pnin, int *pneq)
#endif
{
	int neq = 0, nin = 0;

	while(n-- > 0) {
		if (*L == *U)
			neq++;
		else {
			if (*L > negInfinity)
				nin++;
			if (*U < Infinity)
				nin++;
			}
		L++, U++;
		}
	*pnin = nin;
	*pneq = neq;
	}

 static void
#ifdef KR_headers
bset(n, i, L, U, cpin, cpeq) int n; int i; real *L; real *U; int *cpin; int *cpeq;
#else
bset(int n, int i, real *L, real *U, int *cpin, int *cpeq)
#endif
{
	while(n-- > 0) {
		if (*L == *U)
			*cpeq++ = i;
		else {
			if (*L > negInfinity)
				*cpin++ = -(i+1);
			if (*U < Infinity)
				*cpin++ = i;
			}
		L++, U++, i++;
		}
	}

 static void
#ifdef KR_headers
report(asl, inform, m, nobj, y, f)
	ASL *asl; int inform, m, nobj; real *y; real *f;
#else
report(ASL *asl, int inform, int m, int nobj, real *y, real *f)
#endif
{
	int i, j;
	real F, *pi = pi0, t, t1;
	char msg[256];
	typedef struct { char *msg; int code, wantsol; } Sol_info;
	Sol_info *SI;
	static Sol_info solinfo[] = {
	 { "Optimal", 0, 1 },
	 { "Linear constraints are infeasible", 200, 0 },
	 { "Cannot satisfy the nonlinear constraints", 201, 0 },
	 { "Iteration limit", 400, 1 },
	 { "Line search failure", 501, 1 },
	 { "QP solver failed (d0)", 502, 1 },
	 { "QP solver failed (d1)", 503, 1 },
	 { "Input data are inconsistent (details with iprint > 0)", 202, 0 },
	 { "Overly tight stopping criteria (eps, epseqn)", 504, 1 },
	 { "Penalty parameter got too large", 505, 1 }
	 };

	t1 = -objsign;
	for(i = 0; i < m;) {
		t = y[i];
		j = cperm[i++];
		if (j < 0) {
			j = -(j+1);
			if (i < m && cperm[i] == j) {
				if (t > y[i])
					t = y[i];
				else
					t = -t;
				i++;
				}
			else
				t = -t;
			}
		pi[j] = t*t1;
		}
	if (inform >= 0 && inform < 10) {
		SI = solinfo + inform;
		solve_result_num = SI->code;
		i = sprintf(msg, "%s: %s\n", Oinfo.bsname, SI->msg);
		if (SI->wantsol)
		    if (nobj > 0) {
			F = f[0];
			j = 1;
			if (F < 0. && abs_mode)
				F = -F;
			while(j < nobj) {
				t = f[j++];
				if (t < 0. && abs_mode)
					t = -t;
				if (F < t)
					F = t;
				}
			i += sprintf(msg+i, "Objective %.*g\n", obj_prec(),
				objsign*F);
			}
		    else
			i += sprintf(msg+i, "No objective.\n");
		sprintf(msg+i,
			"%d trial points and %d objective evaluations.",
			nx, nf);
		}
	else {
		sprintf(msg, "%s: Unexpected value of inform = %d", inform);
		solve_result_num = 500;
		}
	write_sol(msg, X0, pi, &Oinfo);
	}

 int
#ifdef KR_headers
main(argc, argv) int argc; char **argv;
#else
main(int argc, char **argv)
#endif
{
	ASL *asl;
	char *stub;
	FILE *F;
	int inform, m, mode, n, neq, neqn, nineq, nineqn, nobj;
	real *f, *g, *y;
	static int mesh_pts[1];
	Not_Used(argc);

	asl = ASL_alloc(ASL_read_fg);
	stub = getstops(argv, &Oinfo);
	F = jac0dim(stub, (ftnlen)strlen(stub));
	LUv = (real*)Malloc(3*sizeof(real)*(n_con + n_var));
	Uvx = LUv + n_var;
	LUrhs = Uvx + n_var;
	Urhsx = LUrhs + n_con;
	X0 = Urhsx + n_con;
	pi0 = X0 + n_var;
	fg_read(F, 0);
	if (objno <= 0 || objno > n_obj)
		objno = 0;
	obj_no = objno - 1;
	nobj = 1;
	if (obj_no >= 0) {
		objsign = objtype[obj_no] ? -1. : 1.;
		nobj = maxchk(((ASL_fg*)asl)->I.obj_de_ + obj_no);
		}
	bcount(nlc, LUrhs, Urhsx, &nineqn, &neqn);
	bcount(n_con, LUrhs, Urhsx, &nineq, &neq);
	m = neq + nineq;
	n = m + n_var + nobj;
	g = (real*)Malloc((n+m+nobj)*sizeof(real) + m*sizeof(int));
	y = g + m;
	f = y + n;
	cperm = (int*)(f + nobj);
	bset(nlc, 0, LUrhs, Urhsx, cperm, cperm + nineq);
	bset(n_con-nlc, nineqn+neqn, LUrhs+nlc, Urhsx+nlc,
		cperm+nineqn, cperm + nineq + neqn);
	mode = 100;
	if (always)
		mode = 200;
	if (!monotone)
		mode += 10;
	mode |= abs_mode;
	if (objno <= 0)
		nobj = 0;
	cfsqp(n_var, nobj, 0, nineqn, nineq, neqn, neq, 0, 0,
		mesh_pts, mode, iprint, maxiter, &inform, Infinity,
		eps, epseqn, 0., LUv, Uvx, X0, f, g, y, obj, constr,
		gradobj, gradcon, asl);
	report(asl, inform, m, nobj, y + n_var, f);
	return 0;
	}
