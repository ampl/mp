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

#include "getstub.h"
#define asl cur_ASL

real macheps =
#ifdef IEEE_MC68k
	2.220446049250313e-16
#else
#ifdef IEEE_8087
	2.220446049250313e-16
#else
#ifdef VAX
	2.77555756156289135e-17
#else
	MACHEPS
#endif
#endif
#endif
	;

typedef /* Subroutine */ int (*Sfun)(fint *, real *, real *, real *);

Cextern int lmqn_(fint *ifail, fint *n, real *x, real *f, real *g,
		real *w, fint *lw, Sfun,
		fint *msglvl, fint *maxit, fint *maxfun,
		real *eta, real *stepmx, real *accrcy, real *xtol);

Cextern int lmqnbc_(fint *ifail, fint *n, real *x, real *f, real *g,
		real *w, fint *lw, Sfun,
		real *low, real *up, fint *ipivot,
		fint *msglvl, fint *maxit, fint *maxfun,
		real *eta, real *stepmx, real *accrcy, real *xtol);

fint maxfun, maxit, msglvl;
int nprob, solprt, x0prt;
real accrcy, eta, stepmx, xtol;

static int nfcall;

static keyword keywds[] = {	/* must be in alphabetical order */
 KW("accrcy", D_val, &accrcy, "accuracy of computed function values"),
 KW("eta", D_val, &eta, "severity of the linesearch"),
 KW("maxfun", L_val, &maxfun, "maximum function evaluations"),
 KW("maxit", L_val, &maxit, "maximum number of inner iterations per step"),
 KW("msglvl", L_val, &msglvl, "printing:  0 = none, 1 = one line per major iteration"),
 KW("nprob", I_val, &nprob, "objective choice: 0 (default) = 1st"),
 KW("solprt", I_val, &solprt, "print final solution"),
 KW("stepmx", D_val, &stepmx, "maximum step in the linesearch"),
 KW("wantsol", WS_val, 0, WSu_desc_ASL+5),
 KW("x0prt", L_val, &x0prt, "print initial x"),
 KW("xtol", D_val, &xtol, "desired accuracy for the solution x*")
 };

 static Option_Info Oinfo =
	{ "tn", "TN", "tn_options", keywds, nkeywds, 1 };

 int
sfun(fint *N, real *X, real *F, real *G)
{
	fint nerror = 0;

	*F = objval(nprob, X, &nerror);
	if (!nerror)
		objgrd(nprob, X, G, &nerror);
	nfcall++;
	return nerror;
	}

 static void
show_x(char *when, fint n, int nb, real *X, real *L, real *U)
{
	fint i;

	if (nb) {
		printf("i\t%s x(i)\tL(i)\t\tU(i)\n", when);
		for(i = 1; i <= n; i++)
			printf("%d\t%9g\t%9g\t%9g\n",
				i, *X++, *L++, *U++);
		}
	else {
		printf("i\t%s x(i)\n", when);
		for(i = 1; i <= n; i++)
			printf("%d\t%g\n", i, *X++);
		}
	}

 void
MAIN__(void)
{
	FILE *nl;
	fint LW, N, ierror;
	fint *ipivot;
	real *G, *W;
	real f;
	int i, nb;
	extern int xargc;
	extern char **xargv;
	char **av = xargv, buf[256], *stub;
	typedef struct { char *msg; int code; } Sol_info;
	static Sol_info solinfo[] = {
		{ "normal return", 0 },
		{ "bug: error in input pars", 500 },
		{ "more than maxfun evaluations", 400 },
		{ "line search did not find lower point (maybe OK)", 100 },
		{ 0, 501 },
		};

	ASL_alloc(ASL_read_fg);
	stub = getstub(&av, &Oinfo);
	nl = jac0dim(stub, (fint)strlen(stub));
	if (n_con) {
		fprintf(Stderr, "Ignoring %d constraints.\n", n_con);
		fflush(Stderr);
		}

	N = n_var;
	LUv = (real *)Malloc(N*(18*sizeof(real)) + N*sizeof(fint));
	Uvx = LUv + N;
	X0 = Uvx + N;
	G = X0 + N;
	W = G + N;
	LW = 14*N;
	ipivot = (fint *)(W + LW);

	fg_read(nl,0);

	for(i = nb = 0; i < N; i++) {
		if (LUv[i] > negInfinity)
			nb++;
		if (Uvx[i] < Infinity)
			nb++;
		}
	msglvl = amplflag ? -3 : 1;
	maxit = N >> 1;
	if (maxit > 50)
		maxit = 50;
	else if (maxit <= 0)
		maxit = 1;
	maxfun = 150 * N;
	eta = 0.25;
	stepmx = 10.;
	accrcy = 100. * macheps;
	xtol = sqrt(accrcy);
	if (getopts(av, &Oinfo))
		return;
	if (x0prt)
		show_x("initial", N, nb, X0, LUv, Uvx);
	ierror = 0;
	f = objval(nprob, X0, &ierror);
	if (ierror) {
		fprintf(Stderr,
			"%s: cannot compute initial objective value\n",
			Oinfo.sname);
		exit(1);
		}
	if (nb)
		lmqnbc_(&ierror, &N, X0, &f, G, W, &LW, sfun, LUv, Uvx,
			ipivot, &msglvl, &maxit, &maxfun, &eta,
			&stepmx, &accrcy, &xtol);
	else
		lmqn_(&ierror, &N, X0, &f, G, W, &LW, sfun,
			&msglvl, &maxit, &maxfun, &eta,
			&stepmx, &accrcy, &xtol);

	if (ierror == -1)
		ierror = 1;
	else if (ierror < 0 || ierror > 3) {
		Sprintf(solinfo[4].msg = buf + 100,
			"bug: unexpected ierror = %ld", (long)ierror);
		ierror = 4;
		}
	solve_result_num = solinfo[ierror].code;
	i = Sprintf(buf, "tn: %s.\nFinal f = ", solinfo[ierror].msg);
	i += g_fmtop(buf+i, f);
	i += Sprintf(buf+i, " (%d evaluations)", nfcall);
	write_sol(buf, X0, 0, &Oinfo);
	if (solprt)
		show_x("final", N, nb, X0, LUv, Uvx);

	}
