/****************************************************************
Copyright (C) 2000 Lucent Technologies
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

Cextern int setulb_(fint *n, fint *m, real *x, real *l, real *u,
	fint *nbd, real *f, real *g, real *factr, real *pgtol, real *wa,
	fint *iwa, char *task, fint *iprint, char *csave, fint *lsave,
	fint *isave, real *dsave, ftnlen task_len, ftnlen csave_len);

 static real factr = 100.;
 static fint iprint = -1;
 static fint m = 10;
 static fint maxfun = 2000;
 static fint maxit = 1000;
 static int nprob = -1;
 static real pgtol = 1e-12;

 static char iprint_desc[] = "print level:\n\
		-1 ==> no printing (default)\n\
		 0 ==> 1 line at last iteration\n\
		 0 < iprint < 99 ==> f and |proj g| every iprint iters\n\
		 100 ==> also changes of active set and final x\n\
		 > 100 ==> details every iteration.\n\
		For iprint > 0, also write file \"iterate.dat\".";

static keyword keywds[] = {	/* must be in alphabetical order */
 KW("factr", D_val, &factr, "factr*macheps = tol. for rel. change in func. values; default 100."),
 KW("iprint", L_val, &iprint, iprint_desc),
 KW("m", L_val, &m, "number of corrections, 3 <= m <= 20 recommended; default = 10."),
 KW("maxfun", L_val, &maxfun, "maximum function and grad. evals (default 2000)."),
 KW("maxit", L_val, &maxit, "maximum iterations (default 1000)."),
 KW("nprob", I_val, &nprob, "objective choice: 1 (default) = 1st."),
 KW("outlev", L_val, &iprint, "synonym for iprint."),
 KW("pgtol", D_val, &pgtol, "tolerance for max(|proj g_i|), default 1e-12."),
 KW("version", Ver_val, 0, "report version"),
 KW("wantsol", WS_val, 0, WSu_desc_ASL+5)
 };

 static char lbfgvers[] = "AMPL/L-BFGS-B\0\nAMPL/LBFGSB Driver Version 20110622\n";

 static Option_Info Oinfo =
	{"lbfgsb", "L-BFGS-B", "lbfgsb_options", keywds, nkeywds, 1., lbfgvers,
	 0,0,0,0,0, 20110622};

 extern char **xargv;

 void
MAIN__(void)
{

	ASL *asl;
	FILE *nl;
	char buf[256], csave[60], *msg, *stub, task[60];
	fint i, isave[44], *iwa, lsave[4], n, *nbd, sr;
	int mlen, negate;
	real *L, *U, dsave[29], f, *g, *wa, *x;
	size_t ilen, wlen;

	asl = ASL_alloc(ASL_read_fg);
	stub = getstops(xargv, &Oinfo);
	nl = jac0dim(stub, (fint)strlen(stub));
	if (nprob < 0)
		nprob = 0;
	else
		--nprob;
	if (nprob < 0 || nprob >= n_obj) {
		sr = 101;
		msg = "No objective.";
		fclose(nl);
		goto done;
		}
	if (n_con) {
		printf("\nIgnoring %d constraints.\n", n_con);
		need_nl = 0;
		}
	if (m <= 0) {
		printf("\nAssuming m = 3\n");
		m = 3;
		}
	n = n_var;
	wlen = n*(2*m+5) + 12*m*(m + 1);
	ilen = 3*n;
	X0 = x = (real*)Malloc((4*n+wlen)*sizeof(real)
				+ (ilen+n)*sizeof(fint));
	g = x + n;
	L = LUv = g + n;
	U = Uvx = L + n;
	wa = U + n;
	iwa = (fint*)(wa + wlen);
	nbd = iwa + ilen;
	fg_read(nl,0);

	if ((negate = objtype[nprob]) && iprint >= 0)
		printf("\nFunction is to be maximized:\n\t%s%d.\n\n",
			"function values are negated in output for iprint=",
			iprint);

	for(i = 0; i < n; i++) {
		if (L[i] <= negInfinity)
			nbd[i] = U[i] >= Infinity ? 0 : 3;
		else
			nbd[i] = U[i] >= Infinity ? 1 : 2;
		}

	strcpy(task, "START");
	memset(task+5, ' ', sizeof(task)-5);
	msg = 0;
 loop:
	setulb_(&n, &m, x, L, U, nbd, &f, g, &factr, &pgtol, wa,
		iwa, task, &iprint, csave, lsave, isave, dsave,
		(ftnlen)sizeof(task), (ftnlen)sizeof(csave));
	switch(task[0]) {
	 case 'A':	/*ABNO*/
			sr = 500;
			break;

	 case 'C':	/*CONV*/
			sr = dsave[12] <= pgtol ? 0 : 100;
			break;

	 case 'F':	/*FG*/
			f = objval(nprob, x, 0);
			objgrd(nprob, x, g, 0);
			if (negate) {
				f = -f;
				for(i = 0; i < n; i++)
					g[i] = -g[i];
				}
			goto loop;

	 case 'N':	/*NEW_X*/
			if (isave[33] > maxfun) {
				sr = 400;
				msg = "more than maxfun evaluations";
				break;
				}
			if (isave[29] > maxit) {
				sr = 401;
				msg = "more than maxit iterations";
				break;
				}
			goto loop;
	 default:
			sr = 501;
	 }
 done:
	solve_result_num = sr;
	if (msg)
		mlen = strlen(msg);
	else {
		msg = task;
		for(i = sizeof(task); i > 0 && task[--i] == ' '; );
		mlen = i;
		}
	i = sprintf(buf, "%s: %.*s\n%ld iterations, %ld %s\n",
			Oinfo.bsname, mlen, msg, isave[29], isave[33],
			"function and gradient evaluations.");
	i += sprintf(buf+i, "Objective = %.*g.\n", obj_prec(),
			negate ? -f: f);
	i += sprintf(buf+i, "Projected gradient maxnorm = %g.\n", dsave[12]);
	write_sol(buf, x, 0, &Oinfo);
	}
