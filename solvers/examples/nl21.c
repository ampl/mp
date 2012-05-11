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

#include "asl.h"
#define asl cur_ASL

typedef void (*U_fp)(void);

Cextern void divset_(fint *L1, fint *IV, fint *LIV, fint *LV, real *V);

Cextern void dn2gb_(fint *n, fint *p, real *x, real *b,
		   U_fp calcr, U_fp calcj,
		   fint *iv, fint *liv, fint *lv, real *v,
		   fint *uiparm, real *urparm, U_fp ufparm);

 static void
calcr(fint *N, fint *P, real *X, fint *NF, real *R,
	fint *UI, real *UR, U_fp *uf)
{
	real *Re;
	fint nerror = 0;

	conval(X, R, &nerror);
	if (nerror)
		*NF = 0;
	else
		for(Re = R + *N; R < Re; UR += 2)
			*R++ -= *UR;
	}

 static void
calcj(fint *N, fint *P, real *X, fint *NF, real *J,
	fint *UI, real *UR, U_fp *uf)
{
	fint nerror = 0;

	jacval(X, J, &nerror);
	if (nerror)
		*NF = 0;
	}

/* The following #defines are subscripts for the v and iv arrays */
/* (of control information and scratch space) passed to dn2gb.   */

#define f (10-1)
#define nfcall (6-1)
#define ngcall (30-1)
#define prunit (21-1)

 void
MAIN__(void)
{
	FILE *nl;
	fint *iv, liv, lv;
	fint N, NZ, P;
	real *rhsLU, *v;
	int i, j;
	extern int xargc;
	extern char **xargv;
	char *stub;
	static fint L1 = 1;
	static char *rvmsg[9] = {
		"X-Convergence", /* 3 */
		"Relative Function Convergence",
		"X- and Relative Function Convergence",
		"Absolute Function Convergence",
		"Singular Convergence",
		"False Convergence",
		"Function Evaluation Limit",
		"Iteration Limit",
		"Unexpected return code"
		};
	char buf[256];

	if (xargc < 2) {
		fprintf(Stderr, "usage: %s stub\n", xargv[0]);
		exit(1);
		}
	stub = xargv[1];

	ASL_alloc(ASL_read_fg);

	amplflag = xargc >= 3 && !strncmp(xargv[2], "-AMPL", 5);

	nl = jac0dim(stub, (fint)strlen(stub));
	if (n_obj) {
		fprintf(Stderr, "Ignoring %d objectives.\n", n_obj);
		fflush(Stderr);
		}

	N = n_con;
	P = n_var;
	NZ = nzc;
	liv = 82 + 4*P;
	lv = 105 + P*(N + 2*P + 21) + 2*N;
	v = (real *)Malloc((lv + P)*sizeof(real) + liv*sizeof(fint));
	X0 = v + lv;
	iv = (fint *)(X0 + P);

	fg_read(nl,0);

	/* Check for valid problem: all equality constraints. */

	for(i = j = 0, rhsLU = LUrhs; i++ < N; rhsLU += 2)
		if (rhsLU[0] != rhsLU[1]) {
			if (j++ > 4) {
				/* Stop chattering if > 4 errors. */
				fprintf(Stderr, "...\n");
				exit(2);
				}
			fprintf(Stderr, "Lrhs(%d) = %g < Urhs(%d) = %g\n",
				i, rhsLU[0], i, rhsLU[1]);
			}
	if (j)
		exit(2);

	dense_j();	/* Tell jacval_ we want a dense Jacobian. */

	divset_(&L1, iv, &liv, &lv, v);	/* set default iv and v values */
	if (amplflag)
		iv[prunit] = 0; /* Turn off printing . */

	dn2gb_(&N, &P, X0, LUv, (U_fp)calcr, (U_fp)calcj,
		iv, &liv, &lv, v, &NZ, LUrhs, (U_fp)calcr);

	j = iv[0] >= 3 && iv[0] <= 10 ? (int)iv[0] - 3 : 8;
	i = Sprintf(buf, "nl21: %s", rvmsg[j]);
	if (j == 8)
		i += Sprintf(buf+i, " %ld", iv[0]);
	i += Sprintf(buf+i,
		"\n%ld function, %ld gradient evaluations",
		iv[nfcall], iv[ngcall]);
	i += Sprintf(buf+i, "\nFinal sum of squares = ");
	g_fmtop(buf+i, 2*v[f]);
	write_sol(buf, X0, 0, 0);
	}
