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

#include "rvmsg.h"
#define asl cur_ASL

Cextern void dn2g_(fint *n, fint *p, real *x,
		  U_fp calcr, U_fp calcj,
		  fint *iv, fint *liv, fint *lv, real *v,
		  fint *uiparm, real *urparm, U_fp ufparm);

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

#define covprt	(14-1)
#define covreq	(15-1)
#define f	(10-1)
#define nfcall	(6-1)
#define ngcall	(30-1)
#define outlev	(19-1)
#define prunit	(21-1)
#define rdreq	(57-1)
#define solprt	(22-1)
#define statpr	(23-1)
#define x0prt	(24-1)

 void
MAIN__(void)
{
	FILE *nl;
	fint *iv, liv, lv;
	fint N, NZ, P;
	real *rhsLU, *v;
	int i, j;
	char **av = xargv, *stub;
	static fint L1 = 1;

	ASL_alloc(ASL_read_fg);
	stub = getstub(&av, &Oinfo);
	if (!stub)
		usage_ASL(&Oinfo, 1);

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
	if (amplflag) {
		iv[covprt] = 0;
		iv[covreq] = 0;
		iv[outlev] = 0;
		iv[rdreq] = 0;
		iv[solprt] = 0;
		iv[statpr] = 0;
		iv[x0prt] = 0;
		}
	vivvals("nl2", "nl2_options", av, iv, v);
	if (amplflag
		&& !iv[covprt]
		&& !iv[covreq]
		&& !iv[outlev]
		&& !iv[rdreq]
		&& !iv[solprt]
		&& !iv[statpr]
		&& !iv[x0prt])
		iv[prunit] = 0;

	for(i = j = 0; i < P;  i++)
		if (LUv[2*i] > negInfinity || LUv[2*i+1] < Infinity) {
			j = 1;
			break;
			}
	if (j)
		dn2gb_(&N, &P, X0, LUv, (U_fp)calcr, (U_fp)calcj,
			iv, &liv, &lv, v, &NZ, LUrhs, (U_fp)calcr);
	else
		dn2g_(&N, &P, X0,       (U_fp)calcr, (U_fp)calcj,
			iv, &liv, &lv, v, &NZ, LUrhs, (U_fp)calcr);

	write_sol(rvmsg("nl2",iv,v,1), X0, 0, &Oinfo);
	}
