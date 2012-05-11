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
#ifdef __cplusplus
extern "C" {
#endif

extern void dn2g_ ANSI((fint *n, fint *p, real *x,
		  U_fp calcr, U_fp calcj,
		  fint *iv, fint *liv, fint *lv, real *v,
		  fint *uiparm, real *urparm, U_fp ufparm));

extern void dn2gb_ ANSI((fint *n, fint *p, real *x, real *b,
		   U_fp calcr, U_fp calcj,
		   fint *iv, fint *liv, fint *lv, real *v,
		   fint *uiparm, real *urparm, U_fp ufparm));

extern void ceval_ ANSI((fint *needfg, real *x, real *c, real *J));

static real *C0;

static fint I1 = 1, I2 = 2;

 static void
#ifdef KR_headers
calcr(N, P, X, NF, R, UI, UR, uf)
	fint *N;
	fint *P;
	real *X;
	fint *NF;
	real *R;
	fint *UI;
	real *UR;
	U_fp *uf;
#else
calcr(fint *N, fint *P, real *X, fint *NF, real *R, fint *UI, real *UR, U_fp *uf)
#endif
{
	real *Re;

	ceval_(&I1, X, R, (real *)0);
	for(Re = R + *N; R < Re; UR += 2)
		*R++ -= *UR;
	}

 static void
#ifdef KR_headers
calcj(N, P, X, NF, J, UI, UR, uf)
	fint *N;
	fint *P;
	real *X;
	fint *NF;
	real *J;
	fint *UI;
	real *UR;
	U_fp *uf;
#else
calcj(fint *N, fint *P, real *X, fint *NF, real *J, fint *UI, real *UR, U_fp *uf)
#endif
{ ceval_(&I2, X, C0, J); }

 extern real boundc_[], x0comn_[];
 extern fint funcom_[];

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

#undef n_con
#undef n_obj
#undef n_var
#undef nzc
#undef LUrhs
#undef LUv
#undef X0
#define n_var	funcom_[0]
#define n_obj	funcom_[1]
#define n_con	funcom_[2]
#define nzc	funcom_[3]
#define densej	funcom_[4]
#define X0 x0comn_

 void
MAIN__(VOID)
{
	fint *iv, liv, lv;
	fint N, NZ, P;
	real *rhsLU, *v;
	int i, j;
	char **av = xargv;
	static fint L1 = 1;
	real *LUv = boundc_ + 1, *LUrhs = LUv + 2*n_var;
	real Inf = boundc_[0], negInf = -Inf;

	Stderr_init_ASL(); /* Normally done by ASL_alloc(). */
	if (!densej && nzc < n_con*n_var) {
		fprintf(Stderr, "Use nlc -d\n");
		exit(1);
		}
	if (n_obj) {
		fprintf(Stderr, "Ignoring %d objectives.\n", n_obj);
		fflush(Stderr);
		}
	if (!n_con) {
		fprintf(Stderr, "No constraint!\n");
		exit(1);
		}

	ASL_alloc(ASL_read_f);
	progname = *av++;
	if (*av && !strcmp(*av,"-AMPL")) {
		amplflag = 1;
		av++;
		}

	N = n_con;
	P = n_var;
	NZ = nzc;
	liv = 82 + 4*P;
	lv = 105 + P*(N + 2*P + 20) + 2*N;
	v = (real *)Malloc((lv + P)*sizeof(real) + liv*sizeof(fint));
	iv = (fint *)(v + lv);

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

	divset_(&L1, iv, &liv, &lv, v);	/* set default iv and v values */
	if (amplflag) {
		iv[covprt] = 0;
		iv[covreq] = 0;
		iv[outlev] = 0;
		iv[rdreq] = 0;
		/*iv[solprt] = 0;*/
		iv[statpr] = 0;
		iv[x0prt] = 0;
		}
	 /* call vivvals to process command-line arguments */
	vivvals("nl2nlc", "nl2_options", av, iv, v);
	C0 = (real *)Malloc(N*sizeof(real));

	for(i = j = 0; i < N;  i++) {
		if (LUv[2*i] <= negInf) {
			LUv[2*i] = negInfinity;
			j++;
			}
		if (LUv[2*i+1] >= Inf) {
			LUv[2*i+1] = Infinity;
			j++;
			}
		}
	if (j < 2*n_var)
		dn2gb_(&N, &P, X0, LUv, (U_fp)calcr, (U_fp)calcj,
			iv, &liv, &lv, v, &NZ, LUrhs, (U_fp)calcr);
	else
		dn2g_(&N, &P, X0,       (U_fp)calcr, (U_fp)calcj,
			iv, &liv, &lv, v, &NZ, LUrhs, (U_fp)calcr);

	}
#ifdef __cplusplus
	}
#endif
