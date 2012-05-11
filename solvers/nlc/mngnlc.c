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

 static int maximize;
 static fint Nprob;
 static real objsign = 1.;

extern int dmngb_ ANSI((fint *n, real *d, real *x, real *b,
		  U_fp calcf, U_fp calcg,
		  fint *iv, fint *liv, fint *lv, real *v,
		  fint *uiparm, real *urparm, U_fp ufparm));

extern real feval_ ANSI((fint *nobj, fint *needfg, real *x, real *g));

static fint I1 = 1, I2 = 2;

 static void
#ifdef KR_headers
calcf(N, X, NF, F, UI, UR, UF)
	fint *N;
	real *X;
	fint *NF;
	real *F;
	fint *UI;
	real *UR;
	U_fp UF;
#else
calcf(fint *N, real *X, fint *NF, real *F, fint *UI, real *UR, U_fp UF)
#endif
{ *F = objsign * feval_(&Nprob, &I1, X, (real *)0); }

 static void
#ifdef KR_headers
calcg(N, X, NF, G, UI, UR, UF)
	fint *N;
	real *X;
	fint *NF;
	real *G;
	fint *UI;
	real *UR;
	U_fp UF;
#else
calcg(fint *N, real *X, fint *NF, real *G, fint *UI, real *UR, U_fp UF)
#endif
{
	real *Ge;

	feval_(&Nprob, &I2, X, G);
	if (maximize)
		for(Ge = G + *N; G < Ge; G++)
			*G = -*G;
	}

 extern real boundc_[], x0comn_[];
 extern fint funcom_[];

/* The following #defines are subscripts for argument iv to dmngb. */
/* It contains knobs and scratch space. */

#define mxfcal (17-1)
#define mxiter (18-1)
#define outlev (19-1)
#define prunit	(21-1)
#define solprt (22-1)
#define statpr	(23-1)
#define x0prt (24-1)

#undef n_con
#undef n_obj
#undef n_var
#undef nzc
#undef objtype
#undef LUrhs
#undef LUv
#undef X0
#define n_var	funcom_[0]
#define n_obj	funcom_[1]
#define n_con	funcom_[2]
#define nzc	funcom_[3]
#define objtype	(funcom_+5)
#define LUv	(boundc_+1)
#define X0 x0comn_

/* v[f] is the objective value: for a maximization problem,	*/
/* negate v[f] so rvmsg() will get the objective's sign right.	*/
#define f (10-1)

 void
MAIN__(VOID)
{
	fint N, liv, lv;
	fint *iv;
	real *D, *v;
	int i;
	char **av = xargv;
	static fint L2 = 2;

	ASL_alloc(ASL_read_f);
	progname = *av++;
	if (*av && !strcmp(*av,"-AMPL")) {
		amplflag = 1;
		av++;
		}

	if (n_con) {
		fprintf(Stderr, "Ignoring %d constraints.\n", n_con);
		fflush(Stderr);
		}

	N = n_var;
	liv = 59 + N;
	lv = 71 + N*(N + 21)/2;
	v = (real *)Malloc((N + lv)*sizeof(real) + liv*sizeof(fint));
	D = v + lv;
	iv = (fint *)(D + N);

	for(i = 0; i < N; i++)
		D[i] = 1.;

	divset_(&L2, iv, &liv, &lv, v);	/* set default iv and v values */
	iv[mxfcal] = iv[mxiter] = 1200;	/* increase defaults */

	if (amplflag) {
		/* Turn off printing (unless requested in $mng_options). */
		iv[outlev] = 0;
		/* iv[solprt] = 0; */
		iv[statpr] = 0;
		iv[x0prt] = 0;
		}
	 /* call vivvals to process command-line arguments */
	vivvals("mngnlc", "mng_options", av, iv, v);

	Nprob = nprob;
	if (nprob >= 0 && nprob < n_obj && objtype[nprob]) {
		maximize = 1;
		objsign = -1.;
		}

	dmngb_(&N, D, X0, LUv, (U_fp)calcf, (U_fp)calcg, iv, &liv, &lv, v,
		iv, v, (U_fp)calcf);
	}
#ifdef __cplusplus
	}
#endif
