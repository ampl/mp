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

#include "asl_pfgh.h"
#define Std_dev	  /* offer std. dev. with inverse Hessian as cov. matrix */
#include "rvmsg.h"
#define asl cur_ASL

 static int maximize;
 static real objsign = 1., *obj_weight;

Cextern int dmnhb_(fint *n, real *d, real *x, real *b,
		  U_fp calcf, U_fp calcgh,
		  fint *iv, fint *liv, fint *lv, real *v,
		  fint *uiparm, real *urparm, U_fp ufparm);

#ifdef Std_dev
static void stddev(fint, real*, real*, real*, real*);
#endif

 static void
calcf(fint *N, real *X, fint *NF, real *F,
	fint *UI, real *UR, U_fp UF)
{
	fint nerror = 0;

	*F = objsign * objval(nprob, X, &nerror);
	if (nerror)
		*NF = 0;
	}

 static void
calcgh(fint *N, real *X, fint *NF, real *G, real *H,
	fint *UI, real *UR, U_fp UF)
{
	fint nerror = 0;
	real *Ge;

	objgrd(nprob, X, G, &nerror);
	if (nerror) {
		*NF = 0;
		return;
		}
	if (maximize) {
		duthes(H, -1, obj_weight, 0);
		for(Ge = G + *N; G < Ge; G++)
			*G = -*G;
		}
	else
		duthes(H, nprob, 0, 0);
	}

/* The following #defines are subscripts for argument iv to dmnhb. */
/* It contains knobs and scratch space. */

#define g0	(48-1)
#define mxfcal	(17-1)
#define mxiter	(18-1)
#define outlev	(19-1)
#define prunit	(21-1)
#define solprt	(22-1)
#define statpr	(23-1)
#define x0	(43-1)
#define x0prt	(24-1)
#define H	(56-1)

/* v[f] is the objective value: for a maximization problem,	*/
/* negate v[f] so rvmsg() will get the objective's sign right.	*/
#define f (10-1)

 void
MAIN__(void)
{
	FILE *nl;
	fint N, liv, lv;
	fint *iv;
	real *D, *v;
	int i;
	char **av = xargv, *stub;
	static fint L2 = 2;

	ASL_alloc(ASL_read_pfgh);
	stub = getstub(&av, &Oinfo);
	if (!stub)
		usage_ASL(&Oinfo,1);

	nl = jac0dim(stub, (fint)strlen(stub));
	if (n_con) {
		fprintf(Stderr, "Ignoring %d constraints.\n", n_con);
		fflush(Stderr);
		}

	N = n_var;
	liv = 59 + 3*N;
	lv = 78 + N*(N + 15);
	v = (real *)Malloc((2*N + lv)*sizeof(real) + liv*sizeof(fint));
	D = v + lv;
	X0 = D + N;
	iv = (fint *)(X0 + N);

	pfgh_read(nl, ASL_findOgroups);

	for(i = 0; i < N; i++)
		D[i] = 1.;

	divset_(&L2, iv, &liv, &lv, v);	/* set default iv and v values */
	iv[mxfcal] = iv[mxiter] = 1200;	/* increase defaults */

	if (amplflag) {
		/* Turn off printing (unless requested in $mnh_options). */
		iv[outlev] = 0;
		iv[solprt] = 0;
		iv[statpr] = 0;
		iv[x0prt] = 0;
		}
	vivvals("mnh", "mnh_options", av, iv, v);
	if (amplflag
		&& !iv[outlev]
		&& !iv[solprt]
		&& !iv[statpr]
		&& !iv[x0prt])
		iv[prunit] = 0;	/* Finish turning off printing. */

	if (nprob >= 0 && nprob < n_obj && objtype[nprob]) {
		obj_weight = (real *)Malloc(n_obj*sizeof(real));
		memset(obj_weight, 0, n_obj*sizeof(real));
		maximize = 1;
		obj_weight[nprob] = objsign = -1.;
		}

	dmnhb_(&N, D, X0, LUv, (U_fp)calcf, (U_fp)calcgh, iv, &liv, &lv, v,
		iv, v, (U_fp)calcf);

	v[f] *= objsign;	/* get objective sign right */
	write_sol(rvmsg("mnh",iv,v,1), X0, 0, &Oinfo);
#ifdef Std_dev
	if ((want_stddev || stddev_file) && iv[0] <= 6) {
		real *v0 = v - 1;
		/* Use the "g0" and "x0" arrays for scratch. */
		stddev(N, v0 + iv[H], X0, v0 + iv[x0], v0 + iv[g0]);
		}
#endif
	}


#ifdef Std_dev

#undef f

/* Some Port library routines... */
/* The #defines map names in the TOMS papers to PORT names. */
#define livmul_ dl7ivm_
#define lsqrt_ dl7srt_
#define v2norm_ dv2nrm_
Cextern void livmul_(fint* ,real *x, real *L, real *y);	/* solve L*x = y */
Cextern void lsqrt_(fint*,fint*,real*,real*,fint*); /* Cholesky factor */
Cextern real v2norm_(fint*, real*); /* compute 2-norm of a vector */

 static void
devout(FILE *F, fint N, real *x, real *z)
{
	fint i;

	fprintf(F, "var\t\tvalue\t\tstddev\n");
	for(i = 0; i < N; i++)
		fprintf(F, "%-8s\t%-8.6g\t%.6g\n", var_name(i), *x++, *z++);
	}

 static void
stddev(fint N, real *h, real *x, real *y, real *z)
{
	fint i, j;
	static fint I1 = 1;
	FILE *F;

	lsqrt_(&I1, &N, h, h, &i); /* overwrite h with Cholesky factor */
	if (i) {
		fprintf(Stderr, "indefinite Hessian\n");
		return;
		}
	for(i = 0; i < N; i++)
		y[i] = 0.;
	while(--i >= 0) {
		y[i] = 1.;
		livmul_(&N, y, h, y);
		j = N - i;
		z[i] = v2norm_(&j, y+i);
		for(j = i; j < N; j++)
			y[j] = 0.;
		}
	if (stddev_file && (F = fopen(stddev_file, "w")))
		devout(F, N, x, z);
	if (want_stddev)
		devout(stdout, N, x, z);
	}
#endif
