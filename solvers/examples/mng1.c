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
Cextern int dmngb_(fint *n, real *d, real *x, real *b,
		  U_fp calcf, U_fp calcg,
		  fint *iv, fint *liv, fint *lv, real *v,
		  fint *uiparm, real *urparm, U_fp ufparm);

 static void
calcf(fint *N, real *X, fint *NF, real *F,
	fint *UI, real *UR, U_fp UF)
{
	fint nerror = 0;

	*F = objval(0, X, &nerror);
	if (nerror)
		*NF = 0;
	}

 static void
calcg(fint *N, real *X, fint *NF, real *G,
	fint *UI, real *UR, U_fp UF)
{
	fint nerror = 0;

	objgrd(0, X, G, &nerror);
	if (nerror)
		*NF = 0;
	}

 void
MAIN__(void)
{
	fint N, *iv, liv, lv;
	real *D, *v;
	int i;
	extern int xargc;
	extern char **xargv;
	char buf[64], *stub;
	FILE *nl;

	if (xargc < 2) {
		fprintf(Stderr, "usage: %s stub\n", xargv[0]);
		exit(1);
		}
	stub = xargv[1];
	ASL_alloc(ASL_read_fg);
	nl = jac0dim(stub, (fint)strlen(stub));
	if (n_con) {
		fprintf(Stderr, "Ignoring %d constraints.\n", n_con);
		fflush(Stderr);
		}

	X0 = (real *)Malloc(n_var*sizeof(real));
	fg_read(nl,0);

	N = n_var;
	liv = 59 + N;
	lv = 71 + N*(N + 21)/2;
	v = (real *)Malloc((N + lv)*sizeof(real) + liv*sizeof(fint));
	D = v + lv;
	iv = (fint *)(D + N);

	for(i = 0; i < N; i++)
		D[i] = 1.;

	iv[0] = 0; /* request all defaults */

	dmngb_(&N, D, X0, LUv, (U_fp)calcf, (U_fp)calcg, iv, &liv, &lv, v,
		iv, v, (U_fp)calcf);
	sprintf(buf, "mngb: return code %ld; final f = %.17g", iv[0], v[9]);
	write_sol(buf, X0, 0, 0);
	}
