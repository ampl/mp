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

/* Print an LP rowwise */

#include "asl.h"

main(int argc, char **argv)
{
	FILE *nl;
	int i;
	char *fmt, *stub;
	cgrad *cg;
	ograd *og;
	real *b, *c;
	ASL *asl;

	if (argc < 2) {
		printf("Usage: %s stub\n", argv[0]);
		return 1;
		}

	asl = ASL_alloc(ASL_read_f);
	stub = argv[1];
	nl = jac0dim(stub, (fint)strlen(stub));

	f_read(nl,0);

	/* Print variable bounds and coefficients of first objective */

	c = (real *)Malloc(n_var*sizeof(real));
	for(i = 0; i < n_var; i++)
		c[i] = 0;
	if (n_obj)
		for(og = Ograd[0]; og; og = og->next)
			c[og->varno] = og->coef;
	printf("\nVariable\tlower bound\tupper bound\tcost\n");
	for(i = 0; i < n_var; i++)
		printf("%8d\t%-8g\t%-8g\t%g\n", i+1,
			LUv[2*i], LUv[2*i+1], c[i]);

	/* Print rows */

	b = LUrhs;
	for(i = 0; i < n_con; i++, b += 2) {
		printf("\nRow %d:", i+1);
		if (b[0] < b[1] && b[0] > negInfinity)
			printf(" %g <=", b[0]);
		fmt = " %g*X%d";
		for(cg = Cgrad[i]; cg; cg = cg->next, fmt = " + %g*X%d")
			printf(fmt, cg->coef, cg->varno+1);
		if (b[1] < Infinity) {
			printf(b[0] < b[1] ? " <= %g" : " == %g", b[1]);
			}
		printf("\n");
		}
	return 0;
	}
