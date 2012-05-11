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

/* Print an LP columnwise */

#include "asl.h"

main(int argc, char **argv)
{
	FILE *nl;
	int i, j, je;
	char *stub;
	ograd *og;
	ASL *asl;

	if (argc < 2) {
		printf("Usage: %s stub\n", argv[0]);
		return 1;
		}

	asl = ASL_alloc(ASL_read_f);
	stub = argv[1];
	nl = jac0dim(stub, (fint)strlen(stub));

	A_vals = (real *)Malloc(nzc*sizeof(real));

	f_read(nl,0);

	/* Print variable bounds */

	printf("\nVariable\tlower bound\tupper bound\n");
	for(i = 0; i < n_var; i++)
		printf("%8d\t%-8g\t%g\n", i+1,
			LUv[2*i], LUv[2*i+1]);

	/* Objectives */

	for(i = 0; i < n_obj; i++) {
		printf("\nObjective %d:\n", i+1);
		for(og = Ograd[i]; og; og = og->next)
			printf("\t%d\t%g\n", og->varno+1, og->coef);
		}

	if (!n_con)
		return 0;	/* bail out if no constraints */

	/* Print columns */

	for(i = 1; i <= n_var; i++) {
		printf("\nColumn %d:\n", i);
		j = A_colstarts[i-1];
		for(je = A_colstarts[i]; j < je; j++)
			printf("\t%d\t%g\n", A_rownos[j]+1, A_vals[j]);
		}

	/* Constraint lower, upper bounds */

	printf("\nConstraint\tlower bound\tupper bound\n");
	for(i = 0; i < n_con; i++)
		printf("%10ld\t%-8g\t%g\n", i+1, LUrhs[2*i], LUrhs[2*i+1]);
	return 0;
	}
