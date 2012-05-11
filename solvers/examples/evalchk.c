/****************************************************************
Copyright (C) 2009 David M. Gay

Permission to use, copy, modify, and distribute this software and its
documentation for any purpose and without fee is hereby granted,
provided that the above copyright notice appear in all copies and that
both that the copyright notice and this permission notice and warranty
disclaimer appear in supporting documentation.

The author disclaims all warranties with regard to this software,
including all implied warranties of merchantability and fitness.
In no event shall the author be liable for any special, indirect or
consequential damages or any damages whatsoever resulting from loss of
use, data or profits, whether in an action of contract, negligence or
other tortious action, arising out of or in connection with the use or
performance of this software.
****************************************************************/

/* Check for numerical errors in evaluating constraints and objectives
 * at the given starting point, setting suffix .numerr to 1 for those
 * that have errors.
 */

#include "getstub.h"

 static SufDecl
suftab[] = {
	{ "numerr", 0, ASL_Sufkind_con  | ASL_Sufkind_outonly },
	{ "numerr", 0, ASL_Sufkind_obj  | ASL_Sufkind_outonly },
	{ "numerr", 0, ASL_Sufkind_prob | ASL_Sufkind_outonly }
	};

 static keyword keywds[] = {	/* must be in alphabetical order */
	{ "version", Ver_val }
	};

 static char *usage[] = {
	"to set suffix .numerr on constaints and objectives",
	"to 0 if the entity can be evaluated at the given starting guess,",
	"to 1 if not, and to 2 if the entity but not its gradient can be",
	"evaluated.  For the current problem, .numerr is set to the number",
	"of objectives and constraints having nonzero .numerr values.",
	0};

 static Option_Info Oinfo = { "evalchk", 0, /*"evalchk_options"*/ 0,
				keywds, nkeywds, 0, "AMPL evaluation checker",
				usage,0,0,0,0, 20090712 };

 int
main(int argc, char **argv)
{
	ASL *asl;
	FILE *f;
	char *b, buf[256], *fmt, msg[256], *stub;
	fint bad;
	int i, j, nbad, nbadc, nbadcomm, nbado, nc, no, *z;
	real  *x;
	char fmtnoname[] = "Error evaluating defined variable %d\n";
	char fmtwithnm[] = "Error evaluating defined variable %d = %s\n";

	asl = ASL_alloc(ASL_read_fg);
	stub = getstops(argv, &Oinfo);
	f = jac0dim(stub, (fint)0);
	suf_declare(suftab, sizeof(suftab)/sizeof(SufDecl));
	X0 = x = M1zapalloc(n_var * sizeof(real));
	fg_read(f,0);
	nc = n_con;
	no = n_obj;
	z = M1zapalloc((nc+no+1)*sizeof(int));
	if (nc)
		suf_iput("numerr", ASL_Sufkind_con, z);
	bad = nbad = nbadc = nbadcomm = nbado = 0;
	asl->p.Xknown(asl, x, &bad);
	if (bad) {
		nbadcomm = 1;
		nbad = -1;
		z += nc;
		if (no) {
			suf_iput("numerr", ASL_Sufkind_obj, z);
			z += no;
			}
		if (i = cv_index) {
			fmt = fmtnoname;
			strcpy(stub_end, ".fix");
			j = 0;
			if (f = fopen(filename, "r")) {
				for(;;) {
#define					next_line fgets(buf,sizeof(buf),f)
					if (!next_line)
						goto eof;
					for(b = buf; *b; b++)
						if (*b == '=') {
							while(++j < i)
								if (!next_line)
									goto eof;
							b = buf;
							while(*b && *b != '=')
								b++;
							if (*b != '=' || b < buf + 2)
								j = 0;
							else
								b[-1] = 0;
							fmt = fmtwithnm;
							goto eof;
							}
					}
	 eof:
				fclose(f);
				}
			sprintf(msg, fmt, i, buf);
			}
		else
			sprintf(msg, "Surprise \"bad\" setting by Xknown -- bug?");
		}
	else {
		for(i = 0; i < nc; ++i, ++z) {
			conival(i,x,&bad);
			if (bad) {
				++nbadc;
				*z = bad;
				bad = 0;
				}
			}
		suf_iput("numerr", ASL_Sufkind_obj, z);
		for(i = 0; i < no; ++i, ++z) {
			objval(i,x,&bad);
			if (bad) {
				++nbado;
				*z = bad;
				bad = 0;
				}
			}
		if (nbad = nbadc + nbado) {
			if (nbadc == 0)
				sprintf(msg, "%d bad objective%s", nbad,
					nbad == 1 ? "" : "s");
			else if (nbado == 0)
				sprintf(msg, "%d bad constraint%s", nbad,
					nbad == 1 ? "" : "s");
			else
				sprintf(msg, "%d bad entities:\n\t"
					"%d constraint%s and %d objective%s",
					nbadc, nbadc == 1 ? "" : "s",
					nbado, nbado == 1 ? "" : "s");
			}
		}
	if (*z = nbad)
		solve_result_num = 500;
	suf_iput("numerr", ASL_Sufkind_prob, z);
	write_sol(msg, 0, 0, &Oinfo);
	return 0;
	}
