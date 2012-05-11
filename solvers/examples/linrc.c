/****************************************************************
Copyright (C) 1997-2000 Lucent Technologies
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

/* Print an LP rowwise, including names from .row and .col files. */

#include "asl.h"

 int
main(int argc, char **argv)
{
	FILE *nl;
	int bc, clen, dw, i, j, k, rlen;
	char gbuf[32];
	char *buf, *progname, *s, *stub;
	cgrad *cg, *cg0;
	ograd *og;
	real *b, *c, t;
	ASL *asl;

	dw = clen = rlen = 0;

	progname = *argv;
	if (argc < 2) {
 usage:
		fprintf(Stderr, "Usage: [-wnnn] %s stub\n\
	nnn = display width (overrides $display_width)\n", progname);
		return 1;
		}

	asl = ASL_alloc(ASL_read_f);
	stub = *++argv;
	if (*stub == '-')
		switch(stub[1]) {
		  case 'w':
			if (stub[2])
				stub += 2;
			else if (--argc < 2)
				goto usage;
			else
				stub = *++argv;
			dw = atoi(stub);
			if (--argc < 2)
				goto usage;
			stub = *++argv;
			break;
		  default:
			goto usage;
		}
	else if (s = getenv("display_width"))
		dw = atoi(s);
	if (dw < 2)
		dw = 79;

	nl = jac0dim(stub, (fint)strlen(stub));

	rlen = maxrownamelen;
	clen = maxcolnamelen;
	if (rlen == 0 || clen == 0) {
		printf("To get .row and .col files, in your AMPl session say\n");
		if (s = getenv("solver"))
			printf("\toption %s_auxfiles rc;\n%s", s,
				"before saying\n\tsolve;\nor\n");
		printf("\toption auxfiles rc;\nbefore saying\n\twrite ...\n");
		if (rlen == 0)
			rlen = strlen(con_name(n_con-1));
		if (clen == 0)
			clen = strlen(var_name(n_var-1));
		}

	f_read(nl,0);

	/* Print constant term in first objective if it is nonzero. */

	if (t = objconst(0)) {
		g_fmt(gbuf, t);
		printf("\nObjective adjustment (added constant) = %s\n", gbuf);
		}

	/* Print variable bounds and coefficients of first objective. */

	c = (real *)Malloc(n_var*sizeof(real));
	for(i = 0; i < n_var; i++)
		c[i] = 0;
	if (n_obj)
		for(og = Ograd[0]; og; og = og->next)
			c[og->varno] = og->coef;
	if (clen < 8)
		clen = 8;
	clen += 2;
	printf("\nVariable%*s%-16s%-16s%s\n",
		clen - 8, "", n_obj ? obj_name(0) : "", "lower bound",
		"upper bound");
	for(i = 0; i < n_var; i++)
		printf("%-*s%-16g%-16g%g\n", clen, var_name(i),
			c[i], LUv[2*i], LUv[2*i+1]);

	/* Print rows. */

	buf = (char *)Malloc(rlen + clen + 20);

	b = LUrhs;
	for(i = 0; i < n_con; i++, b += 2) {
		k = printf("\n%s:", con_name(i));
		if (k > 8)
			putchar('\n');
		if (b[0] < b[1] && b[0] > negInfinity) {
			k = printf("\t%g <=", b[0]) + 7;
			bc = 1;
			}
		else
			k = bc = dw;
		for(cg = cg0 = Cgrad[i]; cg; cg = cg->next) {
			t = cg->coef;
			if (cg == cg0)
				if (t == -1.) {
					j = 1;
					buf[0] = '-';
					t = -t;
					}
				else
					j = 0;
			else {
				j = 2;
				buf[1] = ' ';
				if (t < 0) {
					t = -t;
					buf[0] = '-';
					}
				else
					buf[0] = '+';
				}
			if (t != 1.) {
				g_fmtp(gbuf, t, 6);
				j += Sprintf(buf+j, "%s*", gbuf);
				}
			j += Sprintf(buf+j, var_name(cg->varno));
			if (j + bc + k < dw)
				k += printf("%*s%s", bc, "", buf);
			else if (bc > 1)
				k = printf("\t%s", buf) + 7;
			else
				k = printf("\n\t\t%s", buf) + 13;
			bc = 1;
			}
		if (b[1] < Infinity) {
			g_fmtp(gbuf, b[1], 6);
			j = Sprintf(buf, b[0] < b[1] ? "<= %s" : "= %s",
					gbuf);
			printf(j + k >= dw ? "\n\t\t%s" : " %s", buf);
			}
		printf("\n");
		}
	/* Give return code 1 when invoked by ampl to turn off */
	/* ampl's complaint about not reading a solution file. */
	/* Instead, ampl will say "exit code 1". */
	return (s = *++argv) && !strncmp(s, "-AMPL", 5) ? 1 : 0;
	}
