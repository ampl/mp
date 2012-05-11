/****************************************************************
Copyright (C) 1999 Lucent Technologies
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

/* Test and illustrate fg_write(): copy .nl files. */
/* Invoke with -? (or '-?') for a usage summary. */

/* See comments and code after "if (myflags & addstuff)" for */
/* an example of adding some new linear constraints and objectives. */

#include "asl.h"

enum Myflags {
	Avals = 1,
	verbose = 2,
	addstuff = 4,
	addiguess = 8
	};

 static int
process(int myflags, char *istub, char *ostub, int flags)
{
	ASL *asl;
	FILE *f;
	NewVCO nu, *nup;
	int i, n, rc;
	ograd *og, *og1;

	asl = ASL_alloc(ASL_read_fg);
	return_nofile = rc = 1;
	f = jac0dim(istub, 0);
	rc = 1;
	if (!f) {
		printf("Can open neither %s nor %s.nl.\n", istub, istub);
		goto done;
		}
	nup = 0;
	if (myflags & addstuff) {
		/* Illustrate adding some (silly) variables, constraints,
		 * and objectives: assuming original variables
		 * {i in 1..n_var} x[i] (with n_var >= 2), we add
		 * new variables
		 *	{i in 1 .. 3} xnew[i] >= 0, <= 1
		 * new constraints
		 *	sum{i in 1..3} xnew[i] = 1
		 * and
		 *	-2 <= -2.5*x[2] + sum{i in 1..3} (i+1)*xnew[i] <= 4
		 * and new objectives
		 *	minimize 1.2 + sum{i in 1..3} (10 + i)*xnew[i]
		 * and
		 *	maximize -3.5 + x[1] - sum{i in 1..3} (i+1)*xnew[i]
		 */
		static real x0[3] = { .2, .3, .5}, y0[2] = { -1.5, 3.75 };
		nup = &nu;
		nu.nnv = 3;
		nu.nnc = 2;
		nu.nno = 2;
		n = nu.nnc + nu.nno;
		nu.LUnc = (real*)M1alloc(n*sizeof(ograd*)
				+ (nu.nnv*n + 2)*sizeof(ograd)
				+ nu.nno
				+ (nu.nno + 2*(nu.nnv+nu.nnc))*sizeof(real));
		nu.Unc = nu.LUnc + nu.nnc;
		nu.LUnv = nu.Unc + nu.nnc;
		nu.Unv = 0;
		nu.oc = nu.LUnv + 2*nu.nnv;
		og1 = (ograd*)(nu.oc + nu.nno);
		nu.newc = (ograd**)(og1 + nu.nnv*n + 2);
		nu.newo = nu.newc + nu.nnc;
		nu.ot = (char*)(nu.newo + nu.nno);
		nu.newc[0] = og1;
		for(i = 0; i < 3; i++) {
			og = og1++;
			og->next = og1;
			og->varno = n_var + i;
			og->coef = 1.;
			}
		og->next = 0;
		nu.newc[1] = og = og1++;
		og->next = og1;
		og->varno = 1;
		og->coef = -2.5;
		for(i = 0; i < 3; i++) {
			og = og1++;
			og->next = og1;
			og->varno = n_var + i;
			og->coef = 2. + i;
			}
		og->next = 0;
		nu.ot[0] = 0;	/* minimize */
		nu.ot[1] = 1;	/* maximize */
		nu.oc[0] = 1.2;	/* constant term */
		nu.oc[1] = -3.5;
		nu.newo[0] = og1;
		for(i = 0; i < 3; i++) {
			og = og1++;
			og->next = og1;
			og->varno = n_var + i;
			og->coef = 11. + i;
			}
		og->next = 0;
		nu.newo[1] = og = og1++;
		og->next = og1;
		og->varno = 0;
		og->coef = 4.5;
		for(i = 0; i < 3; i++) {
			og = og1++;
			og->next = og1;
			og->varno = n_var + i;
			og->coef = -(2. + i);
			}
		og->next = 0;
		nu.LUnc[0] = nu.Unc[0] = 1.;
		nu.LUnc[1] = -2;
		nu.Unc[1] = 4.;
		for(i = 0; i < 6; ) {
			nu.LUnv[i++] = 0.;
			nu.LUnv[i++] = 1.;
			}
		nu.x0 = nu.d0 = 0;
		if (myflags & addiguess) {
			nu.x0 = x0;
			nu.d0 = y0;
			}
		}
	if (fg_wread(f, ASL_return_read_err | ASL_keep_all_suffixes))
		{
		*stub_end = 0;
		printf("Trouble reading %s.nl\n", filename);
		goto done;
		}
	n = strlen(ostub);
	if (n > 3 && !strcmp(ostub+n-3, ".nl"))
		n -= 3;
	switch(i = fg_write(ostub, nup, flags)) {
	  case 0:
		if (myflags & verbose)
			printf("%.*s written.\n", n, ostub);
		rc = 0;
		break;
	  case ASL_writeerr_openfail:
		printf("Could not open output file %.*s.nl\n", n, ostub);
		break;
	  case ASL_writeerr_badrops:
		printf("Bug: return ASL_writeerr_badrops from fg_write.\n");
		break;
	  default:
		printf("Bug: unexpected return %d from fg_write.\n", i);
		}
 done:
	ASL_free(&asl);
	return rc;
	}

 int
main(int argc, char **argv)
{
	char *b, *buf, *out_dir, *progname, *s, *t;
	int flags, i, myflags, n, rc;
	size_t L, L1, bs;

	progname = argv[0];
	flags = myflags = n = rc = 0;
	if (argc < 3) {
 usage:
		rc = 1;
 usage0:
		printf("Usage:\n\t%s [options] in_stub out_stub\nor\n\
	%s [options] -d out_dir in_stub [in_stub2 ...]\nor\n\
	%s [options] in_stub1 in_stub2 [in_stub3 ...] out_dir\n\noptions:\n\
	-A ==> read A_vals rather than cgrad structures\n\
	-g ==> write ASCII output (default = binary)\n\
	-k ==> keep derivative info\n\
	-r ==> \\r\\n in header\n\
	-x ==> add extra variables, constraints, objectives\n\
	-y ==> also add initial guesses (implies -x)\n",
			progname, progname, progname);
		return rc;
		}
	out_dir = 0;
 optloop:
	while(++n < argc && *(s = argv[n]) == '-')
		while(*++s) switch(*s) {
			case 'A':
				myflags |= Avals;
				break;
			case '-':
				++n;
				goto optsdone;
			case 'd':
				if (!*++s && !(s = argv[++n]))
					goto usage;
				out_dir = s;
				goto optloop;
			case 'g':
				flags |= ASL_write_ASCII;
				break;
			case 'k':
				flags |= ASL_keep_derivs;
				break;
			case 'r':
				flags |= ASL_write_CR;
				break;
			case 'v':
				myflags |= verbose;
				break;
			case 'x':
				myflags |= addstuff;
				break;
			case 'y':
				myflags |= addstuff | addiguess;
				break;
			case '?':
				goto usage0;
			default:
				goto usage;
			}
 optsdone:
	if (n >= argc)
		goto usage;
	if (out_dir) {
 have_out_dir:
		bs = L1 = strlen(out_dir) + 2;
		for(i = n; s = argv[i]; i++) {
			L = L1 + strlen(s);
			if (bs < L)
				bs = L;
			}
		buf = (char*)Malloc(bs);
		strcpy(buf, out_dir);
		buf[--L1 - 1] = '/';
		while(n < argc) {
			s = argv[n++];
			/* b = basename(s) */
			for(b = t = s; *t;)
				if (*t++ == '/')
					b = t;
			strcpy(buf+L1, b);
			rc |= process(myflags, s, buf, flags);
			}
		}
	else switch(argc - n) {
		case 1:
			goto usage;
		case 2:
			rc = process(myflags, argv[n], argv[n+1], flags);
			break;
		default:
			out_dir = argv[--argc];
			goto have_out_dir;
		}
	return rc;
	}
