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

#ifdef KR_headers
#include "varargs.h"
#else
#include "stdarg.h"
#endif

/* #include "nlp.h" */
#include "getstub.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef BSD
#define Len(x,y) (x, strlen(y))
#else
#define Len(x,y) x
#endif

static ASL *asl;

 static char *outstub, *outstub_end;

 static char *
#ifdef KR_headers
pv(x, buf) real x; char *buf;
#else
pv(real x, char *buf)
#endif
{
	int i, j;
	i = g_fmtp(buf, x, 12);
	for(j = 24 - i; i > 12; j += 12 - i)
		i = g_fmtp(buf, x, j);
	return buf;
	}

 static char *Basename;
 extern char *progname;
 static int objno = 1;

 static FILE *
#ifdef KR_headers
openf1(fname, mode) char *fname, *mode;
#else
openf1(char *fname, char *mode)
#endif
{
	FILE *rv;
	rv = fopen(fname, mode);
	if (!rv) {
		fprintf(Stderr, "%s: can't open %s\n", progname, fname);
		exit(2);
		}
	return rv;
	}

 static FILE *
#ifdef KR_headers
openfi(s, mode)
	char *s;
	char *mode;
#else
openfi(char *s, char *mode)
#endif
{
	strcpy(stub_end,s);
	return openf1(filename, mode);
	}

 static FILE *
#ifdef KR_headers
openfo(s, mode) char *s, *mode;
#else
openfo(char *s, char *mode)
#endif
{
	strcpy(outstub_end,s);
	return openf1(outstub, mode);
	}

 static int
#ifdef KR_headers
cglen(cg) register cgrad *cg;
#else
cglen(register cgrad *cg)
#endif
{
	register int k = 0;
	for(; cg; cg = cg->next)
		k++;
	return k;
	}

 static FILE *fd;

 static void
#ifdef KR_headers
aprintf(va_alist) va_dcl
#else
aprintf(const char *fmt, ...)
#endif
{
	char *s;
	char buf[32];
	va_list ap;
	int i, j;
	double x;

#ifdef KR_headers
	char *fmt;
	va_start(ap);
	fmt = va_arg(ap, char*);
#else
	va_start(ap, fmt);
#endif
	for(;;) {
		for(;;) {
			switch(i = *fmt++) {
				case 0:	  goto done;
				case '%': break;
				default:  putc(i,fd);
					  continue;
				}
			break;
			}
		switch(*fmt++) {
			case 'c':
				i = va_arg(ap, int);
				putc(i,fd);
				continue;
			case 'd':
				i = va_arg(ap, int);
				if (i < 0) {
					putc('-',fd);
					i = -i;
					}
				s = buf;
				do {
					j = i / 10;
					*s++ = i - 10*j + '0';
					}
					while(i = j);
				do {
					i = *--s;
					putc(i,fd);
					}
					while(s > buf);
				continue;
			case 'g':
				x = va_arg(ap, double);
				g_fmt(s = buf, x);
				while(i = *s++)
					putc(i,fd);
				continue;
			default:
				fprintf(Stderr, "aprintf bug: unexpect fmt %s\n",
					fmt-1);
				exit(1);
			}
		}
 done:
	va_end(ap);
	}

 static void
#ifdef KR_headers
bprintf(va_alist) va_dcl
#else
bprintf(const char *fmt, ...)
#endif
{
	char c;
	va_list ap;
	int i;
	double x;

#ifdef KR_headers
	char *fmt;
	va_start(ap);
	fmt = va_arg(ap, char*);
#else
	va_start(ap, fmt);
#endif
	if ((i = *fmt) != '%') {
		fmt++;
		c = i;
		fwrite(&c, 1, 1, fd);
		}

	for(;;) {
		while(*fmt == ' ')
			fmt++;
		if (*fmt++ != '%')
			break;
		switch(*fmt++) {
			case 'c':
				c = va_arg(ap, int);
				fwrite((char *)&c, 1, 1, fd);
				break;
			case 'd':
				i = va_arg(ap, int);
				fwrite((char *)&i, 1, sizeof(int), fd);
				break;
			case 'g':
				x = va_arg(ap, double);
				fwrite((char *)&x, 1, sizeof(double), fd);
				break;
			default:
				fprintf(Stderr, "bprintf bug: unexpect fmt %s\n",
					fmt-1);
				exit(1);
			}
		}
	va_end(ap);
	}

 static void (*pr)ANSI((const char *, ...));

 static void
#ifdef KR_headers
brout(c, n, LU)
	int c;
	int n;
	real *LU;
#else
brout(int c, int n, real *LU)
#endif
{
	int i;
	real L, U;
	pr("%c\n", c);
	for(i = 0; i < n; i++, LU += 2) {
		L = LU[0];
		U = LU[1];
		if (L <= negInfinity)
			if (U >= Infinity)
				pr("3\n", 3);
			else
				pr("1 %g\n", U);
		else
			if (U >= Infinity)
				pr("2 %g\n", L);
			else if (L == U)
				pr("4 %g\n", L);
			else
				pr("0 %g %g\n", L, U);
		}
	}

 static void
#ifdef KR_headers
nl_write(m, n, cg0, LU_bounds, LU_ranges, Max)
	int m;
	int n;
	cgrad **cg0;
	real *LU_bounds;
	real *LU_ranges;
	int Max;
#else
nl_write(int m, int n, cgrad **cg0, real *LU_bounds, real *LU_ranges, int Max)
#endif
{
	cgrad *cg, **cgt, *nextcg;
	int c, i, j, k, kg, nr;
	int *x;
	char *fmt;
	double *lu, *lue;
	size_t zaplen;

	if (binary_nl) {
		c = 'b';
		pr = bprintf;
		fmt = "wb";
		}
	else {
		c = 'g';
		pr = aprintf;
		fmt = "w";
		}

	fd = openfo(".nl", fmt);
	zaplen = (m+1)*sizeof(cgrad *);
	cgt = (cgrad **)Malloc(zaplen + n*sizeof(int));
	x = (int *)(cgt + m + 1);
	memset((char *)cgt, 0, zaplen);
	for(i = k = kg = 0; i < n; i++) {
		for(cg = cg0[i]; cg; cg = cg->next)
			if (--cg->varno < m)
				k++;
			else
				kg++;
		x[i] = k;
		}
	while(--i >= 0)
		for(cg = cg0[i]; cg; cg = nextcg) {
			nextcg = cg->next;
			cg->next = cgt[j = cg->varno];
			cgt[j] = cg;
			cg->varno = i;
			}
	lue = LU_ranges + 2*m;
	for(lu = LU_ranges, nr = 0; lu < lue; lu += 2)
		if (lu[0] > negInfinity
		 && lu[1] < Infinity
		 && lu[0] < lu[1])
			nr++;
	fprintf(fd, "%c3 0 0 0\n", c);
	fprintf(fd, " %d %d 1 %d\n 0 0\n 0 0\n 0 0\n 0 0 %d %d\n 0 0\n",
		n, m, nr, binary_nl ? Arith_Kind_ASL : 0, asl->i.flags);
	fprintf(fd, " %d %d\n 0 0\n 0 0 0 0 0\n", k, kg);
	for(i = 0; i < m; i++) {
		pr("C%d\n", i);
		pr("n%g\n", 0.);
		}
	pr("O%d %d\n", 0, 1-Max);
	pr("n%g\n", objconst(objno));
	brout('r', m, LU_ranges);
	brout('b', n, LU_bounds);
	pr("k%d\n", --n);
	for(i = 0; i < n; i++)
		pr("%d\n", x[i]);
	for(i = 0; i < m; i++) {
		cg = cgt[i];
		pr("J%d %d\n", i, cglen(cg));
		for(; cg; cg = cg->next)
			pr("%d %g\n", cg->varno, cg->coef);
		}
	cg = cgt[m];
	pr("G%d %d\n", 0, cglen(cg));
	for(; cg; cg = cg->next)
		pr("%d %g\n", cg->varno, cg->coef);
	fclose(fd);
	}

 static void
#ifdef KR_headers
mps(m, n, cg0, LU_bounds, LU_ranges)
	int m;
	int n;
	cgrad **cg0;
	real *LU_bounds;
	real *LU_ranges;
#else
mps(int m, int n, cgrad **cg0, real *LU_bounds, real *LU_ranges)
#endif
{
	FILE *f;
	int bounds, c, i, j, ranges, rhs;
	real *LU, *LUe;
	real L, U;
	cgrad **cg, *cgi;
	char buf[24], rbuf[16], *rname;

	f = openfo(".mps", "w");
	g_fmt_E = 'E';

	fprintf(f, "NAME          %s\nROWS\n", Basename);
	LU = LU_ranges;
	ranges = rhs = 0;
	for(i = 1; i <= m; i++, LU += 2) {
		if (LU[0] == LU[1]) {
			c = 'E';
			if (LU[0])
				rhs++;
			}
		else if (LU[0] <= negInfinity) {
			if (LU[1] >= Infinity)
				c = 'N';
			else {
				c = 'L';
				if (LU[1])
					rhs++;
				}
			}
		else {
			c = 'G';
			if (LU[0])
				rhs++;
			if (LU[1] < Infinity)
				ranges++;
			}
		fprintf(f, " %c  R%.4d\n", c, i);
		}
	j = m + 1;
	fprintf(f, " N  DUMMY\nCOLUMNS\n");
	cg = cg0;
	for(i = 1; i <= n; i++)
		for(cgi = *cg++; cgi; cgi = cgi->next) {
			c = cgi->varno;
			if (L = cgi->coef) {
				if (c == j)
					rname = "DUMMY   ";
				else
					sprintf(rname = rbuf, "R%-7.4d", c);
				fprintf(f, "    C%-7.4d  %s  %s\n", i, rname,
					pv(L,buf));
			}
		}
	fprintf(f, "RHS\n");
	if (rhs) {
		LU = LU_ranges;
		for(i = 1; i <= m; i++, LU += 2)
			if (LU[0] > negInfinity) {
				if (LU[0])
					fprintf(f, "    B         R%-7.4d  %s\n",
						i, pv(LU[0],buf));
				}
			else if (LU[1] < Infinity && LU[1])
				fprintf(f, "    B         R%-7.4d  %s\n",
					i, pv(LU[1],buf));
		}
	if (ranges) {
		fprintf(f, "RANGES\n");
		LU = LU_ranges;
		for(i = 1; i <= m; i++) {
			L = *LU++;
			U = *LU++;
			if (L < U && L > negInfinity && U < Infinity)
				fprintf(f, "    RANGE     R%-7.4d  %s\n",
					i, pv(U-L,buf));
			}
		}
	bounds = 0;
	LU = LU_bounds;
	for(LUe = LU + 2*n; LU < LUe; LU += 2)
		if (LU[0] || LU[1] < Infinity) {
			bounds = 1;
			break;
			}
	if (bounds) {
		fprintf(f, "BOUNDS\n");
		for(i = 1, LU = LU_bounds; LU < LUe; i++) {
			L = *LU++;
			U = *LU++;
			if (L <= negInfinity) {
				if (U >= Infinity) {
					fprintf(f, " FR BOUND     C%-7.4d\n", i);
					continue;
					}
				fprintf(f, " MI BOUND     C%-7.4d\n", i);
				fprintf(f, " UP BOUND     C%-7.4d  %s\n",
					i, pv(U,buf));
				continue;
				}
			if (L)
				fprintf(f, " LO BOUND     C%-7.4d  %s\n",
					i, pv(L,buf));
			if (U < Infinity)
				fprintf(f, " UP BOUND     C%-7.4d  %s\n",
					i, pv(U,buf));
			}
		}
	fprintf(f, "ENDATA\n");
	fclose(f);
	}

 static void
#ifdef KR_headers
usage(rc)
	int rc;
#else
usage(int rc)
#endif
{
	static char *ops[] = {
		"-b {write binary (-ob format) outstub.nl or outstub.sol}",
		"-g {write ASCII  (-og format) outstub.nl or outstub.sol}",
		"-m {write MPS format outstub.mps}",
		"-onn {use objective nn: 1 = first; 0 = no objective}",
		"-u {read stub.sol, write outstub.sol}",
		0};
	char **o = ops;
	fprintf(Stderr, "usage: %s [options] stub [outstub]\noptions:\n",
		progname);
	while(*o)
		fprintf(Stderr, "\t%s\n", *o++);
	exit(rc);
	}

 struct dvthead {
	real vbtol;
	fint Options[10];
	int m;
	int n;
	int nextra;
	int maxobj;
	int binary;
	};
 typedef struct dvthead dvthead;

 static void
#ifdef KR_headers
dfmt(op)
	int op;
#else
dfmt(int op)
#endif
{
	cgrad **cg0, **cg01, *cgi, *cgi0, *cgi00, *cgj, **cgnew,
		**cgprev, **cgx, *free_cg;
	ograd *og;
	real *LU, *LUdv, *LUdvi, *LUdrhs, *LUdvxi, *LUrhse, *LUve;
	int Max, aextra, i, m, me, n, n1, nextra;
	FILE *f;
	real t;
	char *dvtype, *dvt;
	static char *AXIN[2] = { "AX", "IN" };
	dvthead dvth;
	char buf[32];

	n = c_vars;
	if (n < o_vars)
		n = o_vars;
	n1 = n + 1;
	m = n_con;
	LUrhse = LUrhs + 2*m;
	LUve = LUv + 2*n;
	aextra = nextra = 0;;
	for (LU = LUrhs, cgx = Cgrad; LU < LUrhse; LU += 2, cgx++)
		if (LU[0] > negInfinity) {
			if (LU[0])
				aextra++;
			if (LU[0] < LU[1] && LU[1] < Infinity) {
				nextra++;
				if (LU[1])
					aextra++;
				for(cgi = *cgx; cgi; cgi = cgi->next)
					aextra++;
				}
			else if (LU[0])
				aextra++;
			}
		else if (LU[1])
			aextra++;
	for (LU = LUv; LU < LUve; LU += 2) {
		if (LU[0] > negInfinity) {
			aextra++;
			nextra++;
			if (LU[0])
				aextra++;
			}
		if (LU[1] < Infinity) {
			aextra++;
			nextra++;
			if (LU[1])
				aextra++;
			}
		}
	me = m + nextra;
	LUdvi = LUdv = (real *)Malloc((me + n)*2*sizeof(real)
				+ me*sizeof(cgrad *)
				+ aextra*sizeof(cgrad)
				+ m);
	LUdvxi = LUdv + 2*m;
	LUdrhs = LUdvxi + 2*nextra;
	free_cg = (cgrad *)(LUdrhs + 2*n);
	cg0 = cg01 = (cgrad **)(free_cg + aextra);
	cgnew = cg0 + m;
	dvt = dvtype = (char *)(cgnew + nextra);
	for (LU = LUrhs, cgx = Cgrad; LU < LUrhse; LU += 2, cgx++) {
		cgi0 = 0;
		for(cgi = cgi00 = *cgx; cgi; cgi = cgj) {
			cgj = cgi->next;
			cgi->next = cgi0;
			cgi0 = cgi;
			cgi->varno++;
			}
		*cg01++ = cgi0;
		if (LU[0] > negInfinity) {
			*LUdvi++ = LU[0] == LU[1] ? negInfinity : 0;
			*LUdvi++ = Infinity;
			if (LU[0] < LU[1] && LU[1] < Infinity) {
				*LUdvxi++ = 0;
				*LUdvxi++ = Infinity;
				cgprev = cgnew++;
				for(cgi = *cgx; cgi; cgi = cgi->next) {
					*cgprev = cgj = free_cg++;
					cgprev = &cgj->next;
					cgj->varno = cgi->varno;
					cgj->coef = -cgi->coef;
					}
				if (LU[1]) {
					cgi = *cgprev = free_cg++;
					cgi->varno = n1;
					cgi->coef = -LU[1];
					cgprev = &cgi->next;
					}
				*cgprev = 0;
				*dvt++ = 2;
				}
			else
				*dvt++ = LU[0] == LU[1] ? 3 : 0;
			if (LU[0]) {
				cgi = cgi00->next = free_cg++;
				cgi->varno = n1;
				cgi->coef = LU[0];
				cgi->next = 0;
				}
			}
		else {
			*dvt++ = 1;
			*LUdvi++ = 0;
			*LUdvi++ = Infinity;
			for(cgi = cgi0; cgi; cgi = cgi->next)
				cgi->coef = -cgi->coef;
			if (LU[1]) {
				cgi = cgi00->next = free_cg++;
				cgi->varno = n1;
				cgi->coef = -LU[1];
				cgi->next = 0;
				}
			}
		}
	for (LU = LUv, i = 1; LU < LUve; LU += 2, i++) {
		if (LU[0] > negInfinity) {
			*LUdvxi++ = 0;
			*LUdvxi++ = Infinity;
			*cgnew++ = cgi = free_cg++;
			cgi->varno = i;
			cgi->coef = 1;
			if (LU[0]) {
				cgi = cgi->next = free_cg++;
				cgi->varno = n1;
				cgi->coef = LU[0];
				}
			cgi->next = 0;
			}
		if (LU[1] < Infinity) {
			*LUdvxi++ = 0;
			*LUdvxi++ = Infinity;
			*cgnew++ = cgi = free_cg++;
			cgi->varno = i;
			cgi->coef = -1;
			if (LU[1]) {
				cgi = cgi->next = free_cg++;
				cgi->varno = n1;
				cgi->coef = -LU[1];
				}
			cgi->next = 0;
			}
		}
	memset(LUdrhs, 0, n*2*sizeof(real));
	if (objno >= 0)
		for(og = Ograd[objno]; og; og = og->next) {
			LU = LUdrhs + 2*og->varno;
			LU[0] = LU[1] = og->coef;
			}
	if (Max = objtype[objno])
		for(LU = LUdv; LU < LUdrhs; LU += 2) {
			t = LU[0];
			LU[0] = -LU[1];
			LU[1] = -t;
			}

	/* Negate columns with lower bound -Infinity, finite upper bound. */
	/* This shouldn't be necessary, but shortens the MPS file */
	/* and may avoid bugs in some solvers. */
	for(cg01 = cg0, LU = LUdv; LU < LUdrhs; LU += 2, cg01++)
		if (LU[0] <= negInfinity && LU[1] < Infinity) {
			t = LU[0];
			LU[0] = -LU[1];
			LU[1] = -t;
			for(cgi = *cg01; cgi; cgi = cgi->next)
				cgi->coef = -cgi->coef;
			}

	if (op != 'm') {
		switch(op) {
			case 'b':
				binary_nl = 1;
				break;
			case 'g':
				binary_nl = 0;
			}
		nl_write(n, me, cg0, LUdv, LUdrhs, Max);
		}
	else {
		mps(n, me, cg0, LUdv, LUdrhs);
		f = openfo(".spc", "w");
		fprintf(f, "BEGIN %s\nROWS %d\nCOLUMNS %d\n",
			Basename, n + 1, me+1);
		fprintf(f, "*true value: COLUMNS %d\n", me);
		fprintf(f, "ELEMENTS %d\nM%sIMIZE\nOBJECTIVE DUMMY\n",
			nzc + aextra, AXIN[Max]);
		fprintf(f, "END %s\n", Basename);
		fclose(f);
		f = openfo(".adj", "w");
		g_fmt(buf, objconst(objno));
		fprintf(f, "'objective' %s\n", buf);
		fclose(f);
		}
	f = openfo(".duw", "wb");
	for(i = 0; i < 10; i++)
		dvth.Options[i] = ampl_options[i];
	dvth.vbtol = ampl_vbtol;
	dvth.m = m;
	dvth.n = n;
	dvth.nextra = nextra;
	dvth.maxobj = Max;
	dvth.binary = binary_nl;
	fwrite(&dvth, sizeof(dvthead), 1, f);
	fwrite(dvtype, m, 1, f);
	fclose(f);
	}

 static void
#ifdef KR_headers
fread1(x, siz, f, what)
	void *x;
	size_t siz;
	FILE *f;
	char *what;
#else
fread1(void *x, size_t siz, FILE *f, char *what)
#endif
{
	if (fread(x, siz, 1, f) != 1) {
		fprintf(Stderr, "Failure reading %s\n", what);
		exit(1);
		}
	}

 static int
#ifdef KR_headers
undo(s, op)
	char *s;
	int op;
#else
undo(char *s, int op)
#endif
{
	char *dve, *dvt, *dvtype, *msg;
	dvthead dvth;
	int i;
	fint L;
	double *x, *xi, *y, *yi, *yext;
	FILE *f;
	Option_Info Oinfo;

	L = strlen(s);
	filename = (char *)Malloc(L + 5);
	stub_end = filename + L;
	strcpy(filename, s);
	f = openfi(".duw", "rb");
	fread1(&dvth, sizeof(dvthead), f, filename);

	yext = (real *)Malloc(dvth.m + dvth.nextra*sizeof(real));
	dvtype = (char *)(yext + dvth.nextra);

	fread1(dvtype, dvth.m, f, filename);
	fclose(f);

	n_var = dvth.m + dvth.nextra;
	n_con = dvth.n;
	if (!(msg = read_soln(&y, &x)))
		return 1;
	dve = dvtype + dvth.m;
	yi = y;
	yext = y + dvth.m;
	for(dvt = dvtype; dvt < dve; dvt++, yi++)
		switch(*dvt) {
			case 1:
				*yi = -*yi;
				break;
			case 2:
				if (*yext > *yi)
					*yi = -*yext;
				yext++;
			}
	if (dvth.maxobj) {
		for(xi = y, dvt = dvtype; dvt < dve; xi++, dvt++)
			if (*dvt != 3)
				*xi = -*xi;
		}
	asl->i.n_var0 = n_var = dvth.n;
	asl->i.n_con0 = n_con = dvth.m;
	filename = outstub;
	stub_end = outstub_end;
	switch(op) {
		case 'b': binary_nl = 1; break;
		case 'g': binary_nl = 0; break;
		default:  binary_nl = dvth.binary;
		}
	for(i = 0; i < 10; i++)
		ampl_options[i] = dvth.Options[i];
	ampl_vbtol = dvth.vbtol;
	Oinfo.wantsol = 9;	/* suppress message echo, force .sol writing */
	i = solve_result_num;
	if (i >= 200 && i < 400) {
		/* interchange "infeasible" and "unbounded" */
		if (i >= 300)
			i -= 100;
		else
			i += 100;
		solve_result_num = i;
		}
	write_sol(msg, x, y, &Oinfo);
	return 0;
	}

#ifdef __cplusplus
	}
#endif

 int
#ifdef KR_headers
main(argc, argv) int argc; char **argv;
#else
main(int argc, char **argv)
#endif
{
	char *s, *stub;
	fint stublen;
	FILE *nl;
	static int Undo, haveo, options = 1, outform = -1;

	progname = *argv;
	asl = ASL_alloc(ASL_read_f);
 top:
	if (--argc <= 0)
		usage(1);
	s = *++argv;
	if (options && *s == '-') for(;;) {
		switch(*++s) {
			case 0:
				goto top;
			case '?':
				usage(0);
			case '-':
				options = 0;
				goto top;
			case 'b':
			case 'g':
			case 'm':
				outform = *s;
				continue;
			case 'o':
				haveo = 1;
				objno = atoi(++s);
				goto top;
			case 'u':
				Undo = 1;
				continue;
			default:
				usage(1);
			}
		}
	stub = s;
	if (argc > 1) {
		--argc;
		s = *++argv;
		}
	if (argc != 1)
		usage(1);
	stublen = strlen(s);
	outstub = (char *)Malloc(stublen + 5);
	outstub_end = outstub + stublen;
	strcpy(outstub, s);
	Basename = basename(s);
	if (Undo)
		return undo(stub, outform);
	else {
		nl = jac0dim(stub, (fint)strlen(stub));
		if (haveo && (objno < 0 || objno > n_obj)) {
			fprintf(Stderr, "-onn must specify 0 <= nn <= %d, not %d\n",
				n_obj, objno);
			exit(1);
			}
		objno--;
		f_read(nl,0);
		dfmt(outform);
		}
	return 0;
	}
