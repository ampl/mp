/****************************************************************
Copyright (C) 1997-1999 Lucent Technologies
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

#include "getstub.h"

 static ASL *asl;

 static fint L[17] = { 2, 0,0,0,0, 3, 1, 0,0,0,0,0, 0, 1, -1, -1, -1 };
 static real DA[4] = { 1., 1., 1e-5, -1. };
 static real objsign = 1.;
 static fint *B, *Beq, *Bne;
 static fint *c_ef, *c_eg, *c_egradf, *c_egradg, *c_egradh, *c_eh;
 static int Nb, maximize, namelen;
 static int nprob = 1;
 static char *stub;

#ifdef __cplusplus
 extern "C" {
#endif

 extern void setup1_
	ANSI((fint *b, real *da, fint *L, char *name, real *x0a, fint namelen));
 extern void setup2_ ANSI((fint*, real*, fint*, real*));
 extern void donlp2_(VOID);
 extern void ef_ ANSI((real *x, real *fx));
 extern void egradf_ ANSI((real *x, real *g));
 extern void eg_ ANSI((fint *i, real *x, real *fx));
 extern void egradg_ ANSI((fint *i, real *x, real *g));
 extern void eh_ ANSI((fint *i, real *x, real *fx));
 extern void egradh_ ANSI((fint *i, real *x, real *g));
 extern void endinf_(VOID);
 extern void errset_ ANSI((fint*));
 extern fint mmax_ ANSI((void));
 extern fint nmax_ ANSI((void));
 extern void outinf_ ANSI((fint *ninfo, real *dinfo, real *x, real *u));

 extern char **xargv;	/* Arguments to Fortran main() */

#ifdef __cplusplus
	}
#endif

 static keyword
keywds[] = {	/* must be in alphabetical order */
 KW("del0", D_val, DA+0, "scaled inequalities <= del0 are initially \"active\"; default 1"),
 KW("epsphi", D_val, DA+3, "tolerance for small penalty function differences; default 1e3*macheps"),
 KW("epsx", D_val, DA+2, "stopping tolerance on KKT conditions; default 1e-5"),
 KW("intakt", L_val, L+14, "\"interactive\" output to screen"),
 KW("maxfwd", IA_val, voffset_of(ASL,p.maxfwd_), "# of partials to forward recur; default = 5"),
 KW("maxit", L_val, L+15, "iteration limit"),
 KW("meu", L_val, L+0, "Fortran unit for special events"),
 KW("nprob", I_val, &nprob, "objective number (1 = first)"),
 KW("nreset", L_val, L+4, "limit on number of small steps; default n (# of variables)"),
 KW("numsm", L_val, L+16, "limit on small penalty function differences; default max(n,10)"),
 KW("objno", I_val, &nprob, "synonym for nprob"),
 KW("outlev", L_val, L+12, "1 = summarize iters; bits 2,4,8 increase file output"),
 KW("prou", L_val, L+5, "Fortran unit for \"normal\" output"),
 KW("silent", L_val, L+6, "1 (default) = no output; 0 = write to meu and prou"),
 KW("tau0", D_val, DA+1, "bound on sum of constraint violations for evaluations; default 1"),
 KW("wantsol", WS_val, 0, WS_desc_ASL+5)
 };

 static char donlpvers[] = "AMPL/DONLP2\0\nAMPL/DONLP2 Driver Version 20020506\n";

 static Option_Info Oinfo = { "donlp2", "DONLP2", "donlp2_options", keywds, nkeywds,
				1, donlpvers, 0,0,0,0,0, 20020506 };
 static void
#ifdef KR_headers
nbcount(lu, n, NE, NNE) real *lu; int n; fint *NE; fint *NNE;
#else
nbcount(real *lu, int n, fint *NE, fint *NNE)
#endif
{
	real *lue;
	fint ne, nne;

	ne = nne = 0;
	for(lue = lu + 2*n; lu < lue; lu += 2) {
		if (lu[0] > negInfinity) {
			if (lu[1] == lu[0]) {
				ne++;
				continue;
				}
			nne++;
			}
		if (lu[1] < Infinity)
			nne++;
		}
	*NE = ne;
	*NNE = nne;
	}

 static void
#ifdef KR_headers
bset(be, bne, lu, n, noff) fint *be, *bne; real *lu; int n, noff;
#else
bset(fint *be, fint *bne, real *lu, int n, int noff)
#endif
{
	fint i, k;

	i = noff;
	for(k = i + n; i < k; lu += 2, i++) {
		if (lu[0] > negInfinity) {
			if (lu[1] == lu[0]) {
				*be++ = i;
				continue;
				}
			*bne++ = i;
			}
		if (lu[1] < Infinity)
			*bne++ = -1 - i;
		}
	}

 static void
ncheck(fint Nvar)
{
	char buf[256];
	fint nx = nmax_();
	if (Nvar > nx) {
		sprintf(buf, "%s %ld,\n%s %ld.\n%s %s",
		 "DONLP2 was compiled with NX (max. number of variables) =",
		 (long)nx, "but N (current number of variables) =",
		 (long)Nvar,
		 "Adjust O8PARA.INC, recompile donlp2.f and setup1.f,",
		 "and relink donlp2.");
		solve_result_num = 500;
		write_sol(buf, 0, 0, &Oinfo);
		if (amplflag)
			return;
		exit(1);
		}
	}

 static void
nbcheck(fint Ng, fint Nh)
{
	char buf[600];
	int i;
	fint mx = mmax_();
	if (Ng + Nh > mx) {
		i = sprintf(buf, "%s %ld,\n%s %ld\n%s %ld.",
		 "DONLP2 was compiled with NRESM (max. number of constraints) =",
		 (long)mx, "but NG (current number of inequalities) =",
		 (long)Ng, "and NH (current number of equalities) =", (long)Nh);
		i += sprintf(buf+i,"\n%s\n%s\n%s %d constraints.",
			 "Bounds on variables contribute to these numbers;",
			 "range constraints contribute two inequalities.",
		"Before adjustments for bounds and ranges, there were",
			 n_con);
		sprintf(buf+i, "\n%s %s",
		 "Adjust O8PARA.INC, recompile donlp2.f and setup1.f,",
		 "and relink donlp2.");
		solve_result_num = 500;
		write_sol(buf, 0, 0, &Oinfo);
		if (amplflag)
			return;
		exit(1);
		}
	}

 static SufDecl
suftab[] = {
	{ "vsclog10", 0, ASL_Sufkind_var | ASL_Sufkind_real }
	};

#define Nvar	L[1]
#define Ng	L[2]
#define Nh	L[3]
#define Lineq	L[7]
#define Linge	L[8]
#define Nsbeq	L[9]
#define Nsbge	L[10]

 void
MAIN__(VOID)
{
	FILE *nl;
	fint nleq, nlge;

	asl = ASL_alloc(ASL_read_fg);
	if (!(stub = getstops(xargv, &Oinfo)))
		return;

	suf_declare(suftab, sizeof(suftab)/sizeof(SufDecl));
	nl = jac0dim(stub, (fint)strlen(stub));

	ncheck(Nvar = n_var);

	X0 = (real *)Malloc((n_var+n_con)*sizeof(real));
	if (n_con)
		pi0 = X0 + n_var;

	fg_read(nl,0);
	dense_j();

	if (--nprob >= 0 && nprob < n_obj && objtype[nprob]) {
		maximize = 1;
		objsign = -1;
		}
	if ((namelen = strlen(stub)) > 40)
		namelen = 40;

	nbcount(LUrhs, nlc, &nleq, &nlge);
	nbcount(LUrhs+2*nlc, n_con-nlc, &Lineq, &Linge);
	nbcount(LUv, n_var, &Nsbeq, &Nsbge);
	if (Nb = Lineq + Linge + Nsbeq + Nsbge + nleq + nlge) {
		Ng = Linge + Nsbge + nlge;
		Nh = Lineq + Nsbeq + nleq;
		nbcheck(Ng, Nh);
		B = (fint*)Malloc(Nb*sizeof(fint));
		Bne = B + Nh;
		bset(B, Bne, LUrhs, n_con, n_var);
		bset(B + Lineq + nleq, Bne + Linge + nlge, LUv, n_var, 0);
		Beq = B - 1;
		--Bne;
		Lineq += Nsbeq;
		Linge += Nsbge;
		}
	if (!L[4])
		L[4] = n_var;
	donlp2_();
	endinf_();
	}

 void
setup0_(VOID)
{
	setup1_(B, DA, L, stub, X0, namelen);
	}

 void
setup_(VOID)
{ setup2_(L+15, DA+2, L+15, DA+3); }

 void
#ifdef KR_headers
setup3_(cef, cegradf, ceg, cegradg, ceh, cegradh, xsc)
	fint *cef, *cegradf, *ceg, *cegradg, *ceh, *cegradh; real *xsc;
#else
setup3_(fint *cef, fint *cegradf, fint *ceg, fint *cegradg, fint *ceh, fint *cegradh, real *xsc)
#endif
{
	SufDesc *s;
	real *r, *re, rlast, rlastp10;
	c_ef = cef;
	c_egradf = cegradf;
	c_eg = ceg - 1;
	c_egradg = cegradg - 1;
	c_eh = ceh - 1;
	c_egradh = cegradh - 1;

	rlast = 0;
	if (s = suf_get("vsclog10", ASL_Sufkind_var|ASL_Sufkind_input)) {
		r = s->u.r;
		for(re = r + n_var; r < re; r++, xsc++)
			if (*r) {
				if (*r != rlast)
					rlastp10 = pow(10., rlast = *r);
				*xsc = rlastp10;
				}
		}
	}

 void
#ifdef KR_headers
ef_(X, F) real *X; real *F;
#else
ef_(real *X, real *F)
#endif
{
	static fint I0 = 0;
	fint ierror = 0;
	++*c_ef;
	*F = objsign * objval(nprob, X, &ierror);
	if (ierror)
		errset_(&I0);
	}

 void
#ifdef KR_headers
egradf_(X, G) real *X; real *G;
#else
egradf_(real *X, real *G)
#endif
{
	real *Ge;

	++*c_egradf;
	objgrd(nprob, X, G, (fint*)0);
	if (maximize)
		for(Ge = G + n_var; G < Ge; G++)
			*G = -*G;
	}

 static void
#ifdef KR_headers
IGchk(i, who) int i; char *who;
#else
IGchk(int i, char *who)
#endif
{
	if (i < 1 || i > Ng) {
		fprintf(Stderr, "%s(I = %d) has I outside [1, %ld]\n",
			who, i, Ng);
		exit(1);
		}
	}

 static void
#ifdef KR_headers
IHchk(i, who) int i; char *who;
#else
IHchk(int i, char *who)
#endif
{
	if (i < 1 || i > Nh) {
		fprintf(Stderr, "%s(I = %d) has I outside [1, %ld]\n",
			who, i, Nh);
		exit(1);
		}
	}

 void
#ifdef KR_headers
eg_(I, X, F) fint *I; real *X; real *F;
#else
eg_(fint *I, real *X, real *F)
#endif
{
	int i, J;
	fint ierror = 0, k;

	IGchk(i = (int)*I, "EG");
	c_eg[i]++;
	if ((J = (int)Bne[i]) < 0) {
		J = -1 - J;
		if (J < n_var)
			*F = LUv[2*J+1] - X[J];
		else {
			J -= n_var;
			*F = LUrhs[2*J+1] - conival(J, X, &ierror);
			if (ierror) {
 botch:
				k = J + Nh + 1;
				errset_(&k);
				}
			}
		}
	else {
		if (J < n_var)
			*F = X[J] - LUv[2*J];
		else {
			J -= n_var;
			*F = conival(J, X, &ierror) - LUrhs[2*J];
			if (ierror)
				goto botch;
			}
		}
	}

 void
#ifdef KR_headers
eh_(I, X, F) fint *I; real *X; real *F;
#else
eh_(fint *I, real *X, real *F)
#endif
{
	int i, J;
	fint ierror = 0, k;

	IHchk(i = (int)*I, "EH");
	c_eh[i]++;
	if ((J = (int)Beq[i]) < n_var)
		*F = X[J] - LUv[2*J];
	else {
		J -= n_var;
		*F = conival(J, X, &ierror) - LUrhs[2*J];
		if (ierror) {
			k = J + 1;
			errset_(&k);
			}
		}
	}

 void
#ifdef KR_headers
egradg_(I, X, G) fint *I; real *X; real *G;
#else
egradg_(fint *I, real *X, real *G)
#endif
{
	int i, J;
	real *Ge;

	IGchk(i = (int)*I, "EGRADG");
	c_egradg[i]++;
	if ((J = (int)Bne[i]) < 0) {
		J = -1 - J;
		if (J < n_var) {
			memset(G, 0, n_var*sizeof(real));
			G[J] = -1.;
			}
		else {
			J -= n_var;
			congrd(J, X, G, (fint*)0);
			for(Ge = G + n_var; G < Ge; G++)
				*G = -*G;
			}
		}
	else {
		if (J < n_var) {
			memset(G, 0, n_var*sizeof(real));
			G[J] = 1.;
			}
		else {
			J -= n_var;
			congrd(J, X, G, (fint*)0);
			}
		}
	}

 void
#ifdef KR_headers
egradh_(I, X, G) fint *I; real *X; real *G;
#else
egradh_(fint *I, real *X, real *G)
#endif
{
	int i, J;

	IHchk(i = (int)*I, "EGRADH");
	c_egradh[i]++;
	if ((J = (int)Beq[i]) < n_var) {
		memset(G, 0, n_var*sizeof(real));
		G[J] = 1.;
		}
	else {
		J -= n_var;
		congrd(J, X, G, (fint*)0);
		}
	}

 static fint
#ifdef KR_headers
nev(n, c) fint n; fint *c;
#else
nev(fint n, fint *c)
#endif
{
	/* sum c[1] ... c[n] */
	fint *ce = c + n;
	n = 0;
	while(c < ce)
		n += *++c;
	return n;
	}

 void
#ifdef KR_headers
outinf_(ni, di, x, u) fint *ni; real *di; real *x; real *u;
#else
outinf_(fint *ni, real *di, real *x, real *u)
#endif
{
	int i, n, nint;
	fint *b, *be, k;
	real t, t1, tneg, *y;
	typedef struct { char *msg; int code; } Sol_info;
	static Sol_info solinfo[19] = {
	  { "unknown termination reason", 501 },
	  { "cannot evaluate constraints", 510 },
	  { "cannot evaluate objective", 511 },
	  { "dual extended QP failure: singular working set", 520 },
	  { "infeasible QP", 521 },
	  { "no descent in QP", 522 },
	  { "tiny QP step from infeasible point", 523 },
	  { "nondescent direction from QP", 524 },
	  { "reached maxit steps", 400 },
	  { "no acceptable step size", 525 },
	  { "tiny correction from QP at infeasible point", 200 },
	  { "Success! KKT conditions satisfied", 0 },
	  { "tiny step", 527 },
	  { "nearly feasible with small directional deriv.", 100 },
	  { "relaxed KKT conditions satisfied: singular point", 101 },
	  { "slow primal progress: singular or ill-conditioned problem?", 540 },
	  { "more than NRESET small primal corrections", 528 },
	  { "small QP correction: nearly feasible and singular", 101 },
	  { "NUMSM (default max(N,10)) small penalty function differences", 529 }
	  };
	char buf[512];

	k = ni[0];
	if (k < 0 || k > 18)
		k = 0;
	solve_result_num = solinfo[k].code;
	i = sprintf(buf, "DONLP2: %s\nF = %.*g", solinfo[k].msg, obj_prec(),
		objsign*di[0]);
	i += sprintf(buf+i, "\n%ld iters; %ld function, %ld gradient evals\n",
		ni[1], *c_ef, *c_egradf);
	i += sprintf(buf+i, "%ld component constraint, %ld grad evals\n",
		nev(Ng,c_eg) + nev(Nh,c_eh),
		nev(Ng,c_egradg) + nev(Nh,c_egradh));
	i += sprintf(buf+i, "final objective scaling = %g\n", di[1]);
	i += sprintf(buf+i, "norm(grad(f)) = %g, Lagrangian violation = %g\n",
		di[2], di[3]);
	i += sprintf(buf+i, "feas. violation = %g, dual feas. violation = %g",
		di[4], di[5]);
	if (nint = nlogv + niv + nlvbi + nlvci + nlvoi)
		sprintf(buf+i, "\nIntegrality of %d variables ignored.", nint);
	if (y = pi0) {
		tneg = -objsign;
		b = B;
		be = b + Nh;
		n = n_var;
		while(b < be) {
			if ((i = *b++) >= n)
				y[i-n] = *u;
			u++;
			}
		be = B + Nb;
		while(b < be) {
			if ((i = *b++) >= 0) {
				if (b < be && *b == -1 - i) {
					if (i >= n) {
						t = u[0];
						t1 = u[1];
						if (t < -t1)
							t = -t1;
						y[i-n] = t*objsign;
						}
					b++;
					u += 2;
					continue;
					}
				if (i >= n)
					y[i-n] = objsign**u;
				}
			else {
				i = -1 - i;
				if (i >= n)
					y[i-n] = tneg**u;
				}
			u++;
			}
		}
	write_sol(buf, x, y, &Oinfo);
	}

 void
solchk_(VOID) {}
