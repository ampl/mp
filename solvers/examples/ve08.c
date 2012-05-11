/****************************************************************
Copyright (C) 1997-1998 Lucent Technologies
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
#ifdef PSEDREAD
#include "asl_pfg.h"
#undef fg_read
#define fg_read(a,b) pfg_read_ASL((ASL*)asl,a,b)
#define psinfo asl->P
 static ASL_pfg *asl;
#else
 static ASL *asl;
#endif
real macheps =
#ifdef IEEE_MC68k
	2.220446049250313e-16
#else
#ifdef IEEE_8087
	2.220446049250313e-16
#else
#ifdef VAX
	2.77555756156289135e-17
#else
	MACHEPS
#endif
#endif
#endif
	;

#ifdef __cplusplus
extern "C" {
#endif

 extern void f_exit(void);

#ifdef Timing
int nTimes;
static real ftime;
#endif

typedef fint logical;
typedef real (*D_fp)(fint*);

typedef void (*range_fp)(fint *k, fint *mode, real *w1,
	real *w2, fint *ndimk, fint *nsubk, fint *ns);

typedef void (*elfnct_fp)(fint *k, real *x, real *fx,
	real *gx, fint *ndimk, fint *ns, fint *jflag,
	real *fmax, real *fnoise);

extern void ve08ad_(fint *n, fint *ns, fint *invar, fint *nvar,
	elfnct_fp Elfnct, range_fp Range, D_fp xlower, D_fp xupper,
	real *x, real *fx, real *epsil, real *ebound,
	fint *ngr, fint *nit, logical *fknown, real *flowbd,
	real *relpr, real *difgrd, logical *restrt,
	logical *testgx, logical *hesdif, real *stmax,
	real *stepl, fint *istate, fint *ipdevc,
	fint *ipfreq, fint *ipwhat, fint *lwk, real *wk,
	fint *info, fint *iflag);

 /* PORT 3 routines */

 extern void dl7itv_(fint *N, real *X, real *L, real *Y);
 extern void dl7ivm_(fint *N, real *X, real *L, real *Y);
 extern void dl7tvm_(fint *N, real *X, real *L, real *Y);
 extern void dl7vml_(fint *N, real *X, real *L, real *Y);
 extern void dq7rgs_(fint *ierr, fint *ipivot, fint *L, fint *N, fint *NN,
	fint *NOPIVK, fint *P, real *Q, real *R, real *W);

#ifdef __cplusplus
	}
#endif

fint hesdif, ipwhat, maxfun[2], maxit[2], outlev, showstats = 1, solprt,
	testgx, x0prt;
int nprob = 1;
static long ranges = 1, useRP;
real difgrd, ebound, gtol, stepl[2], stmax;

#ifdef PSEDREAD
static fint I0 = 0, I1 = 1;
#endif
static int heslen, nfcall;

static keyword keywds[] = { /* must be alphabetical */
 KW("difgrd", D_val, &difgrd, "rel. step size for finite diffs for first gradient"),
 KW("ebound", D_val, &ebound, "tolerance for bounds considered active"),
 KW("gtol", D_val, &gtol, "stopping tolerance for projected gradient norm"),
 KW("hesdif", L_val, &hesdif, "compute initial Hessians by finite diffs?"),
#ifdef PSEDREAD
 KW("ignore", IA_val, voffset_of(ASL_pfg,i.xknown_ignore), "ignore xknown calls"),
#endif
 KW("ipwhat", L_val, &ipwhat, "what to print"),
 KW("maxfun", L_val, &maxfun[0], "maximum function evaluations"),
 KW("maxfwd", IA_val, voffset_of(ASL,p.maxfwd_), "# of partials to forward recur; default = 5"),
 KW("maxit", L_val, &maxit[0], "maximum number of inner iterations per step"),
#ifdef PSEDREAD
 KW("merge", IA_val, voffset_of(ASL_pfg,P.merge), "do not merge elements with the same internal variables"),
#endif
 KW("nprob", I_val, &nprob, "objective choice: 1 (default) = 1st"),
 KW("outlev", L_val, &outlev, "how often to print"),
 KW("ranges", L_val, &ranges, "use ranges"),
 KW("showstats", L_val, &showstats, "show timing and memory statistics"),
 KW("solprt", L_val, &solprt, "print final solution"),
 KW("stepl", D_val, &stepl[0], "maximum step length for first step"),
 KW("stmax", D_val, &stmax, "maximum step length"),
 KW("testgx", L_val, &testgx, "test accuracy of gradients?"),
 KW("use_rp", L_val, &useRP, "use R and P in ranges"),
 KW("wantsol", WS_val, 0, WSu_desc_ASL+5),
 KW("x0prt", L_val, &x0prt, "print initial x")
 };

 static Option_Info Oinfo = {
#ifdef PSEDREAD
	"v8", "V8",
#else
	"ve08", "VE08",
#endif
	"ve08_options", keywds, nkeywds, 1 };

 static void
Elfnct(fint *kp, real *X, real *F,
	real *G, fint *ndimk, fint *ns, fint *jflag,
	real *fmax, real *fnoise)
{
	fint nerror = 0;
#ifdef Timing
	real t0 = xectim_();
#endif

#ifdef PSEDREAD /*{{*/
	ps_func *f = psinfo.ops + nprob;
	fint k = *kp - 1;
	psb_elem *b;
	int i, deriv, nintv;
	int *v, *v0, *ve;
	cexp *c;
	expr *e;
	real t;
	linpart *L, *Le;
	cde *d;
	ograd *og;
	range *r;
	linarg *la, **lap, **lap0, **lape;
	register real *Adjoints;

	Not_Used(ns);
	Not_Used(fmax);
	Not_Used(ndimk);

	errno = 0;
	want_deriv = deriv = *jflag == 2;
	if (k == f->nb) {
		og = Ograd[nprob];
		t = og->coef * *X;
		if (deriv) {
			*G = og->coef;
			while(og = og->next) {
				t += og->coef * *++X;
				*++G = og->coef;
				}
			}
		else
			while(og = og->next)
				t += og->coef * *++X;
		*F = t;
		}
	else {
		b = f->b + k;
		r = b->U;
		v = v0 = r->ui;
		ve = v + r->nv;
		lap0 = 0;
		if (nintv = r->nintv) {
			lap = lap0 = r->lap;
			lape = lap + nintv;
			}
		if (asl->i.x_known)
			goto known_x;
		while(v < ve)
			var_e[*v++].v = *X++;
		if (lap = lap0) {
			do {
				la = *lap++;
				og = la->nz;
				t = var_e[og->varno].v*og->coef;
				while(og = og->next)
					t += var_e[og->varno].v*og->coef;
				la->v->v = t;
				}
				while(lap < lape);
			}
		if (v = b->ce) {
			ve = v + *v;
			do {
				cv_index = i = *++v;
				c = cexps + i;
				e = c->e;
				t = (*e->op)(e C_ASL);
				if (L = c->L)
					for(Le = L + c->nlin; L < Le; L++)
						t += L->fac *
							((expr_v*)L->v.vp)->v;
				psinfo.vp[i]->v = t;
				if (c->funneled && deriv)
					funnelset(c->funneled);
				}
				while(v < ve);
			cv_index = 0;
			}
 known_x:
		d = &b->D;
		e = d->e;
		*F = (*e->op)(e C_ASL);
		if (deriv) {
			Adjoints = adjoints;
			if (lap = lap0)
				do Adjoints[(*lap++)->v->a] = 0;
					while(lap < lape);
			v = v0 = r->ui;
			ve = v + r->nv;
			while(v < ve)
				Adjoints[*v++] = 0;
			if (i = d->zaplen) {
				memset((char *)adjoints_nv1, 0, i);
				derprop(d->d);
				}
			if (lap = lap0)
				do  {
					t = Adjoints[(la = *lap++)->v->a];
					og = la->nz;
					do Adjoints[og->varno] += t*og->coef;
						while(og = og->next);
					}
					while(lap < lape);
			while(v0 < ve)
				*G++ = Adjoints[*v0++];
			}
		}
#else /*}{ PSEDREAD */
	Not_Used(kp);
	Not_Used(ns);
	Not_Used(fmax);
	Not_Used(ndimk);

	*F = objval(nprob, X, &nerror);
	if (!nerror && *jflag == 2)
		objgrd(nprob, X, G, &nerror);
#endif /*}} PSEDREAD */
	*fnoise = 0;
	if (nerror)
		*jflag = -nerror;
#ifdef Timing
	t0 -= xectim_();
	ftime -= t0;
#endif
	nfcall++;
	}

 static void
show_x(char *when, long n, int nb, real *X, real *LU)
{
	long i;

	if (nb) {
		printf("i\t%s x(i)\tL(i)\t\tU(i)\n", when);
		for(i = 1; i <= n; i++, LU += 2)
			printf("%d\t%9g\t%9g\t%9g\n",
				i, *X++, LU[0], LU[1]);
		}
	else {
		printf("i\t%s x(i)\n", when);
		for(i = 1; i <= n; i++)
			printf("%d\t%g\n", i, *X++);
		}
	}

 static real
xlower(fint *I)
{
	return LUv[(*I-1)<<1];
	}

 static real
xupper(fint *I)
{
	return LUv[(*I<<1) - 1];
	}

#ifdef PSEDREAD
 static real **Qmats, *Xscratch;
 static fint *Qrank;

 static void
Range_setup(ps_func *f, int N, int nog, fint *istate)
{
	range *r;
	psb_elem *b, *be;
	int *vr, *vre;
	fint i, m, n, ni, nr, ns, nx;
	fint *Nq, *Nqe, *Iscratch;
	real *q, **qp, *rmat, *w, *z;
	linarg *la, **lap, **lape;
	ograd *og;

	nr = nog + 1;
	nx = nog;
	ns = f->nb;
	b = f->b;
	be = b + ns++;
	Nq = Qrank = (fint*)Malloc(ns*sizeof(fint));
	Nqe = Nq + ns;
	do *Nq++ = -1;
		while(Nq < Nqe);
	if (!ranges)
		return;
	ni = ns;
	Nq = Qrank;
 next_R:
	while(b < be) {
		Nq++;
		if ((r = (b++)->U) && (n = r->n) > 0 && n < r->nv) {
			m = r->nv;
			if (istate) {
				/* No range mapping for fixed elements. */
				/* Only possible with presolve disabled. */
				vr = r->ui;
				vre = vr + m;
				while(vr < vre)
					if (!istate[*vr++])
						goto next_R;
				}
			Nq[-1] = 0;
			if (nx < n)
				nx = n;
			ni += n;
			nr += n*m + (n*(n+1) >> 1);
			}
		}
	w = Xscratch = (real*)Malloc((nr+nx)*sizeof(real));
	z = w + nx;
	Qmats = qp = (real **)Malloc(ns*sizeof(real*));
	Iscratch = (fint *)Malloc(N*sizeof(fint));
	b = f->b;
	Nqe = Qrank;
	while(b < be) {
		r = (b++)->U;
		*qp++ = z;
		if (*(Nq = Nqe++) || !r)
			continue;
		n = r->nv;
		vr = r->ui;
		vre = vr + n;
		i = 0;
		do Iscratch[*vr++] = i++;
			while(vr < vre);
		m = r->n;
		lap = r->lap;
		lape = lap + m;
		q = z;
		memset(q, 0, m*n*sizeof(real));
		do {
			la = *lap++;
			for(og = la->nz; og; og = og->next)
				z[Iscratch[og->varno]] = og->coef;
			z += n;
			} while(lap < lape);
		rmat = z;
		z += m*(m+1) >> 1;
		*Iscratch = 1;	/* dq7rgs need not initialize Iscratch */
		dq7rgs_(Nq, Iscratch, &I1, &n, &n, &I0, &m, q, rmat, w);
		if (*Nq)
			--*Nq;
		else
			*Nq = m;
		}
	*qp = q = z;
	if ((n = nog) > 0) {
		for(og = Ograd[nprob]; og; og = og->next)
			*z++ = og->coef;
		*Iscratch = 1;
		dq7rgs_(Nqe, Iscratch, &I1, &n, &n, &I0, &I1, q, z, w);
		}
	}

 static void
matmul(fint m, fint n, real *y, real *A, real *x)	/* y := A*x */
{
	real *A1, t, *x1, *xe, *ye;

	xe = x + n;
	ye = y + m;
	while(y < ye) {
		x1 = x;
		A1 = A++;
		t = *A1 * *x1++;
		while(x1 < xe) {
			A1 += m;
			t += *A1 * *x1++;
			}
		*y++ = t;
		}
	}

 static void
mattmul(fint m, fint n, real *y, real *A, real *x)	/* y := A'*x */
{
	real t, *x1, *xe, *ye;

	xe = x + m;
	ye = y + n;
	while(y < ye) {
		x1 = x;
		t = *A++ * *x1++;
		while(x1 < xe)
			t += *A++ * *x1++;
		*y++ = t;
		}
	}
#endif /* PSEDREAD */

 static void
Range(fint *kp, fint *mode, real *w1, real *w2, fint *ndimk,
	fint *nsubk, fint *ns)
{
	fint m = *ndimk;
#ifdef PSEDREAD /*{{ */
	fint k = *kp - 1;
	fint n;
	real *Q, *R, *w;

	if ((n = Qrank[k]) > 0) {
		*nsubk = n;
		Q = Qmats[k];
		R = Q + m*n;
		switch(*mode) {
		  case 1:
			mattmul(m, n, w2, Q, w1);
			if (useRP)
				dl7vml_(&n, w2, R, w2);
			return;
		  case 2:
			w = w1;
			if (useRP)
				dl7ivm_(&n, w = Xscratch, R, w1);
			matmul(m, n, w2, Q, w);
			return;
		  case 3:
			w = w1;
			if (useRP)
				dl7tvm_(&n, w = Xscratch, R, w1);
			matmul(m, n, w2, Q, w);
			return;
		  case 4:
			mattmul(m, n, w2, Q, w1);
			if (useRP)
				dl7itv_(&n, w2, R, w2);
			return;
		  }
		}
#else  /*ifndef PSEDREAD }{ */
	Not_Used(kp);
	Not_Used(mode);
#endif /*PSEDREAD }}*/
	Not_Used(ns);
	*nsubk = m;
	while(--m >= 0)
		*w2++ = *w1++;
	}

#ifdef PSEDREAD /*{*/
 static void
funnelzap(funnel *f)
{
	funnel *f1;

	while(f1 = f) {
		f = f1->next;
		f1->next = 0;
		}
	}

 static int
compar(const void *a, const void *b)
{ return *(int*)a - *(int*)b; }
#endif /*} PSEDREAD */

 void
MAIN__(void)
{
	fint LW[2], M, N, NE, NO, MXROW, MXCOL, NQ, NS, NV, NZ;
	fint iflag, info, nb;
	fint *istate;
	static fint fknown, ipdevc = 6, restrt;
	real *B, *W, *X;
	real Inf, f, negInf;
	int i, j;
	extern int xargc;
	extern char **xargv;
	char **av = xargv, msg[160], *stub;
	fint *nvar, *invar;
	FILE *nl;
	typedef struct { char *msg; int code, iflag; } Sol_info;
	Sol_info *SI;
	static Sol_info solinfo[] = {
		{ "Solved: norm2(grad) < gtol", 0, 1 },
		{ "Too many function evaluations", 400, 2 },
		{ "Tiny step and function improvement", 100, 3 },
		{ "Iteration limit", 401, 18 },
		{ "Unbounded", 300, 26 },
		{ "Linesearch step got too small", 101, 28 },
		{ "Tiny function change", 102, 29 },
		{ "Tiny search direction", 103, 35},
		{ "Bug", 500, 999 }
		};
#ifdef PSEDREAD /*{*/
	ps_func *pf;
	psb_elem *b, *be;
	ograd *og;
	int nog, nog0, *vi, *vie;
	fint *invi, *nvi;
	range *r;
	fint *istchk = 0;
#endif /*} PSEDREAD */
#ifdef Timing
	real tend, tsetup, tstart;
	char *sbrk0;
	extern char *sbrk(int);
	long edagmem, totmem;

	showstats = 1;
	tstart = xectim_();
	sbrk0 = sbrk(0);
#endif

#ifdef PSEDREAD
	asl = (ASL_pfg *)ASL_alloc(ASL_read_pfg);
#else
	asl = ASL_alloc(ASL_read_fg);
#endif
	stub = getstub(&av, &Oinfo);
	nl = jacdim(stub, &M, &N, &NO, &NZ, &MXROW, &MXCOL,
		(long)strlen(stub));
	if (!nl)
		return;
	if (M) {
		fprintf(Stderr, "%s: ignoring %ld constraints\n",
			progname, M);
		fflush(Stderr);
		}
	X0 = (real*)Malloc(N*sizeof(real));
	fg_read(nl,ASL_GJ_zerodrop|ASL_find_default_no_groups);
#ifdef Timing
	edagmem = sbrk(0) - sbrk0;
#endif
	gtol = sqrt(macheps);
	stepl[0] = stmax = -1;
	maxit[0] = N >> 1;
	if (maxit[0] < 50)
		maxit[0] = 50;
	maxfun[0] = 150 * N;
	if (amplflag) {
		showstats = 0;
		outlev = -3;
		}
	if (getopts(av, &Oinfo))
		exit(1);
	if (nprob <= 0 || nprob > n_obj) {
		fprintf(Stderr, "nprob = %d must lie in [1, %d].\n",
			nprob, n_obj);
		exit(1);
		}
	--nprob;
#ifdef PSEDREAD /*{{*/
	funnelzap(f_b);
	funnelzap(f_o);
	pf = psinfo.ops + nprob;
	NS = pf->nb;
	nog = 0;
	for(og = Ograd[nprob]; og; og = og->next)
		nog++;
	b = pf->b;
	be = b + NS;
	heslen = nog*(nog+1) >> 1;
	if (nog0 = nog)
		NS++;
	while(b < be)
		nog += (b++)->U->nv;
	nvi = nvar = (fint *)Malloc((nog+NS+1)*sizeof(fint));
	invi = invar = nvar + (NS+1);
	nvar[0] = i = 1;
	b = pf->b;
	while(b < be) {
		r = (b++)->U;
		vi = r->ui;
		*++nvi = i += j = r->nv;
		if (j > 1)
			qsort(vi, j, sizeof(fint), compar);/* !!for debugging*/
		heslen += j*(j+1) >> 1;
		vie = vi + j;
		do *invi++ = *vi++ + 1;
			while(vi < vie);
		}
	if (nog0) {
		*++nvi = i + nog0;
		for(og = Ograd[nprob]; og; og = og->next)
			*invi++ = og->varno + 1;
		}
#else /*}{ ifndef PSEDREAD */
	NS = 1;
	nvar = (fint *)Malloc((N+2)*sizeof(fint));
	invar = nvar + 2;
	for(i = 0; i < N; i++)
		invar[i] = i + 1;
	nvar[0] = 1;
	nvar[1] = N + 1;
	heslen = N*(N+1)/2;
#endif /*}} PSEDREAD */

	NV = nvar[NS] - 1;
	NE = 0;
	for(i = 0; i < NS; i++) {
		j = nvar[i+1] - nvar[i];
		if (NE < j)
			NE = j;
		}
	NQ = N + N;
	if (NQ < NV)
		NQ = NV;

	LW[0] = NS + NV + 4*N + NQ + 3*NE + heslen;

	W = (real *)Malloc(LW[0]*sizeof(real) + (N+NS)*sizeof(fint));
	istate = (fint *)(W + LW[0]);

	X = X0;
	B = LUv;
	Inf = Infinity;

	negInf = -Inf;
	for(i = nb = 0; i < N; i++, B += 2) {
		if (B[0] == B[1]) {
			istate[i] = 0;
			X[i] = B[0];
#ifdef PSEDREAD
			istchk = istate;
#endif
			}
		else if (B[0] > negInf || B[1] < Inf) {
			nb++;
			istate[i] = 1;
			}
		else
			istate[i] = -1;
		}
	for(j = N + NS; i < j; i++)
		istate[i] = 1;
#ifdef PSEDREAD
	Range_setup(psinfo.ops + nprob, (int)N, nog0, istchk);
	if (showstats) {
		printf("%d initial elements\n%d distinct linear terms\n\
%d duplicate linear terms in different elements\n\
%d duplicate linear terms in the same element\n%d adjusted elements\n",
			psinfo.ns0, psinfo.nlttot, psinfo.ndupdt,
			psinfo.ndupst, NS);
		fflush(stdout);
		}
#endif
#ifdef Timing
	tsetup = xectim_() - tstart;
#endif
	if (x0prt)
		show_x("initial", N, nb, X, LUv);
	ve08ad_(&N, &NS, invar, nvar, Elfnct, Range, xlower, xupper,
		X, &f, &gtol, &ebound, maxfun, maxit, &fknown, &negInf,
		&macheps, &difgrd, &restrt, &testgx, &hesdif, &stmax,
		stepl, istate, &ipdevc, &outlev, &ipwhat, LW, W,
		&info, &iflag);

	i = (int) iflag;
	for(SI = solinfo; SI->iflag != i; SI++) {
		if (SI->iflag == 999) {
			SI->code += i;
			break;
			}
		}
	solve_result_num = SI->code;
	i = Sprintf(msg, "%s: %s.\nFinal f = ", Oinfo.sname, SI->msg);
	i += g_fmtop(msg+i, f);
	if (SI->iflag > 18)
		i += Sprintf(msg+i, " (iflag = %ld, info = %ld)",
			(long)iflag, (long)info);
	Sprintf(msg+i, "\n%ld evaluations, %ld iters.\n",
		(long)maxfun[1], (long)maxit[1]);
	write_sol(msg, X, 0, &Oinfo);
	if (solprt)
		show_x("final", N, nb, X, LUv);

#ifdef Timing
	if (showstats) {
	tend = xectim_();
	printf("setup time %g seconds\n", tsetup);
	printf("nfcall = %d (%g seconds, %g per call)\n", nfcall, ftime,
		nfcall ? ftime/nfcall : 0.);
	printf("Total time %g seconds\n", tend - tstart);
	}
	if (showstats) {
		totmem = sbrk(0) - sbrk0;
		printf(
	"malloced %ld bytes ( %ld for edag: %d derps, amax = %d )\n",
			totmem, edagmem, nderps, amax);
		}
#endif
	f_exit();
	}
