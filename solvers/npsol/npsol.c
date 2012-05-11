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

#include "nlp.h"
#include "getstub.h"
#include "signal.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef int (*c_fp)ANSI((fint*,fint*,fint*,fint*,fint*,real*,real*,real*,fint*));
typedef int (*o_fp)ANSI((fint*,fint*,real*,real*,real*,fint*));

static real *gs, objsign;
static fint NERROR = -1;
static int nobj = -12345;
#define asl cur_ASL

extern int npsol_ ANSI((fint *n, fint *nclin, fint *ncnln, fint *nrowa,
	fint *nrowuj, fint *nrowr, real *a, real *bl, real *bu,
	c_fp confun, o_fp objfun, fint *inform, fint *iter,
	fint *istate, real *c, real *ujac, real *clamda, real *objf,
	real *ugrad, real *r, real *x, fint *iw, fint *leniw,
	real *w, fint *lenw));

extern int npoptn_ ANSI((char*,fint));
extern char npsol_version[];

#ifdef __cplusplus
	}
#endif

 static int
#ifdef KR_headers
objfun(MODE, N, X, F, G, NSTATE) fint *MODE, *N, *NSTATE; real *X, *F, *G;
#else
objfun(fint *MODE, fint *N, real *X, real *F, real *G, fint *NSTATE)
#endif
{
	real *ge;

	Not_Used(NSTATE); /* quiet non-use warning */
	if (*MODE != 1)
		*F = objsign * objval(nobj, X, &NERROR);
	if (*MODE > 0) {
		objgrd(nobj, X, G, &NERROR);
		if (objsign < 0)
			for(ge = G + *N; G < ge; G++)
				*G = -*G;
		}
	return 0;
	}

 static int
#ifdef KR_headers
obj0fun(MODE, N, X, F, G, NSTATE) fint *MODE, *N, *NSTATE; real *X, *F, *G;
#else
obj0fun(fint *MODE, fint *N, real *X, real *F, real *G, fint *NSTATE)
#endif
{
	/* variant for no objective */
	real *ge;

	if (*MODE != 1)
		*F = 0;
	if (*MODE > 0)
		for(ge = G + *N; G < ge; G++)
			*G = 0;
	return 0;
	}

 static int
#ifdef KR_headers
confun(MODE, NCNLN, N, NROWJ, NEEDC, X, C, CJAC, NSTATE)
	fint *MODE, *NCNLN, *N, *NROWJ, *NEEDC, *NSTATE; real *X, *C, *CJAC;
#else
confun(fint *MODE, fint *NCNLN, fint *N, fint *NROWJ, fint *NEEDC,
	real *X, real *C, real *CJAC, fint *NSTATE)
#endif
{
	int i, n, nl, nrj;
	real *cj, *cje, *g;

	NSTATE = NSTATE; NCNLN = NCNLN;
	n = (int)*N;
	nl = nlc;
	if (*MODE != 1)
		for(i = 0; i < nl; i++)
			if (NEEDC[i])
				C[i] = conival(i, X, &NERROR);
	if (*MODE > 0) {
		nrj = *NROWJ;
		cje = CJAC + n*nrj;
		for(i = 0; i < nl; i++)
			if (NEEDC[i]) {
				congrd(i, X, g = gs, &NERROR);
				for(cj = CJAC+i; cj < cje; cj += nrj)
					*cj = *g++;
				}
		}
	return 0;
	}

 static char *
#ifdef KR_headers
Dval(oi, kw, v) Option_Info *oi; keyword *kw; char *v;
#else
Dval(Option_Info *oi, keyword *kw, char *v)
#endif
{
	char buf[128], *rv, *s, *se;

	Not_Used(oi);

	strtod(v,&rv);
	if (rv > v) {
		se = buf + sizeof(buf) - 1;
		strcpy(buf, kw->desc);
		s = buf + strlen(buf);
		strcpy(s, " = ");
		s += 3;
		while(v < rv && s < se)
			*s++ = *v++;
		*s = 0;
		npoptn_(buf, (fint)(s - buf));
		}
	return rv;
	}

 static char *
#ifdef KR_headers
Kval(oi, kw, v) Option_Info *oi; keyword *kw; char *v;
#else
Kval(Option_Info *oi, keyword *kw, char *v)
#endif
{
	Not_Used(oi);
	npoptn_(kw->desc, (fint)strlen(kw->desc));
	return v;
	}

 static keyword keywds[] = {
	KW("cdiff", Dval, 0, "Central Difference Interval"),
	KW("cend", Dval, 0, "Stop Constraint Check"),
	KW("crashtol", Dval, 0, "Crash Tolerance"),
	KW("cstart", Dval, 0, "Start Constraint Check"),
	KW("derlev", Dval, 0, "Derivative Level"),
	KW("diff", Dval, 0, "Difference Interval"),
	KW("feastol", Dval, 0, "Feasibility Tolerance"),
	KW("fprec", Dval, 0, "Function Precision"),
	KW("infbound", Dval, 0, "Infinite Bound Size"),
	KW("infstep", Dval, 0, "Infinite Step Size"),
	KW("itlim", Dval, 0, "Major Iterations Limit"),
	KW("lftol", Dval, 0, "Linear Feasibility Tolerance"),
	KW("lstol", Dval, 0, "Linesearch Tolerance"),
	KW("maxit", Dval, 0, "Major Iterations Limit"),
	KW("mitlim", Dval, 0, "Minor Iterations Limit"),
	KW("mprlev", Dval, 0, "Minor Print Level"),
	KW("nftol", Dval, 0, "Nonlinear Feasibility Tolerance"),
	KW("objno", I_val, &nobj, "Ojbective number: 1 = first objective"),
	KW("oend", Dval, 0, "Stop Objective Check"),
	KW("ostart", Dval, 0, "Start Objective Check"),
	KW("prlev", Dval, 0, "Major Print Level"),
	KW("tol", Dval, 0, "Optimality Tolerance"),
	KW("vc", Kval, 0, "Verify Constraint Gradients"),
	KW("vco", Kval, 0, "Verify Gradients"),
	KW("verify", Kval, 0, "Verify Gradients"),
	KW("vo", Kval, 0, "Verify Objective Gradients"),
	KW("wantsol", WS_val, 0, WS_desc_ASL+5)
	};

 static fint
#ifdef KR_headers
npkey(s, slen) char *s; fint slen;
#else
npkey(char *s, fint slen)
#endif
{
	npoptn_(s, slen);
	return 0;	/* Assume it's OK: there's no way to tell... */
	}

 static Option_Info Oinfo = { "npsol", "NPSOL", "npsol_options",
				keywds, nkeywds, 1, npsol_version, 0, npkey };

static Jmp_buf Jb;

#ifndef Sig_ret_type
#define Sig_ret_type void
#endif
 static Sig_ret_type
#ifdef KR_headers
catchfpe(n) int n;
#else
catchfpe(int n)
#endif
{
	Not_Used(n);
	report_where(asl);
	printf("\nFloating point error.\n");
	fflush(stdout);
	need_nl = 0;
	longjmp(Jb.jb,2);
	}

#undef asl

 void
MAIN__(VOID)
{
	fint M, N, N1, NO, MXROW, MXCOL, NZ;
	fint INFORM, ITER, L, LENIW, LENW, NCLIN, NCNLN, NROWA, NROWJ, NROWR;
	fint *ISTATE, *IW;
	real OBJF;
	real *A, *BL, *BU, *C, *CJAC, *CLAMDA, *GRAD, *R, *W;
	real *a, *bl, *bu, *lu, *lue;
	int i, nv;
	extern char **xargv;
	cgrad *gr, **gr0, **gr1;
	char *b, *stub;
	char buf[100];
	ASL *asl;
	typedef struct { char *msg; int code, wantsol; } Sol_info;
	Sol_info *SI;
	static Sol_info solinfo[] = {
	  { /* 0 */ "optimal solution", 0, 1 },
	  { /* 1 */ "first-order optimal but not converged", 100, 1 },
	  { /* 2 */ "linear constraints infeasible", 200, 0 },
	  { /* 3 */ "nonlinear constraints infeasible", 201, 0 },
	  { /* 4 */ "major iteration limit", 400, 1 },
	  { /* 5 */ "Bug: INFORM = 5", 501, 0 },
	  { /* 6 */ "stuck (INFORM = 6)", 502, 1 },
	  { /* 7 */ "Bug: bad objective or constraints (INFORM = 7)", 503, 0 },
	  { /* 8 */ "Bug: INFORM = 8", 504, 0 },
	  { /* 9 */ "Bug: invalid input parameter (INFORM = 9)", 505, 0 },
	  { /* 10 */ "unexpected INFORM value", 506, 0 },
	  { /* 11 */ "solution aborted", 507, 0 }
	  };

	asl = ASL_alloc(ASL_read_fg);
	stub = getstops(xargv, &Oinfo);
	if (jac1dim(stub, &M, &N, &NO, &NZ, &MXROW, &MXCOL, (fint)strlen(stub)))
		exit(4);
	if (nobj == -12345)
		nobj = n_obj > 0;
	else if (nobj < 0 || nobj > n_obj) {
		fprintf(Stderr,
		"Bad objno = %d; must be between 0 (for no objective) and %d\n",
			nobj, n_obj);
		exit(1);
		}
	--nobj;
	dense_j();
	npoptn_("NOLIST", 6L);
	npoptn_("VE L -1", 7L);
	npoptn_("MA PR 0", 7L);
	objsign = nobj >= 0 && objtype[nobj] ? -1. : 1.;
	nv = (int)N;
	NCNLN = nlc;
	NCLIN = M - NCNLN;
	NROWA = NCLIN > 0 ? NCLIN : 1;
	NROWJ = NCNLN > 0 ? NCNLN : 1;
	NROWR = N;
	LENIW = 3*N + NCLIN + 2*NCNLN;
	N1 = N + NCLIN + NCNLN;

	if (!NCNLN)
		LENW = NCLIN ? (2*N + 20)*N + 11*NCLIN
			     : 20*N;
	else
		LENW = (2*N + NCLIN + 2*NCNLN + 20)*N + 11*NCLIN + 21*NCNLN;
	L = (LENW + (NROWA+NROWJ+NROWR+1)*N + 3*N1 + NCNLN)*sizeof(real)
			+ (LENIW + N1)*sizeof(fint);
	A = (real *)M1alloc(L);
	bl = BL = A + NROWA*N;
	bu = BU = BL + N1;
	C = BU + N1;
	CJAC = C + NCNLN;
	CLAMDA = CJAC + NROWJ*N;
	GRAD = CLAMDA + N1;
	R = GRAD + N;
	W = R + NROWR*N;
	IW = (fint *)(W + LENW);
	ISTATE = IW + LENIW;
	if (NCLIN) {
		memset((char *)A, 0, NROWA*N*sizeof(real));
		gr0 = Cgrad + nlc;
		gr1 = gr0 + NCLIN;
		for(a = A; gr0 < gr1; a++)
			for(gr = *gr0++; gr; gr = gr->next)
				a[NROWA*gr->varno] = gr->coef;
		}
	if (NCNLN)
		memset((char *)CJAC, 0, NROWJ*N*sizeof(real));

	lu = LUv;
	lue = lu + 2*N;
	while(lu < lue) {
		*bl++ = *lu++;
		*bu++ = *lu++;
		}

	lu = LUrhs + 2*nlc;
	lue = LUrhs + 2*M;
	while(lu < lue) {
		*bl++ = *lu++;
		*bu++ = *lu++;
		}

	lu = LUrhs;
	lue = LUrhs + 2*nlc;
	while(lu < lue) {
		*bl++ = *lu++;
		*bu++ = *lu++;
		}
	gs = (real *)M1alloc(N*sizeof(real));

	want_deriv = 1;

	err_jmp1 = &Jb;
	if (setjmp(Jb.jb))
		i = 11;
	else {
		signal(SIGFPE, catchfpe);
		errno = 0;
		npsol_(&N, &NCLIN, &NCNLN, &NROWA, &NROWJ, &NROWR,
			A, BL, BU,
			confun, nobj >= 0 ? objfun : obj0fun,
			&INFORM, &ITER, ISTATE,
			C, CJAC, CLAMDA, &OBJF, GRAD, R, X0,
			IW, &LENIW, W, &LENW);
		if ((i = (int)INFORM) < 0 || i > 10)
			i = 10;
		}
	SI = solinfo + i;
	solve_result_num = SI->code;
	b = buf + Sprintf(buf, "NPSOL: %s.\n%ld iterations", SI->msg, ITER);
	if (SI->wantsol)
		Sprintf(b, ", objective %.*g", obj_prec(), objsign*OBJF);

	/* Adjust variables */

	if ((i = n_con - nlc) && nlc) {
		memcpy((char *)bl, (char *)(CLAMDA+nv+i),
			nlc*sizeof(real));
		memcpy((char *)(bl + nlc), (char *)(CLAMDA+nv),
			i*sizeof(real));
		}
	else
		bl = CLAMDA + nv;

	if (objsign < 0)
		for(bu = bl + M; bu-- > bl; )
			*bu = -*bu;

	write_sol(buf, X0, bl, &Oinfo);
	}
