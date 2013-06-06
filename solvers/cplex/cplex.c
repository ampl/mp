/****************************************************************
Copyright (C) 1997-2001 Lucent Technologies
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

/* Variant that duplicates columns for ranges when presenting the
   dual problem to CPLEX.  (I thought this wasn't necessary, but
   CPLEX seems to botch the solve without it.)
 */

/* Compile with -DNO_CPLEX_MIP to disable MIP keywords. */

/* Compile with -DNO_BARRIER to disable barrier keywords. */

#ifndef NO_CPLEX_MIP
#undef CPLEX_MIP
#define CPLEX_MIP
#endif

#ifndef NO_BARRIER
#undef BARRIER
#define BARRIER
#endif

#include "limits.h"

#define OBJ_ADJ

#include "cplex.h"
#define STDIO_H_included
#include "nlp.h"
#include "getstub.h"
#include "avltree.h"
#include "signal.h"

/*  {Modifications for CPLEX 7.1 */
#ifndef CPX_PARAM_WORKDIR
#define CPX_PARAM_WORKDIR CPX_PARAM_NODEFILEDIR
#endif
#ifndef CPX_PARAM_PDSWITCH
#define CPX_PARAM_PDSWITCH 1063
#endif
/* End Modifications for CPLEX 7.1} */

/* {Modifications for CPLEX 8 */
#undef Const
#undef BARRIER_FOR_QP
#ifdef CPXERR_Q_NEG_ZERO_COMP	/* CPLEX 7.1 */
#define CPX_STAT_INFEASIBLE CPX_INFEASIBLE
#define CPX_STAT_UNBOUNDED CPX_UNBOUNDED
#define CPXCENVptr CPXENVptr
#define Const /*nothing*/
#define BARRIER_FOR_QP
#else
#define Const const
#endif
/* End Modifications for CPLEX 8} */

/* {Hidden params of CPLEX 6.5.1 */
#ifndef CPX_PARAM_OLDPRICING
#define CPX_PARAM_OLDPRICING 1054
#endif
#ifndef CPX_PARAM_EPSAGG
#define CPX_PARAM_EPSAGG 1055
#endif
/* End hidden params} */

#ifndef Sig_ret_type
#define Sig_ret_type void
#define SigRet /*nothing*/
#endif

#ifndef SigRet
#define SigRet /*nothing*/
#endif

#ifndef Sig_func_type
typedef void sig_func_type ANSI((int));
#endif

#ifdef _WIN32
/* Win32 is too hard. */
#undef  KEEP_BANNER
#define KEEP_BANNER
#else
#ifndef KEEP_BANNER
#include "fcntl.h"
#include "unistd.h"
#endif
#endif


static ASL *asl;
 static double Times[4];
#if CPX_VERSION >= 12050000
 static double DTimes[4];
 static int DTimes_failed, num_cores;
#endif


typedef struct cpxlp cpxlp;

 typedef struct
mint_values {
	int L;
	int U;
	int val;
	} mint_values;

 typedef struct
mdbl_values {
	double L;
	double U;
	double val;
	} mdbl_values;

 enum { /* sf_mint f values */
	set_crossover	=  0,
	set_dualthresh	=  1,
	set_netopt	=  2,
	set_objno	=  3,
	set_conpri	=  4,
	set_objpri	=  5,
	set_prestats	=  6,
	set_sos2	=  7,
	set_timing	=  8,
	set_iis		=  9,
	set_mipststat	= 10,
	set_mipstval	= 11,
	set_namernd	= 12,
	set_sos		= 13,
	set_mipcuts	= 14,
	set_round	= 15,
	set_rays	= 16,
	set_mipbasis	= 17,
	set_basis_cond	= 18,
	set_retmipgap	= 19,
	set_feasopt	= 20,
	set_feasoptobj	= 21,
	set_lazy	= 22,
	set_populate	= 23,
	set_pooldual	= 24,
	set_resolve	= 25,
	set_cutstats	= 26,
	set_incompat	= 27
	};
#ifdef CPX_PARAM_FEASOPTMODE /* >= 9.2b */
#define Uselazy
#endif

 enum { /* sf_mdbl f values */
	set_dual_ratio = 0
	};

 static mint_values
mint_val[28] = {
	/* set_crossover */	{0, 2, 1},
	/* set_dualthresh */	{-0x7fffffff, 0x7fffffff, 0},
	/* set_netopt */	{0, 2, 1},
	/* set_objno */		{0, 0/*n_obj*/,	1},
	/* set_conpri */	{0, 0x7fffffff, 1},
	/* set_objpri */	{0, 0x7fffffff, 2},
	/* set_prestats */	{0, 1, 0},
	/* set_sos2 */		{0, 1, 1},
#if CPX_VERSION >= 12050000
	/* set_timing */	{0, 0x3f, 0},
#else
	/* set_timing */	{0, 3, 0},
#endif
	/* set_iis */		{0, 3, 0},
	/* set_mipststat */	{0, 1, 1},
	/* set_mipstval */	{0, 2, 1},
	/* set_namernd */	{0, 1, 0},
	/* set_sos */		{0, 1, 1},
	/* set_mipcuts */	{-1, 2, 0},
	/* set_round */		{0, 15, 1},
	/* set_rays */		{0, 3, 3},
	/* set_mipbasis */	{0, 1, -1},
	/* set_basis_cond */	{0, 1, 0},
	/* set_retmipgap */	{0, 7, 0},
	/* set_feasopt */	{0, 2, 0},
	/* set_feasoptobj */	{1, 3, 1},
	/* set_lazy */		{0, 3, 3},
	/* set_populate */	{0, 2, 0},
	/* set_pooldual */	{0, 1, 0},
	/* set_resolve */	{0, 1, 1},
	/* set_cutstats */	{0, 1, 0},
	/* set_incompat */	{0, 2, 1}
	};

 static mdbl_values
mdbl_val[] = {
	/* set_dual_ratio */	{1., 1e30, 3.}
	};

#define crossover	mint_val[0].val
#define dual_thresh	mint_val[1].val
#define use_netopt	mint_val[2].val
#define objno		mint_val[3].val
#define conpri		mint_val[4].val
#define objpri		mint_val[5].val
#define prestats	mint_val[6].val
#define sos2		mint_val[7].val
#define time_flag	mint_val[8].val
#define want_iis	mint_val[9].val
#define mipststat	mint_val[10].val
#define mipstval	mint_val[11].val
#define Nameround	mint_val[12].val
#define sos		mint_val[13].val
/* #define mipcuts	mint_val[14].val */ /* only set_mipcuts is used */
#define Round		mint_val[15].val
#define rays		mint_val[16].val
#define mipbasis	mint_val[17].val
#define basis_cond	mint_val[18].val
#define retmipgap	mint_val[19].val
#define want_feasopt	mint_val[20].val
#define feasoptobj	mint_val[21].val
#define uselazy		mint_val[22].val
#define populate	mint_val[23].val
#define pooldual	mint_val[24].val
#define Resolve		mint_val[25].val
#define cutstats	mint_val[26].val
#define Incompat	mint_val[27].val
#define dual_ratio	mdbl_val[0].val

 static int hybmethod = CPX_ALG_PRIMAL;
 static int netiters = -1;
 static CPXFILEptr Logf;
 static char cplex_version[] = "AMPL/CPLEX with bad license\0\nAMPL/CPLEX Driver Version 20130531\n";
 static char *baralgname, *endbas, *endsol, *endtree, *endvec, *logfname;
 static char *paramfile, *poolstub, *pretunefile, *pretunefileprm;
 static char *startbas, *startsol, *starttree, *startvec, *tunefile, *tunefileprm;
 static char *tunefix, *tunefixfile, *workfiledir, *wrtfname;
 static int bestnode, breaking, costsens, lpoptalg, mbas, method;
 static int nbas, netopting, objadj, objsen, relax, zap_lpcbf;
 static int aggtries, net_status, net_nodes, net_arcs;
 static double obj_adj;
 static char *algname = "";
 static real intwarn_tol = 1e-9;

 typedef struct { char *msg; int code; int wantobj; }
Sol_info;
 static Sol_info solinfo[] = {
	 { "unrecoverable failure or licensing problem", 500, 0 },
	 { "optimal solution", 000, 1 },
	 { "infeasible problem", 200, 0 },			/* 2 */
	 { "unbounded problem", 300, 0 },
	 { "phase II objective limit exceeded", 400, 1 },	/* 4 */
	 { "iteration limit in phase II", 401, 1 },
	 { "iteration limit in phase I", 402, 1 },		/* 6 */
	 { "time limit in phase II", 403, 1 },
	 { "time limit in phase I", 404, 1 },			/* 8 */
	 { "infeasible with phase II singularities", 201, 0 },
	 { "infeasible with phase I singularities", 202, 0 },
	 { "optimal with unscaled infeasibilities", 110, 1 },	/* 11 */
	 { "aborted in phase II", 501, 1 },
	 { "aborted in phase I", 502, 0 },
	 { "aborted in barrier, dual infeasible", 503, 0 },	/* 14* */
	 { "aborted in barrier, primal infeasible", 504, 0 },
	 { "aborted in barrier, primal and dual infeasible", 505, 0 },
	 { "aborted in barrier, primal and dual feasible", 506, 0 },
	 { "aborted in crossover", 507, 0 },			/* 18 */

	 { "converged, dual feasible, primal infeasible", 204, 1 }, /*32 --> 19*/
	 { "converged, primal feasible, dual infeasible", 301, 1 },
	 { "converged, primal and dual infeasible", 205, 1 },
	 { "primal objective limit reached", 405, 1 },		/* 35 --> 22 */
	 { "dual objective limit reached", 406, 1 },
	 { "primal has unbounded optimal face", 001, 1 },
	 { "best solution found, primal-dual feasible", 100, 1 },/* 38 --> 25*/
	 { "best solution found, primal infeasible", 206, 1 },
	 { "best solution found, dual infeasible", 302, 1 },
	 { "best solution found, primal-dual infeasible", 207, 1 }, /* 41 --> 28 */
	 { "solution found, numerical difficulties", 508, 1},
	 { "solution found, inconsistent equations", 509, 1},  /* 43 --> 30*/

	 { "infeasible or unbounded in presolve", 208, 0 },	/*1101*/
	 { "not licensed to solve MIP problems", 571, 0},	/* 32 */
	 { "not licensed to use the barrier algorithm", 572, 0 },
	 { "not licensed to solve MIQP or QCP problems", 573, 0 },
	 { "locally optimal solution of indefinite QP", 130, 1 }
#ifdef CPLEX_MIP
	 ,{ "optimal integer solution", 002, 1 },		/* 101 --> 36 */
	 { "optimal integer solution within mipgap or absmipgap", 003, 1 },
	 { "integer infeasible", 220, 0 },			/* 103 */
	 { "mixed-integer solutions limit", 420, 1 },
	 { "node limit with integer solution", 421, 1 },	/* 105 */
	 { "node limit with no integer solution", 410, 0 },
	 { "time limit with integer solution", 422, 1 },
	 { "time limit with no integer solution", 411, 0 },	/* 108 */
	 { "unrecoverable failure with integer solution", 520, 1 },
	 { "unrecoverable failure with no integer solution", 510, 0 },
	 { "treememory limit with integer solution", 423, 1 },
	 { "treememory limit with no integer solution", 412, 0 },/* 112 */
	 { "aborted, integer solution exists", 521, 1 },
	 { "aborted, no integer solution", 511, 0 },		/* 114 */
	 { "integer optimal with unscaled infeasibilities", 111, 1 },
	 { "out of memory, no tree; solution may exist", 523, 1 },/* 116 */
	 { "out of memory, no tree; no integer solution", 512, 1 },
#ifdef CPX_STAT_OPTIMAL	/* >= 8.0 */
	 { "integer unbounded ray", 320, 0},	/* 118 */
	 { "integer infeasible or unbounded in presolve", 209, 0 }, /* 119 --> 54 */
#else
	 { "node file limit with integer solution", 424, 1 },
	 { "node file limit with no integer solution", 413, 1 }, /* 119 --> 54 */
#endif /*CPX_STAT_OPTIMAL*/
	 { "optimal (non-)integer solution", 102, 1 },
	 { "optimal (non-)integer solution within mipgap or absmipgap", 103, 1 },
	 { "feasible relaxed sum in feasopt", 121, 1}, /* --> 57 */
	 { "optimal relaxed sum in feasopt", 122, 1},
	 { "feasible relaxed infeasibility count in feasopt", 123, 1},
	 { "optimal relaxed infeasibility count in feasopt", 124, 1},
	 { "feasible relaxed quadratic penalty in feasopt", 125, 1},
	 { "optimal relaxed quadratic penalty in feasopt", 126, 1} /* --> 62 */
#endif /*CPLEX_MIP*/
	};

 typedef struct
Cbfinfo {
	int disp, mipdisp, np[4], nx[2], nint, pres;
	char *xkind;
	} Cbfinfo;
 static Cbfinfo cbi;

/* Stuff to cope with current CPLEX political correctness... */

 static CPXENVptr Env;

#define disconnectchannel(a)	CPXdisconnectchannel(Env,a)
#define addfuncdest(a,b,c)	CPXaddfuncdest(Env,a,b,c)

 typedef struct
dims {
	double	*c;
	double	*x;
	double	*y;
	int	*cstat;
	int	*rstat;
	SufDesc	*csd;
	SufDesc *rsd;
	char	*rtype;
	/* Stuff for netopt... */
	int	*cs;
	int	*rs;
	/* for solution pools */
	int	npool;
	} dims;

 static char *
strcpy1(char *t, const char *s)
{
	while((*t = *s++))
		t++;
	return t;
	}

 static void
badretfmt(int rc, char *fmt, ...)
{
	va_list ap;
	char buf[4200], *s;
	int k;

	va_start(ap, fmt);
	k = vsprintf(buf, fmt, ap) + 1;
	if (rc) {
		solve_result_num = rc;
		memcpy(s = (char*)M1alloc(k), buf, k);
		asl->i.uinfo = s;
		}
	else
		fprintf(Stderr, "%s\n", buf);
	}

 static void
badret(char *what, int i, int rc)
{
	badretfmt(rc, "%s failed; error code %d.", what, i);
	if (rc)
		exit(1);
	}

#ifdef CPX_CON_INDICATOR /* >= 9.2b */
#define ALLOW_CLP ASL_allow_CLP

 typedef struct
Cpx_Info {
	ASL *asl;
	CPXLPptr cpx;
	CPXENVptr env;
	} Cpx_Info;

 static int
add_indic(void *v, int iv, int compl, int sense, int nz, int *ig, real *g, real rhs)
{
	return CPXaddindconstr(Env, (CPXLPptr)v, iv, compl, nz, rhs, "LGE"[sense], ig, g, 0);
	}

#else /*!CPX_CON_INDICATOR*/
#define ALLOW_CLP 0
#endif /*CPX_CON_INDICATOR*/

#ifdef BASDEBUG
 int basdebug;

 static void
Basedebug(char *name, int *cs, int *rs) /*!!!!*/
{
	FILE *f;
	int i;

	if (!basdebug)
		return;
	if (!(f = fopen(name,"w"))) {
		fprintf(Stderr, "Basedebug: cannot open \"%s\".\n", name);
		return;
		}
	for(i = 0; i < nbas; i++)
		fprintf(f, "cs[%d] = %d\n", i, cs[i]);
	for(i = 0; i < mbas; i++)
		fprintf(f, "rs[%d] = %d\n", i, rs[i]);
	fclose(f);
	}
#else
#define Basedebug(a,b,c)
#endif

 static void
bascopy(CPXLPptr cpx, dims *d, int n, int m, int *cmap, int *ncs,
	int *nrs, int *rmap)
{
	int *cs, i, j, *rs;

#ifdef BASDEBUG
	if (basdebug & 1)
		return;
#endif
	if ((cs = d->cs))
		rs = d->rs;
	else {
		cs = (int*)Malloc((nbas + mbas + 1)*sizeof(int));
		rs = cs + nbas;
		}
	for (i = 0; i < nbas; ++i)
		cs[i] = CPX_AT_LOWER;
	for(i = 0; i < mbas; ++i)
		rs[i] = CPX_BASIC;
	for(i = 0; i < m; ++i)
		if ((j = rmap[i]) >= 0)
			rs[j] = nrs[i];
	for(i = 0; i < n; ++i) {
		if ((j = cmap[i]) >= 0)
			cs[j] = ncs[i];
		else {
			j = -(j+1);
			rs[j] = ncs[i];
			}
		}
	Basedebug("net.bas",cs,rs);
	CPXcopybase(Env, cpx, cs, rs);
	if (!d->cs)
		free(cs);
	}

 static int
netopt(CPXLPptr cpx, dims *d, int *st, int *nodes, int *arcs, int *its)
{
	CPXNETptr net;
	int i, j, k, m, n;
	int *cmap, *ncs, *nrs, *rmap, *s;

	k = 0;
	if (!(net = CPXNETcreateprob(Env,&k,"embedded")) || k)
		return 1;
	m = mbas + 1; /* CPXNETextract may add a dummy node */
	n = mbas + nbas;
	cmap = (int*)Malloc(2*(m+n)*sizeof(int));
	rmap = cmap + n;
	ncs = rmap + m;
	nrs = ncs + n;
	if (CPXNETextract(Env, net, cpx, cmap, rmap)
	 || (*nodes = CPXNETgetnumnodes(Env,net)) <= 0) {
		free(cmap);
		return 1;
		}
	*arcs = CPXNETgetnumarcs(Env,net);
	CPXgetintparam(Env, CPX_PARAM_SIMDISPLAY, &k);
	if (!k) {
		netopting = 1; /* suppress netopt termination msg */
		CPXsetintparam(Env, CPX_PARAM_SIMDISPLAY, 1);
		}
	if ((s = d->cs)) {
		n = *arcs;
		for(i = 0; i < n; i++)
			ncs[i] = (j = cmap[i]) >= 0 ? s[j] : CPX_BASIC;
		n = *nodes;
		s = d->rs;
		for(i = 0; i < n; i++)
			nrs[i] = (j = rmap[i]) >= 0 ? s[j] : CPX_BASIC;
		if ((i = CPXNETcopybase(Env, net, ncs, nrs)))
			fprintf(Stderr, "Return %d from CPXNETcopybase.\n", i);
		}
	if ((i = CPXNETprimopt(Env, net))) {
		fprintf(Stderr, "Return %d from CPXNETprimopt\n", i);
		*st = 11;
#ifdef Student_Edition
#ifndef CPXERR_RESTRICTED_VERSION
#define CPXERR_RESTRICTED_VERSION CPXERR_PROMOTION_VERSION
#endif
		if (i == CPXERR_RESTRICTED_VERSION)
			*st = 12;
#endif
		}
	else {
		i = CPXNETgetstat(Env, net);
		if (i >= 1 && i <= 8) {
			*st = i;
			if (i == 4)
				goto bug;
			}
		else if (i >= 12 && i <= 13)
			*st = i - 3;
		else if (i == 19)
			*st = 4;
		else {
 bug:
			fprintf(Stderr,
				"Return %d from CPXNETgetstat\n", i);
			*st = 11;
			}
		*its = CPXNETgetitcnt(Env,net);
#ifdef BASDEBUG
		CPXNETbasewrite(Env,net,"net.bw");
#endif
		if (!CPXNETgetbase(Env, net, ncs, nrs)) {
			bascopy(cpx, d, *arcs, *nodes, cmap, ncs, nrs, rmap);
#ifdef BASDEBUG
			CPXmbasewrite(Env, cpx, "m.bw");
#endif
			}
		}
	netopting = 0;
	if (!k)
		CPXsetintparam(Env, CPX_PARAM_SIMDISPLAY, 0);
	free(cmap);
	CPXNETfreeprob(Env, &net);
	return 0;
	}

 typedef int (CPXPUBLIC *Optalg)(CPXCENVptr,  cpxlp*);
 static Optalg Optimize;

 static jmp_buf Jb;

#ifdef CPXERR_IN_INFOCALLBACK /* CPLEX 11 */
 static int breaknow;
#endif

 Sig_ret_type
intcatch(int n)
{
	printf("\n<BREAK> (cplex)\n", n);
	fflush(stdout);
#ifdef CPXERR_IN_INFOCALLBACK
	breaknow = 1;
#endif
	if (++breaking > 3)
		longjmp(Jb, 2);
	signal(SIGINT, intcatch);
	SigRet;
	}

 static int CPXPUBLIC
breakcallback(CPXCENVptr e, void *c, int w, void *h)
{
	Not_Used(e);
	Not_Used(c);
	Not_Used(w);
	Not_Used(h);
	return breaking > 1;
	}

#ifdef BARRIER

 static int CPXPUBLIC
Optimize2(CPXCENVptr e, cpxlp *c)
{
	return CPXhybbaropt(e, c, hybmethod);
	}

 static int CPXPUBLIC
Optimizebar(CPXCENVptr e, cpxlp *c)
{
	int rv;

	if (!(rv = CPXbaropt(e, c))) {
#ifndef NO_DEPRECATED
		int s;
		if (endvec && (s = CPXvecwrite(e, c, endvec)))
			printf("\n*** return %d from CPXvecwrite.\n", s);
#endif
		}
	return rv;
	}

 static void
set_baropt(VOID)
{
	baralgname = "";
	if (crossover) {
		Optimize = Optimize2;
		if (crossover == 1) {
			algname = "crossover primal ";
			hybmethod = CPX_ALG_PRIMAL;
			}
		else {
			algname = "crossover dual ";
			hybmethod = CPX_ALG_DUAL;
			}
		}
	else {
		Optimize = Optimizebar;
		algname = "";
		}
	}

 static int
startcomp(ASL *asl, int n0, int m, int nextra, int *ka, int *kal, int *ia,
	double *a, double *b, double *c,
	double **cdualp, double **cprimp, double **rdualp, double **rprimp)
{
	double *r, t, *x, *y;
	int i, j, k;
	int n = n0 - objadj;

	/* Why can't CPXcopystart do this calculation? */

	if (!X0 && !pi0)
		return 0;
	if (method > 0) {
		x = X0;
		y = pi0;
		}
	else {
		x = pi0;
		y = X0;
		}
	if (x) {
		if (method < 0 && (i = nextra + nranges) > 0) {
			r = (double *)M1alloc(n0*sizeof(double));
			memcpy(r, x, (j = n - i)*sizeof(double));
			memset(r+j, 0, i*sizeof(double)); /* Can we do better? */
			x = r;
			if (objadj)
				x[n] = obj_adj;
			}
		else if (objadj) {
			r = (double *)M1alloc(n0*sizeof(double));
			memcpy(r, x, n*sizeof(double));
			x = r;
			x[n] = obj_adj;
			}
		*cprimp = x;
		*rprimp = r = (double *)M1alloc(m*sizeof(double));
		memcpy(r, b, m*sizeof(double));
		for(i = 0; i < n; i++)
			if ((t = x[i])) {
				j = ka[i];
				for(k = j + kal[i]; j < k; j++)
					r[ia[j]] -= t*a[j];
				}
		}
	else
		*cprimp = *rprimp = 0;
	if (y) {
		if (method > 0 && objadj) {
			r = (double *)M1alloc(m*sizeof(double));
			memcpy(r, y, (m-1)*sizeof(double));
			r[m-1] = obj_adj;
			y = r;
			}
		*rdualp = y;
		*cdualp = r = (double*)M1alloc(n0*sizeof(double));
		for(i = 0; i < n0; i++) {
			t = c[i];
			j = ka[i];
			for(k = j + kal[i]; j < k; j++)
				t -= y[ia[j]]*a[j];
			r[i] = t;
			}
		}
	else
		*cdualp = *rdualp = 0;
	return 1;
	}
#endif /*BARRIER*/

#ifdef CPX_PARAM_FEASOPTMODE /* >= 9.2b */
 static int
feasopt(CPXENVptr e,  cpxlp *cpx)
{
	char *grptype;
	double *L, *U, *lb, *rhs, *rng, *ub;
	int *grpbeg, *grpind, i, j, k, n, nqc;

	i = want_feasopt + 2*feasoptobj - 3;
	if (CPXsetintparam(e, CPX_PARAM_FEASOPTMODE, i))
		printf("Could not set CPX_PARAM_FEASOPTMODE to %d\n", i); /* bug */
	if (!(nqc = nlc)) {
		n = 2*(mbas + nbas);
		lb = (double*)Malloc(n*sizeof(double));
		for(i = 0; i < n; i++)
			lb[i] = 1.; /* preferences */
		ub = lb + nbas;
		rhs = ub + nbas;
		rng  = rhs + mbas;
		i = CPXfeasopt(e, cpx, rhs, rng, lb, ub);
		}
	else {
		L = LUv;
		U = Uvx;
		k = nbas;
		n = mbas;
		for(i = 0; i < k; i++) {
			if (L[i] > negInfinity)
				n++;
			if (U[i] < Infinity)
				n++;
			}
		lb = (double*)Malloc(n*(sizeof(double)+2*sizeof(int)+1));
		grpbeg = (int*)(lb + n);
		grpind = grpbeg + n;
		grptype = (char*)(grpind + n);
		for(i = 0; i < n; i++) {
			lb[i] = 1.; /* preferences */
			grpbeg[i] = i;
			}
		for(i = j = 0; i < k; i++) {
			if (L[i] > negInfinity) {
				grpind[j] = i;
				grptype[j++] = CPX_CON_LOWER_BOUND;
				}
			if (U[i] < Infinity) {
				grpind[j] = i;
				grptype[j++] = CPX_CON_UPPER_BOUND;
				}
			}
		k = n_con - nqc;
		for(i = 0; i < k; i++) {
			grpind[j] = i;
			grptype[j++] = CPX_CON_LINEAR;
			}
		k = nqc;
		for(i = 0; i < k; i++) {
			grpind[j] = i;
			grptype[j++] = CPX_CON_QUADRATIC;
			}
		i = CPXfeasoptext(e, cpx, n, n, lb, grpbeg, grpind, grptype);
		}
	free(lb);
	return i;
	}
#endif /*CPX_PARAM_FEASOPTMODE*/

#ifdef MSDOS
#undef Stand_alone
#endif
#ifdef Student_Edition
#undef Stand_alone
#endif

 extern int cpxmain ANSI((int, char**));
typedef char *sfunc ANSI((char *, char *, int));

#ifdef CPLEX_MIP
#ifdef CPX_PARAM_IISIND /* version < 9.2b */

typedef int (CPXPUBLIC *Treeio)(CPXCENVptr,CPXLPptr,Const char*);

 static void
treeio(cpxlp *cpx, char *fname, char *what, Treeio tio)
{
	if ((*tio)(Env,cpx,fname)) {
		printf("Could not %s tree file \"%s\".\n", what, fname);
		need_nl = 0;
		}
	}
#endif
#endif

 enum { /* sf_known f values */
	set_barrier	= 0,
	set_dual	= 1,
	set_dualopt	= 2,
	set_max		= 3,
	set_min		= 4,
	set_primal	= 5,
	set_primalopt	= 6,
	set_relax	= 7,
	set_sens	= 8,
	set_bestnode	= 9,
	set_autoopt	= 10,
	set_siftopt	= 11,
	set_concurrentopt=12,
	set_bestbound	= 13
	};

 static char *
sf_known(Option_Info *oi, keyword *kw, char *v)
{
	switch(Intcast kw->info) {
#ifdef BARRIER
	 case set_barrier:
		set_baropt();
		break;
#endif

	 case set_dual:
		method = -1;
		break;

	 case set_dualopt:
		Optimize = CPXdualopt;
		algname = "dual ";
		baralgname = "";
		break;

	 case set_max:
		objsen = -1;
		break;

	 case set_min:
		objsen = 1;
		break;

	 case set_primal:
		method = 1;
		break;

	 case set_primalopt:
		Optimize = CPXprimopt;
		algname = baralgname = "";
		break;

	 case set_relax:
		relax = 1;
		break;

	 case set_sens:
		costsens = 1;
		break;

	 case set_bestbound:
		bestnode |= 1;
		break;

	 case set_bestnode:
		bestnode |= 2;
		break;

	 case set_autoopt:
		lpoptalg = CPX_ALG_AUTOMATIC;
 use_lpopt:
		Optimize = CPXlpopt;
		CPXsetintparam(Env, CPX_PARAM_LPMETHOD, lpoptalg);
		break;

#ifdef CPX_STAT_OPTIMAL	/* >= 8.0 */

	 case set_siftopt:
		lpoptalg = CPX_ALG_SIFTING;
		goto use_lpopt;

	 case set_concurrentopt:
		lpoptalg = CPX_ALG_CONCURRENT;
		goto use_lpopt;

#endif /*CPX_STAT_OPTIMAL*/
	 }
	return v;
	}

 static void
baddval(Option_Info *oi, keyword *kw, double t, double L, double U)
{
	printf("rejecting %s %g; must be between %g and %g\n",
		kw->name, t, L, U);
	badopt_ASL(oi);
	}

 static void
badival(Option_Info *oi, keyword *kw, int t, int L, int U)
{
	printf("rejecting %s %d; must be between %d and %d\n",
		kw->name, t, L, U);
	badopt_ASL(oi);
	}

 static char*
sf_mdbl(Option_Info *oi, keyword *kw, char *v)
{
	char *rv;
	double t;
	int i = Intcast kw->info;
	mdbl_values *m = mdbl_val + i;

	if (*v == '?' && v[1] <= ' ') {
		printf("%s=%g\n", kw->name, m->val);
		oi->option_echo &= ~ASL_OI_echothis;
		return v + 1;
		}
	t = strtod(v, &rv);
	if (rv == v) {
		printf("Expected a floating-point value for %s, not \"%s\"\n",
			kw->name, v);
		badopt_ASL(oi);
		return v;
		}
	if (t < m->L || t > m->U) {
		baddval(oi,kw,t,m->L,m->U);
		return rv;
		}
	m->val = t;
	return rv;
	}

 static char *
sf_mint(Option_Info *oi, keyword *kw, char *v)
{
	int t;
	char *rv;
	int i = Intcast kw->info;
	mint_values *m = mint_val + i;

	if (*v == '?' && v[1] <= ' ') {
		printf("%s=%d\n", kw->name, m->val);
		oi->option_echo &= ~ASL_OI_echothis;
		return v + 1;
		}
	t = (int)strtol(v, &rv, 10);
	if (rv == v) {
		printf("Expected an integer value for %s, not \"%s\"\n",
			kw->name, v);
		badopt_ASL(oi);
		return v;
		}
	if (t < m->L || t > m->U) {
		badival(oi,kw,t,m->L,m->U);
		return rv;
		}
	m->val = t;
#ifdef BARRIER
	if (i == set_crossover)
		set_baropt();
#endif
#ifdef CPLEX_MIP
	if (i == set_mipcuts) {
		static int op[] = {
			CPX_PARAM_CLIQUES,
			CPX_PARAM_COVERS,
			CPX_PARAM_DISJCUTS,
			CPX_PARAM_FLOWCOVERS,
			CPX_PARAM_FLOWPATHS,
			CPX_PARAM_FRACCUTS,
			CPX_PARAM_GUBCOVERS,
			CPX_PARAM_IMPLBD,
			CPX_PARAM_MIRCUTS
#ifdef CPX_PARAM_ZEROHALFCUTS
			,CPX_PARAM_ZEROHALFCUTS
#endif
			};
		int f;
		for(f = 0; f < sizeof(op)/sizeof(int); f++)
			CPXsetintparam(Env, op[f], t);
		}
#endif
	return rv;
	}

#ifdef CPXERR_PARAM_INCOMPATIBLE /*{*/
 static char *
incompatible(Option_Info *oi, keyword *kw, char *rv, char *v)
{
	oi->option_echo &= ~ASL_OI_echothis;
	if (Incompat) {
		printf("%s \"%s%s%.*s\" as incompatible "
			"with earlier parameter settings.\n",
			Incompat == 1 ? "Ignoring" : "Rejecting",
			kw->name, oi->eqsign, (int)(rv-v), v);
		if (Incompat == 2)
			badopt_ASL(oi);
		}
	return rv;
	}
#endif /*}*/

 static char *
sf_int(Option_Info *oi, keyword *kw, char *v)
{
	int f, k, t, z[3];
	char *rv;
	const char *what = kw->name;

	f = Intcast kw->info;
	if (*v == '?' && v[1] <= ' ') {
		CPXgetintparam(Env, f, &t);
		printf("%s=%d\n", what, t);
		oi->option_echo &= ~ASL_OI_echothis;
		return v + 1;
		}
	t = (int)strtol(v, &rv, 10);
	if (rv == v) {
		printf("Expected an integer value for %s, not \"%s\"\n",
			what, v);
		badopt_ASL(oi);
		}
	else if ((k = CPXsetintparam(Env, f, t))) {
#ifdef CPXERR_PARAM_INCOMPATIBLE
		if (k == CPXERR_PARAM_INCOMPATIBLE)
			return incompatible(oi, kw, rv, v);
#endif
		z[2] = 0;
		z[1] = 1;
		CPXinfointparam(Env, f, z, z+1, z+2);
		badival(oi, kw, t, z[1], z[2]);
		}
	return rv;
	}

 static char *
sf_ipar(Option_Info *oi, keyword *kw, char *v)
{
	int f, t, z[3];
	char *rv, *v0 = v;
	const char *what = kw->name;

	f = (int)strtol(v, &rv, 10);
	if (rv == v) {
		printf("Expected an integer parameter number after %s%s, not \"%s\"\n",
			what, oi->eqsign, v);
		badopt_ASL(oi);
		}
	for(v = rv; *v <= ' '; v++)
		if (!*v) {
			printf("Expected an integer value after %s%s\n", what, oi->eqsign);
			badopt_ASL(oi);
			}
	if (*v == '=')
		while(*++v <= ' ')
			if (!*v) {
				printf("Expected an integer value after %s%s%s\n",
					what, oi->eqsign, v0);
				badopt_ASL(oi);
				}
	if (*v == '?' && v[1] <= ' ') {
		CPXgetintparam(Env, f, &t);
		printf("%s=%d=%d\n", what, f, t);
		oi->option_echo &= ~ASL_OI_echothis;
		return v + 1;
		}
	t = (int)strtol(v, &rv, 10);
	if (rv == v) {
		printf("Expected an integer value for %s %d, not \"%s\"\n",
			what, f, v);
		badopt_ASL(oi);
		}
	else if (CPXsetintparam(Env, f, t)) {
		z[2] = 0;
		z[1] = 1;
		CPXinfointparam(Env, f, z, z+1, z+2);
		printf("rejecting %s=%d=%d; assigned value must be between %d and %d\n",
			what, f, t, z[1], z[2]);
		badopt_ASL(oi);
		}
	return rv;
	}

/* In case some reason surfaces to distinguish the old setzzzpar()	*/
/* calls from the old setzzzind calls (which had sfunc value sf_int and	*/
/* sf_int1, respectively, in the keywds table, we use the following	*/
/* #define and retain appearances of sf_int1 in the keywds table.	*/

#define sf_int1 sf_int

 static char *
sf_dbl(Option_Info *oi, keyword *kw, char *v)
{
	double t, z[3];
	char *rv;
	const char *what = kw->name;
	int f, k;

	f = Intcast kw->info;
	if (*v == '?' && v[1] <= ' ') {
		CPXgetdblparam(Env, f, &t);
		printf("%s=%g\n", what, t);
		oi->option_echo &= ~ASL_OI_echothis;
		return v + 1;
		}
	t = strtod(v, &rv);
	if (rv == v) {
		printf("Expected a numeric value for %s, not \"%s\"\n",
			what, v);
		badopt_ASL(oi);
		}
	else if ((k = CPXsetdblparam(Env, f, t))) {
#ifdef CPXERR_PARAM_INCOMPATIBLE
		if (k == CPXERR_PARAM_INCOMPATIBLE)
			return incompatible(oi, kw, rv, v);
#endif
		z[2] = 0;
		z[1] = 1;
		CPXinfodblparam(Env, f, z, z+1, z+2);
		printf("rejecting %s %g; must be between %g and %g\n",
			what, t, z[1], z[2]);
		badopt_ASL(oi);
		}
	return rv;
	}

/* sf_int2 and sf_dbl2 are variants of sf_int and sf_dbl for  */
/* synonyms, so getparinfo() will not deal with the synonyms. */

 static char *
sf_int2(Option_Info *oi, keyword *kw, char *v)
{ return sf_int(oi,kw,v); }

 static char *
sf_dbl2(Option_Info *oi, keyword *kw, char *v)
{ return sf_dbl(oi,kw,v); }


 static char *
sf_dpar(Option_Info *oi, keyword *kw, char *v)
{
	double t, z[3];
	int f;
	char *rv, *v0 = v;
	const char *what = kw->name;

	f = (int)strtol(v, &rv, 10);
	if (rv == v) {
		printf("Expected an integer parameter number after %s%s, not \"%s\"\n",
			what, oi->eqsign, v);
		badopt_ASL(oi);
		}
	for(v = rv; *v <= ' '; v++)
		if (!*v) {
			printf("Expected a floating-point value after %s%s\n", what, oi->eqsign);
			badopt_ASL(oi);
			}
	if (*v == '=')
		while(*++v <= ' ')
			if (!*v) {
				printf("Expected a floating-point value after %s%s%s\n",
					what, oi->eqsign, v0);
				badopt_ASL(oi);
				}
	if (*v == '?' && v[1] <= ' ') {
		CPXgetdblparam(Env, f, &t);
		printf("%s=%d=%g\n", what, f, t);
		oi->option_echo &= ~ASL_OI_echothis;
		return v + 1;
		}
	t = strtod(v, &rv);
	if (rv == v) {
		printf("Expected a floating-point value for %s %d, not \"%s\"\n",
			what, f, v);
		badopt_ASL(oi);
		}
	else if (CPXsetdblparam(Env, f, t)) {
		z[2] = 0;
		z[1] = 1;
		CPXinfodblparam(Env, f, z, z+1, z+2);
		printf("rejecting %s=%d=%g; assigned value must be between %g and %g\n",
			what, f, t, z[1], z[2]);
		badopt_ASL(oi);
		}
	return rv;
	}


 static char **file_name[19] = { &endbas, &endtree, &startbas, &starttree,
				&startsol, &endsol, &logfname, &wrtfname,
				&workfiledir, &poolstub, &paramfile,
				&pretunefile, &pretunefileprm, &tunefile,
				&tunefileprm, &tunefix, &tunefixfile,
				&startvec, &endvec };

enum {	/* sf_char f values */
	set_endbas	= 0,
	set_endtree	= 1,
	set_startbas	= 2,
	set_starttree	= 3,
	set_startsol	= 4,
	set_endsol	= 5,
	set_logname	= 6,
	set_wrtfname	= 7,
	set_workfiledir = 8,
	set_poolstub	= 9,
	set_paramfile	= 10,
	set_pretunefile	= 11,
	set_pretunefileprm = 12,
	set_tunefile	= 13,
	set_tunefileprm	= 14,
	set_tunefix	= 15,
	set_tunefixfile	= 16,
	set_startvector	= 17,
	set_endvector	= 18
	};

 static char *
sf_char(Option_Info *oi, keyword *kw, char *v)
{
	char *rv, *t;
	int f, q;
	size_t n;

	if (!*v) {
		printf("rejecting %s: no following file name\n", kw->name);
		badopt_ASL(oi);
		return v;
		}
	f = Intcast kw->info;
	if (*v == '?' && v[1] <= ' ') {
		if (!(t = *file_name[f]))
			t = "";
		printf("%s=%s\n", kw->name, t);
		oi->option_echo &= ~ASL_OI_echothis;
		return v + 1;
		}
	q = 0;
	if (*v == '\'' || *v == '"') {
		q = *v++;
		for(rv = v; *rv != q; ++rv) {
			if (!*rv) {
				printf("Missing final %c in file name for %s.\n",
					q, kw->name);
 badname:
				badopt_ASL(oi);
				return rv+1;
				}
			}
		if (v == rv) {
			printf("Empty file name for %s.\n", kw->name);
			goto badname;
			}
		if (*++rv > ' ') {
			printf("Unexpected '%c' after closing $c of file name for %s.\n",
				*rv, q, kw->name);
			goto badname;
			}
		n = rv - v;
		}
	else {
		for(rv = v; *++rv > ' ';);
		n = rv - v + 1;
		}
	t = M1alloc(n);
	strncpy(t, v, --n);
	t[n] = 0;
	*file_name[f] = t;
	if (f == set_logname) {
		if (!(Logf = CPXfopen(t,"w"))) {
			printf("Cannot open logfile \"%s\"\n", t);
			badopt_ASL(oi);
			}
		else
			CPXsetlogfile(Env, Logf);
		}
	return rv;
	}

 static Option_Info Oinfo;

 static char *
sf_par(Option_Info *oi, keyword *kw, char *v)
{
	FILE *f;
	char buf[4096], *s;

	if (*v == '?' && v[1] <= ' ') {
		printf("# %s is processed when specified.\n", kw->name);
		return v + 1;
		}
	paramfile = 0;
	v = sf_char(oi, kw, v);
	if (!paramfile)
		return v;
	if (!(f = fopen(paramfile, "r"))) {
		printf("Cannot open %s \"%s\"\n", kw->name, paramfile);
		badopt_ASL(oi);
		return v;
		}
	while(fgets(buf, sizeof(buf), f))
		for(s = buf; *s; s = get_opt_ASL(&Oinfo, s));
	fclose(f);
	if (Oinfo.n_badopts) {
		if (amplflag)
			badretfmt(564, "Error reading paramfile \"%s\".", paramfile);
		exit(1);
		}
	return v;
	}

#if CPX_VERSION >= 1000 /*{*/
 static char *
sf_parm(Option_Info *oi, keyword *kw, char *v)
{
	if (*v == '?' && v[1] <= ' ') {
		printf("# %s is processed when specified.\n", kw->name);
		return v + 1;
		}
	paramfile = 0;
	v = sf_char(oi, kw, v);
	if (!paramfile)
		return v;
	if (CPXreadcopyparam(Env, paramfile)) {
		if (amplflag)
			badretfmt(564, "Error reading paramfileprm \"%s\".", paramfile);
		exit(1);
		}
	return v;
	}
#endif /*}*/

#define VP (Char*)

 struct option_word {
	char *name;
	sfunc *sf;
	int ival;
	};
 typedef struct option_word option_word;

 static keyword keywds[] = {	/* must be in alphabetical order */

	/* Undocumented keywords start with underscore... */

	{ "_aggsort",	sf_int,		VP 1061 /*CPX_PARAM_PREAGGSORT*/},
	{ "_aggtolerance", sf_dbl,	VP 1055 /*CPX_PARAM_EPSAGG*/},
	{ "_cancel",	sf_int,		VP 1071 /*CPX_PARAM_PRECANCEL*/},
	{ "_clique",	sf_int,		VP 1072 /*CPX_PARAM_PRECLIQUE*/},
	{ "_cliquetablesize", sf_dbl,	VP 2064 /*CPX_PARAM_CLIQUETABLESZ*/},
	{ "_domination", sf_int,	VP 2038 /*CPX_PARAM_COLDOMIND*/},
	{ "_effslack",	sf_int,		VP 1042 /*CPX_PARAM_EFFSLACKIND*/},
	{ "_factormem",	sf_dbl,		VP 3021 /*CPX_PARAM_BARFACTMEM*/},
	{ "_flip",	sf_int,		VP 1051 /*CPX_PARAM_FLIPIND*/},
	{ "_hfeasibility", sf_dbl,	VP 1050 /*CPX_PARAM_EPRHS_H*/},
	{ "_hoptimality", sf_dbl,	VP 1049 /*CPX_PARAM_EPOPT_H*/},
	{ "_insubtree",	sf_int,		VP 2063 /*CPX_PARAM_INSUBTREE*/},
	{ "_kernel",	sf_int,		VP 3020 /*CPX_PARAM_BARKERNEL*/},
	{ "_knapcoeff",	sf_int,		VP 1059 /*CPX_PARAM_KNAPCOERED*/},
	{ "_localcovers", sf_int,	VP 2061 /*CPX_PARAM_LOCALCOVERS*/},
	{ "_logparams",	sf_int,		VP 1075 /*CPX_PARAM_LOGPARAMS*/},
	{ "_memfact",	sf_dbl,		VP 1045 /*CPX_PARAM_PREMEMFACT*/},
	{ "_memsave",	sf_int,		VP 1060 /*CPX_PARAM_PREMEMSAVE*/},
	{ "_minstuck",	sf_int,		VP 3023 /*CPX_PARAM_BARMINSTUCK*/},
	{ "_oldpricing", sf_int,	VP 1054 /*CPX_PARAM_OLDPRICING*/},
	{ "_oldqpfactor", sf_int,	VP 3024 /*CPX_PARAM_OLDQPFACTOR*/},
	{ "_oldratio",	sf_int,		VP 1068 /*CPX_PARAM_OLDRATIO*/},
	{ "_orderthreads", sf_int,	VP 3022 /*CPX_PARAM_BARORDERTHREADS*/},
	{ "_primalstart", sf_dbl,	VP 3005 /*CPX_PARAM_BARPSTART*/},
	{ "_probe",	sf_int,		VP 1070 /*CPX_PARAM_PREPROBE*/},
	{ "_recurseheur", sf_int,	VP 2062 /*CPX_PARAM_RECURSEHEUR*/},
#ifdef CPX_PARAM_IISIND /* version < 9.2b */
	{ "_rowsdense",	sf_int,		VP 3015 /*CPX_PARAM_BARROWSDEN*/},
#endif
	{ "_splitrow",	sf_int,		VP 1079 /*CPX_PARAM_PRESPLITROW*/},
	{ "_svbound",	sf_int,		VP 1069 /*CPX_PARAM_SVBNDSTR*/},

	/* Documented keywords... */

#ifdef CPLEX_MIP
	{ "absmipgap",	sf_dbl,		VP CPX_PARAM_EPAGAP },
#endif
	{ "advance",	sf_int1,	VP CPX_PARAM_ADVIND },
#ifdef CPLEX_MIP
	{ "aggcutlim",	sf_int,		VP CPX_PARAM_AGGCUTLIM },
#endif
	{ "aggfill",	sf_int2,	VP CPX_PARAM_AGGFILL },
	{ "agglim",	sf_int,		VP CPX_PARAM_AGGFILL },
	{ "aggregate",	sf_int1,	VP CPX_PARAM_AGGIND },
	{ "aggtol",	sf_dbl,		VP CPX_PARAM_EPSAGG },
	{ "autoopt",	sf_known,	VP set_autoopt },
	{ "autopt",	sf_known,	VP set_autoopt },
#ifdef CPLEX_MIP
#ifdef CPX_PARAM_AUXROOTTHREADS
	{ "auxrootthreads", sf_int,	VP CPX_PARAM_AUXROOTTHREADS },
#endif
	{ "backtrack",	sf_dbl,		VP CPX_PARAM_BTTOL },
#endif
#ifdef BARRIER
	{ "baralg",	sf_int,		VP CPX_PARAM_BARALG },
	{ "barcorr",	sf_int,		VP CPX_PARAM_BARMAXCOR },
	{ "bardisplay",	sf_int2,	VP CPX_PARAM_BARDISPLAY },
	{ "bargrowth",	sf_dbl,		VP CPX_PARAM_BARGROWTH },
	{ "bariterlim",	sf_int,		VP CPX_PARAM_BARITLIM },
	{ "barobjrange", sf_dbl,	VP CPX_PARAM_BAROBJRNG },
	{ "baropt",	sf_known,	VP set_barrier },
#ifdef CPX_PARAM_BAROOC
	{ "baroutofcore",sf_int,	VP CPX_PARAM_BAROOC },
#endif
	{ "barstart",	sf_int,		VP CPX_PARAM_BARSTARTALG },
	{ "barstartalg",sf_int,		VP CPX_PARAM_BARSTARTALG },
#ifdef CPX_PARAM_BARTHREADS
	{ "barthreads",	sf_int,		VP CPX_PARAM_BARTHREADS },
#endif
#ifdef CPX_PARAM_BARVARUP
	{ "barvarup",	sf_dbl,		VP CPX_PARAM_BARVARUP },
#endif
#endif /* BARRIER */
#ifdef BASDEBUG
	{ "basdebug", I_val,		VP &basdebug },
#endif
	{ "basis_cond",	sf_mint,	VP set_basis_cond },
	{ "basisinterval", sf_int,	VP CPX_PARAM_BASINTERVAL },
#ifdef CPLEX_MIP
	{ "bbinterval",	sf_int,		VP CPX_PARAM_BBINTERVAL },
	{ "bestbound",	sf_known,	VP set_bestbound },
	{ "bestnode",	sf_known,	VP set_bestnode },
	{ "boundstr",	sf_int,		VP CPX_PARAM_BNDSTRENIND },
	{ "branch",	sf_int,		VP CPX_PARAM_BRDIR },
	{ "branchdir",	sf_int,		VP CPX_PARAM_BRDIR },
	{ "cliquecuts",	sf_int2,	VP CPX_PARAM_CLIQUES },
	{ "cliques",	sf_int,		VP CPX_PARAM_CLIQUES },
#endif
	{ "clocktype",	sf_int,		VP CPX_PARAM_CLOCKTYPE },
#ifdef CPLEX_MIP
	{ "coeffreduce", sf_int1,	VP CPX_PARAM_COEREDIND },
#endif
#ifdef BARRIER
	{ "comptol",	sf_dbl,		VP CPX_PARAM_BAREPCOMP },
#endif
	{ "concurrent",	sf_known,	VP set_concurrentopt },
	{ "concurrentopt",sf_known,	VP set_concurrentopt },
#ifdef CPX_PARAM_CONFLICTDISPLAY
	{ "conflictdisplay", sf_int,	VP CPX_PARAM_CONFLICTDISPLAY },
#endif
#ifdef CPLEX_MIP
	{ "covercuts",	sf_int2,	VP CPX_PARAM_COVERS },
	{ "covers",	sf_int,		VP CPX_PARAM_COVERS },
#endif
	{ "crash",	sf_int,		VP CPX_PARAM_CRAIND },
#ifdef BARRIER
	{ "crossover",	sf_mint,	VP set_crossover },
#endif
#ifdef CPLEX_MIP
	{ "cutpass",	sf_int,		VP CPX_PARAM_CUTPASS },
	{ "cutsfactor",	sf_dbl,		VP CPX_PARAM_CUTSFACTOR },
#if CPX_VERSION >= 1100
	{ "cutstats",	sf_mint,	VP set_cutstats },
#endif
#endif
#ifdef BARRIER
	{ "dense",	sf_int2,	VP CPX_PARAM_BARCOLNZ },
	{ "densecol",	sf_int,		VP CPX_PARAM_BARCOLNZ },
#endif
	{ "dependency",	sf_int1,	VP CPX_PARAM_DEPIND },
#ifdef CPX_PARAM_DETTILIM
	{ "dettimelim",	sf_dpar,	VP CPX_PARAM_DETTILIM },
#endif
	{ "dgradient",	sf_int,		VP CPX_PARAM_DPRIIND },
#ifdef CPLEX_MIP
	{ "disjcuts",	sf_int,		VP CPX_PARAM_DISJCUTS },
#endif
	{ "display",	sf_int2,	VP CPX_PARAM_SIMDISPLAY },
	{ "doperturb",	sf_int1,	VP CPX_PARAM_PERIND },
	{ "dparam",	sf_dpar,	0 },
	{ "dual",	sf_known,	VP set_dual },
	{ "dualopt",	sf_known,	VP set_dualopt },
	{ "dualratio",	sf_mdbl,	VP set_dual_ratio },
	{ "dualthresh",	sf_mint,	VP set_dualthresh },
#ifdef CPX_PARAM_EACHCUTLIM
	{ "eachcutlim",	sf_int,		VP CPX_PARAM_EACHCUTLIM },
#endif
	{ "endbasis",	sf_char,	VP set_endbas },
#ifdef CPLEX_MIP
#ifdef CPX_PARAM_IISIND /* version < 9.2b */
	{ "endtree",	sf_char,	VP set_endtree },
#endif
#endif
	{ "endsol",	sf_char,	VP set_endsol },
#if defined(BARRIER) && !defined(NO_DEPRECATED)
	{ "endvector", 	sf_char,	VP set_endvector },
#endif
	{ "feasibility", sf_dbl,	VP CPX_PARAM_EPRHS },
#ifdef CPX_PARAM_FEASOPTMODE /* >= 9.2b */
	{ "feasopt",	sf_mint,	VP set_feasopt },
	{ "feasoptobj",	sf_mint,	VP set_feasoptobj },
#endif
	{ "file",	sf_char,	VP set_wrtfname },
#ifdef CPX_PARAM_FINALFACTOR
	{ "finalfactor", sf_int,	VP CPX_PARAM_FINALFACTOR },
#endif
#ifdef CPLEX_MIP
	{ "flowcuts",	sf_int,		VP CPX_PARAM_FLOWCOVERS },
	{ "flowpathcuts",sf_int,	VP CPX_PARAM_FLOWPATHS },
#ifdef CPX_PARAM_FPHEUR
	{ "fpheur",	sf_int,		VP CPX_PARAM_FPHEUR },
#endif
#ifndef NO_CPLEX66 /* for versions prior to CPLEX 6.6 */
	{ "fraccand",   sf_int,	        VP CPX_PARAM_FRACCAND },
	{ "fraccuts",   sf_int,	        VP CPX_PARAM_FRACCUTS },
	{ "fracpass",   sf_int,	        VP CPX_PARAM_FRACPASS },
	{ "fractionalcuts", sf_int,	VP CPX_PARAM_FRACCUTS },
#endif
#endif
#ifdef BARRIER
	{ "growth",	sf_dbl,		VP CPX_PARAM_BARGROWTH }, /*== bargrowth*/
#endif
#ifdef CPLEX_MIP
	{ "gubcuts",	sf_int,		VP CPX_PARAM_GUBCOVERS },
	{ "heurfreq",	sf_int,		VP CPX_PARAM_HEURFREQ },
#ifdef CPX_PARAM_HEURISTIC
	{ "heuristic",	sf_int,		VP CPX_PARAM_HEURISTIC },
#endif
	{ "heuristicfreq", sf_int,	VP CPX_PARAM_HEURFREQ },
#endif
	{ "iisfind",	sf_mint,	VP set_iis },
#ifdef CPLEX_MIP
	{ "impliedcuts", sf_int,	VP CPX_PARAM_IMPLBD },
#endif
#ifdef CPXERR_PARAM_INCOMPATIBLE
	{ "incompat", sf_mint,		VP set_incompat },
#endif
#ifdef CPLEX_MIP
	{ "integrality", sf_dbl,	VP CPX_PARAM_EPINT },
	{ "intwarntol", D_val,		VP &intwarn_tol },
#endif
	{ "iparam",	sf_ipar,	0 },
	{ "iterations",	sf_int,		VP CPX_PARAM_ITLIM },
	{ "iterlim",	sf_int,		VP CPX_PARAM_ITLIM },
#ifdef Uselazy
	{ "lazy",	sf_mint,	VP set_lazy },
#endif
#ifdef CPLEX_MIP
#ifdef CPX_PARAM_LBHEUR
	{ "lbheur",	sf_int,		VP CPX_PARAM_LBHEUR },
#endif
#endif
	{ "limitperturb",sf_int2,	VP CPX_PARAM_PERLIM },
	{ "logfile",	sf_char,	VP set_logname },
#ifdef CPLEX_MIP
	{ "lowercutoff", sf_dbl,	VP CPX_PARAM_CUTLO },
#endif
	{ "lowerobj",	sf_dbl,		VP CPX_PARAM_OBJLLIM },
	{ "lowerobjlim", sf_dbl,	VP CPX_PARAM_OBJLLIM },
	{ "lpdisplay",  sf_int2,	VP CPX_PARAM_SIMDISPLAY },
	{ "lpiterlim",  sf_int2,	VP CPX_PARAM_ITLIM },
	{ "lptimelim",	sf_dbl2,	VP CPX_PARAM_TILIM },
	{ "markowitz",	sf_dbl,		VP CPX_PARAM_EPMRK },
	{ "maximize",	sf_known,	VP set_max },
#ifdef CPX_PARAM_MCFCUTS
	{ "mcfcuts",	sf_int,		VP CPX_PARAM_MCFCUTS },
#endif
#ifdef CPX_PARAM_MEMORYEMPHASIS
	{ "memoryemphasis", sf_int,	VP CPX_PARAM_MEMORYEMPHASIS },
#endif
	{ "minimize",	sf_known,	VP set_min },
#ifdef CPLEX_MIP
	{ "mipalg",	sf_int,		VP CPX_PARAM_SUBALG },
	{ "mipalgorithm", sf_int,	VP CPX_PARAM_SUBALG },
	{ "mipbasis",	sf_mint,	VP set_mipbasis },
	{ "mipcrossover",sf_int,	VP CPX_PARAM_BARCROSSALG },
	{ "mipcuts",	sf_mint,	VP set_mipcuts },
	{ "mipdisplay",	sf_int2,	VP CPX_PARAM_MIPDISPLAY },
	{ "mipemphasis",sf_int,		VP CPX_PARAM_MIPEMPHASIS },
	{ "mipgap",	sf_dbl,		VP CPX_PARAM_EPGAP },
	{ "mipinterval", sf_int,	VP CPX_PARAM_MIPINTERVAL },
#ifdef CPX_PARAM_MIPKAPPASTATS
	{ "mipkappa",	sf_int,		VP CPX_PARAM_MIPKAPPASTATS },
#endif
	{ "mipordertype",sf_int2,	VP CPX_PARAM_MIPORDTYPE },
#ifdef CPX_PARAM_MIPSEARCH
	{ "mipsearch",	sf_int,		VP CPX_PARAM_MIPSEARCH },
#endif
	{ "mipsolutions", sf_int,	VP CPX_PARAM_INTSOLLIM },
	{ "mipstart",	sf_mint,	VP set_mipstval },
	{ "mipstartalg", sf_int,	VP CPX_PARAM_STARTALG },
	{ "mipstartstatus", sf_mint,	VP set_mipststat },
	{ "mipstartvalue", sf_mint,	VP set_mipstval },
	{ "mipsubalg",	sf_int,		VP CPX_PARAM_SUBALG },
#ifdef CPX_PARAM_MIPTHREADS
	{ "mipthreads",	sf_int,		VP CPX_PARAM_MIPTHREADS },
#endif
#ifdef CPX_PARAM_MIQCPSTRAT
	{ "miqcpstrat",	sf_int,		VP CPX_PARAM_MIQCPSTRAT },
#endif
	{ "mircuts",	sf_int,		VP CPX_PARAM_MIRCUTS },
#endif
	{ "nameround",	sf_mint,	VP set_namernd },
	{ "netdisplay",	sf_int,		VP CPX_PARAM_NETDISPLAY },
	{ "netfeasibility", sf_dbl,	VP CPX_PARAM_NETEPRHS },
	{ "netfind",	sf_int,		VP CPX_PARAM_NETFIND },
	{ "netfinder",	sf_int,		VP CPX_PARAM_NETFIND },
	{ "netiterations", sf_int,	VP CPX_PARAM_NETITLIM },
	{ "netopt",	sf_mint,	VP set_netopt },
	{ "netoptimality", sf_dbl,	VP CPX_PARAM_NETEPOPT },
	{ "netpricing",	sf_int,		VP CPX_PARAM_NETPPRIIND },
#ifdef CPLEX_MIP
	{ "node",	sf_int2,	VP CPX_PARAM_NODELIM },
	{ "nodefile",	sf_int,		VP CPX_PARAM_NODEFILEIND },
	{ "nodefiledir",sf_char,	VP set_workfiledir },
#ifdef CPX_PARAM_NODEFILELIM
	{ "nodefilelim", sf_dbl,	VP CPX_PARAM_NODEFILELIM },
	{ "nodefilesize",sf_dbl,	VP CPX_PARAM_NODEFILELIM },
#endif
	{ "nodelim",	sf_int2,	VP CPX_PARAM_NODELIM },
	{ "nodes",	sf_int,		VP CPX_PARAM_NODELIM },
	{ "nodesel",	sf_int,		VP CPX_PARAM_NODESEL },
	{ "nodeselect",	sf_int,		VP CPX_PARAM_NODESEL },
#endif
#ifdef CPX_PARAM_NUMERICALEMPHASIS
	{ "numericalemphasis",	sf_int,	VP CPX_PARAM_NUMERICALEMPHASIS },
#endif
#ifdef CPLEX_MIP
	{ "objdifference", sf_dbl,	VP CPX_PARAM_OBJDIF },
#endif
	{ "objno",	sf_mint,	VP set_objno },
	{ "oldpricing",	sf_int,		VP CPX_PARAM_OLDPRICING },
	{ "optimality",	sf_dbl,		VP CPX_PARAM_EPOPT },
	{ "optimize",	sf_known,	VP set_primalopt },
#ifdef BARRIER
	{ "ordering",	sf_int,		VP CPX_PARAM_BARORDER },
#endif
#ifdef CPLEX_MIP
	{ "ordertype",	sf_int,		VP CPX_PARAM_MIPORDTYPE },
#endif
	{ "outlev",	sf_int2,	VP CPX_PARAM_SIMDISPLAY },
#ifdef CPX_PARAM_PARALLELMODE
	{ "parallelmode", sf_int,	VP CPX_PARAM_PARALLELMODE },
#endif
	{ "paramfile",	sf_par,		VP set_paramfile },
#if CPX_VERSION >= 1000
	{ "paramfileprm", sf_parm,	VP set_paramfile },
#endif
	{ "pdswitch",	sf_int,		VP CPX_PARAM_PDSWITCH },
	{ "perturb",	sf_int1,	VP CPX_PARAM_PERIND },
	{ "perturbation", sf_dbl,	VP CPX_PARAM_EPPER },
	{ "perturbconst", sf_dbl,	VP CPX_PARAM_EPPER },
	{ "perturblim",	sf_int,		VP CPX_PARAM_PERLIM },
	{ "perturblimit", sf_int,	VP CPX_PARAM_PERLIM },
	{ "pgradient",	sf_int,		VP CPX_PARAM_PPRIIND },
#ifdef CPLEX_MIP
	{ "plconpri",	sf_mint,	VP set_conpri },
	{ "plobjpri",	sf_mint,	VP set_objpri },
#ifdef CPX_PARAM_POLISHAFTEREPAGAP
	{ "polishafter_absmipgap", sf_dbl, VP CPX_PARAM_POLISHAFTEREPAGAP },
	{ "polishafter_intsol",	sf_int, VP CPX_PARAM_POLISHAFTERINTSOL },
	{ "polishafter_mipgap",	sf_dbl, VP CPX_PARAM_POLISHAFTEREPGAP },
	{ "polishafter_nodes",	sf_int, VP CPX_PARAM_POLISHAFTERNODE },
	{ "polishafter_time",	sf_dbl, VP CPX_PARAM_POLISHAFTERTIME },
#endif
#ifdef  CPX_PARAM_POLISHAFTERDETTIME
	{ "polishafter_timedet", sf_dbl, VP CPX_PARAM_POLISHAFTERDETTIME },
#endif
#ifdef CPX_PARAM_POLISHTIME
	{ "polishtime",	sf_dbl,		VP CPX_PARAM_POLISHTIME },
#endif
#ifdef CPX_PARAM_POPULATELIM
	{ "poolagap",	sf_dbl,		VP CPX_PARAM_SOLNPOOLAGAP },
	{ "poolcapacity", sf_int2,	VP CPX_PARAM_SOLNPOOLCAPACITY },
	{ "pooldual",	sf_mint,	VP set_pooldual },
	{ "poolgap",	sf_dbl,		VP CPX_PARAM_SOLNPOOLGAP },
	{ "poolintensity", sf_int,	VP CPX_PARAM_SOLNPOOLINTENSITY },
	{ "poolreplace", sf_int,	VP CPX_PARAM_SOLNPOOLREPLACE },
	{ "poolstub",	sf_char,	VP set_poolstub },
	{ "populate",	sf_mint,	VP set_populate },
	{ "populatelim", sf_int,	VP CPX_PARAM_POPULATELIM },
#endif
#endif /*CPLEX_MIP*/
#ifdef CPX_PARAM_PRECOMPRESS
	{ "precompress", sf_int,	VP CPX_PARAM_PRECOMPRESS },
#endif
	{ "predual",	sf_int,		VP CPX_PARAM_PREDUAL },
	{ "prelinear",	sf_int,		VP CPX_PARAM_PRELINEAR },
#ifdef CPX_PARAM_PREPASS
	{ "prepass",	sf_int,		VP CPX_PARAM_PREPASS },
#endif
	{ "prereduce",	sf_int,		VP CPX_PARAM_REDUCE },
#ifdef CPLEX_MIP
	{ "prerelax",	sf_int1,	VP CPX_PARAM_RELAXPREIND },
#endif
	{ "presolve",	sf_int1,	VP CPX_PARAM_PREIND },
	{ "presolvedual",sf_int,	VP CPX_PARAM_PREDUAL },
#ifdef CPLEX_MIP
	{ "presolvenode",sf_int,	VP CPX_PARAM_PRESLVND },
#endif
	{ "prestats",	sf_mint,	VP set_prestats },
#ifdef CPX_TUNE_TILIM
	{ "pretunefile", sf_char,	VP set_pretunefile },
	{ "pretunefileprm", sf_char,	VP set_pretunefileprm },
#endif
	{ "pricing",	sf_int,		VP CPX_PARAM_PRICELIM },
	{ "primal",	sf_known,	VP set_primal },
	{ "primalopt",	sf_known,	VP set_primalopt },
#ifdef CPLEX_MIP
	{ "priorities",	sf_int1,	VP CPX_PARAM_MIPORDIND },
	{ "probe",	sf_int,		VP CPX_PARAM_PROBE },
#ifdef CPX_PARAM_PROBETIME
	{ "probetime", sf_dbl,		VP CPX_PARAM_PROBETIME },
#endif
#ifdef CPX_PARAM_PROBEDETTIME
	{ "probetimedet", sf_dbl,	VP CPX_PARAM_PROBEDETTIME },
#endif
#endif /*CPLEX_MIP*/
#ifdef CPX_PARAM_BARQCPEPCOMP
	{ "qcpconvergetol", sf_dbl,	VP CPX_PARAM_BARQCPEPCOMP },
#endif
	{ "rays",	sf_mint,	VP set_rays },
	{ "readbasis",	sf_char,	VP set_startbas },
	{ "readsol",	sf_char,	VP set_startsol },
#if defined(BARRIER) && !defined(NO_DEPRECATED)
	{ "readvector",	sf_char,	VP set_startvector },
#endif
	{ "refactor",	sf_int,		VP CPX_PARAM_REINV },
#ifdef CPLEX_MIP
	{ "relax",	sf_known,	VP set_relax },
	{ "relaxpresolve", sf_int,	VP CPX_PARAM_RELAXPREIND },
	{ "relobjdif",	sf_dbl2,	VP CPX_PARAM_RELOBJDIF },
	{ "relobjdiff",	sf_dbl,		VP CPX_PARAM_RELOBJDIF },
	{ "relpresolve", sf_int,	VP CPX_PARAM_RELAXPREIND },
#ifdef CPX_PARAM_REPAIRTRIES
	{ "repairtries", sf_int,	VP CPX_PARAM_REPAIRTRIES },
#endif
#ifdef CPX_PARAM_REPEATPRESOLVE
	{ "repeatpresolve", sf_int,	VP CPX_PARAM_REPEATPRESOLVE },
#endif
#if CPX_VERSION >= 12030000
	{ "reqconvex",	sf_int,		VP CPX_PARAM_SOLUTIONTARGET },
#endif
	{ "resolve",	sf_mint,	VP set_resolve },
	{ "return_mipgap", sf_mint,	VP set_retmipgap },
#ifdef CPX_PARAM_RINSHEUR
	{ "rinsheur",	sf_int,		VP CPX_PARAM_RINSHEUR },
#endif
#ifdef CPX_PARAM_HEURISTIC
	{ "rootheuristic", sf_int,	VP CPX_PARAM_HEURISTIC },
#endif
	{ "round",	sf_mint,	VP set_round },
#endif /*CPLEX_MIP*/
	{ "scale",	sf_int,		VP CPX_PARAM_SCAIND },
#ifdef CPX_PARAM_RANDOMSEED
	{ "seed", sf_int,		VP CPX_PARAM_RANDOMSEED },
#endif
	{ "sensitivity", sf_known,	VP set_sens },
	{ "siftingopt",	sf_known,	VP set_siftopt },
	{ "siftopt",	sf_known,	VP set_siftopt },
#ifdef CPX_PARAM_SIMTHREADS
	{ "simthreads",	sf_int,		VP CPX_PARAM_SIMTHREADS },
#endif
	{ "singular",	sf_int,		VP CPX_PARAM_SINGLIM },
	{ "singularlim", sf_int,	VP CPX_PARAM_SINGLIM },
#ifdef CPLEX_MIP
	{ "solutionlim", sf_int,	VP CPX_PARAM_INTSOLLIM },
	{ "sos",	sf_mint,	VP set_sos },
	{ "sos2",	sf_mint,	VP set_sos2 },
#ifdef CPX_PARAM_LANDPCUTS
	{ "splitcuts",	sf_int,		VP CPX_PARAM_LANDPCUTS },
#endif
	{ "startalg",	sf_int,		VP CPX_PARAM_STARTALG },
	{ "startalgorithm", sf_int,	VP CPX_PARAM_STARTALG },
	{ "startbasis",	sf_char,	VP set_startbas },
	{ "startsol",	sf_char,	VP set_startsol },
#ifdef CPX_PARAM_IISIND /* version < 9.2b */
	{ "starttree",	sf_char,	VP set_starttree },
#endif
#if defined(BARRIER) && !defined(NO_DEPRECATED)
	{ "startvector", sf_char,	VP set_startvector },
#endif
#ifdef CPLEX_MIP
	{ "strongcand",	sf_int,		VP CPX_PARAM_STRONGCANDLIM},
	{ "strongit",	sf_int,		VP CPX_PARAM_STRONGITLIM},
#ifdef CPX_PARAM_STRONGTHREADLIM
	{ "strongthreads", sf_int,	VP CPX_PARAM_STRONGTHREADLIM},
#endif
#endif
	{ "subalg",	sf_int,		VP CPX_PARAM_SUBALG },
	{ "subalgorithm", sf_int,	VP CPX_PARAM_SUBALG },
#endif
#ifdef CPX_PARAM_SUBMIPNODELIM
	{ "submipnodelim", sf_int,	VP CPX_PARAM_SUBMIPNODELIM },
#endif
#ifdef CPLEX_MIP
#ifdef CPX_PARAM_SYMMETRY
	{ "symmetry",	sf_int,		VP CPX_PARAM_SYMMETRY },
#endif
#endif
#ifdef CPX_PARAM_THREADS
	{ "threads",	sf_int,		VP CPX_PARAM_THREADS },
#endif
	{ "time",	sf_dbl,		VP CPX_PARAM_TILIM },
	{ "timelimit",	sf_dbl,		VP CPX_PARAM_TILIM },
	{ "timing",	sf_mint, 	VP set_timing },
	{ "tranopt",	sf_known,	VP set_dualopt },
#ifdef CPLEX_MIP
	{ "treelimit",	sf_dbl2,	VP CPX_PARAM_TRELIM },
	{ "treememlim",	sf_dbl2,	VP CPX_PARAM_TRELIM },
	{ "treememory",	sf_dbl,		VP CPX_PARAM_TRELIM },
#endif
#ifdef CPX_TUNE_TILIM
	{ "tunedisplay",sf_int,		VP CPX_PARAM_TUNINGDISPLAY },
	{ "tunefile",	sf_char,	VP set_tunefile },
	{ "tunefileprm",sf_char,	VP set_tunefileprm },
	{ "tunefix",	sf_char,	VP set_tunefix },
	{ "tunefixfile",sf_char,	VP set_tunefixfile },
	{ "tunerepeat",	sf_int,		VP CPX_PARAM_TUNINGREPEAT },
	{ "tunetime",	sf_dbl,		VP CPX_PARAM_TUNINGTILIM },
#ifdef CPX_PARAM_TUNINGDETTILIM
	{ "tunetimedet", sf_dbl,	VP CPX_PARAM_TUNINGDETTILIM },
#endif
#endif
#ifdef CPLEX_MIP
	{ "uppercutoff", sf_dbl,	VP CPX_PARAM_CUTUP },
#endif
	{ "upperobj",	sf_dbl,		VP CPX_PARAM_OBJULIM },
	{ "upperobjlim",sf_dbl,		VP CPX_PARAM_OBJULIM },
#ifdef CPLEX_MIP
	{ "varsel",	sf_int,		VP CPX_PARAM_VARSEL },
	{ "varselect",	sf_int,		VP CPX_PARAM_VARSEL },
#endif
	{ "version",	Ver_val,	0 },
	{ "wantsol",	WS_val,		0 },
	{ "workfiledir",sf_char,	VP set_workfiledir },
#ifdef CPX_PARAM_WORKMEM
	{ "workfilelim",sf_dbl,		VP CPX_PARAM_WORKMEM },
#endif
	{ "writebasis",	sf_char,	VP set_endbas },
	{ "writeprob",	sf_char,	VP set_wrtfname },
	{ "writesol",	sf_char,	VP set_endsol },
#if defined(BARRIER) && !defined(NO_DEPRECATED)
	{ "writevector", sf_char,	VP set_endvector },
#endif
	{ "xxxstart",	sf_int,		VP CPX_PARAM_XXXIND }
#ifdef CPX_PARAM_ZEROHALFCUTS
	,{ "zerohalfcuts", sf_int,	VP CPX_PARAM_ZEROHALFCUTS }
#endif
	};

 static Option_Info Oinfo = { "cplex", 0, "cplex_options",
				keywds, nkeywds, 0, cplex_version,
				0,0,0,0,0, 20130531 };

 static void
badlic(int rc, int status)
{
	char buf[4096];

	if (!CPXgeterrorstring(Env, status, buf))
		Sprintf(buf,
		"CPLEX licensing problem: error code %d from CPXopenCPLEX.",
			status);
	badretfmt(rc, "%s", buf);
	}

 static void
nonlin(int n, int rc, char *what)
{
	if (n) {
		badretfmt(rc, "%s contains %s.\n", filename, what);
		exit(4);
		}
	}

 static unsigned
namemem(int n)
{
	int k, L;
	unsigned int rv;

	L = 3;
	k = 9;
	rv = 0;
	while(n > k) {
		rv += k*L++;
		n -= k;
		k *= 10;
		}
	return rv + n*L;
	}

 static char *
namegen(char *stub, char **name, char *s, int n)
{
	int i, j, k, k0, L;
	k0 = 10;
	for(i = 1; i <= n; i++) {
		if (k0 <= i)
			k0 *= 10;
		*name++ = s;
		*s++ = *stub;
		for(j = i, k = k0; k /= 10; ) {
			L = j / k;
			j -= L * k;
			*s++ = '0' + L;
			}
		*s++ = 0;
		}
	return s;
	}

 static unsigned int
brlen(char **name, int n)
{
	char *s, *s0;
	unsigned int x = 0;

 top:
	while(n > 0)
		for(s = s0 = name[--n];;)
			switch(*s++) {
				case 0: goto top;
				case '[': case ']':
					x += strlen(s0) + 1;
					goto top;
				}
	return x;
	}

 static char*
br_round(char **name, int n, char *s)
{
	char *t, *t0;

 top:
	while(n > 0)
		for(t = t0 = name[--n];;)
			switch(*t++) {
				case 0: goto top;
				case '[': case ']':
					name[n] = s;
					for(;;) switch(*s++ = *t0++) {
						case 0: goto top;
						case '[': s[-1] = '('; break;
						case ']': s[-1] = ')';
						}
				}
	return s;
	}

 static void
make_names(char ***cnamep, char ***rnamep)
{
	char buf[24], **cname, **rname, *s;
	int cnlen, i, j, n, rnlen;
	unsigned int cstorsz, rstorsz, rsz;

	cnlen = rnlen = 0;
	if (method > 0) {
		cnlen = maxcolnamelen;
		rnlen = maxrownamelen;
		}
	cstorsz = cnlen > 0 ? 0 : namemem(nbas);
	rstorsz = rnlen > 0 ? 0 : namemem(mbas);
	*cnamep = cname = (char **)Malloc((mbas+nbas)*sizeof(char *)
				+ cstorsz + rstorsz);
	*rnamep = rname = cname + nbas;
	s = (char*)(rname + mbas);
	if (cnlen > 0) {
		for(i = 0, n = n_var; i < n; i++)
			cname[i] = var_name(i);
		while(i < nbas) {
			j = sprintf(buf, "_svar[%d]", i+1) + 1;
			memcpy(cname[i++] = (char*)mem(j), buf, j);
			}
		}
	else
		s = namegen("x", cname, s, nbas);
	if (rnlen > 0) {
		for(i = 0, n = n_con; i < n; i++)
			rname[i] = con_name(i);
		while(i < mbas) {
			j = sprintf(buf, "_scon[%d]", i+1) + 1;
			memcpy(rname[i++] = (char*)mem(j), buf, j);
			}
		}
	else
		namegen("c", rname, s, mbas);
	if (Nameround) {
		/* This would be better done by running "tr" separately. */
		rsz = 0;
		if (cnlen > 0)
			rsz = brlen(cname, nbas);
		if (rnlen > 0)
			rsz += brlen(rname, mbas);
		if (!rsz)
			return;
		s = (char*)M1alloc(rsz);
		if (cnlen > 0)
			s = br_round(cname, nbas, s);
		if (rnlen > 0)
			br_round(rname, mbas, s);
		}
	}

#ifdef CPLEX_MIP /*{*/
 static int
nzeros(int *p, int n)
{
	int *pe = p + n;
	n = 0;
	while(p < pe)
		if (*p++)
			n++;
	return n;
	}

 static void
mip_priorities(ASL *asl, cpxlp *cpx)
{
	int baddir, badpri, i, j, k, listsize, nnames, pvk, sk;
	int *ci, *colindex, *dir, *direction, *p, *priority;
	int *start, *num, *pri;
	SufDesc *dp, *dd;

	baddir = badpri = i = listsize = 0;
	dd = suf_get("direction", ASL_Sufkind_var);
	dir = dd->u.i;
	dp = suf_get("priority", ASL_Sufkind_var);
	direction = 0;
	if ((p = dp->u.i) || dir) {
		nnames = n_var;
		k = 2;
		if (p && dir) {
			k = 3;
			for(; i < nnames; i++)
				if (p[i] || dir[i])
					listsize++;
			}
		else if (p)
			listsize = nzeros(p, nnames);
		else
			listsize = nzeros(dir, nnames);
		if (!listsize)
			return;
		priority = 0;
		colindex = (int*)Malloc(k*listsize*sizeof(int));
		i = k = 0;
		if (dir && p) {
			priority = colindex + listsize;
			direction = priority + listsize;
			for(; i < nnames; i++)
				if (p[i] || dir[i]) {
					colindex[k] = i;
					if ((j = p[i]) < 0) {
						badpri++;
						j = 0;
						}
					priority[k] = j;
					switch(dir[i]) {
					  case -1:
						j = CPX_BRANCH_DOWN;
						break;
					  case 1:
						j= CPX_BRANCH_UP;
						break;
					  default:
						baddir++;
						/* no break */
					  case 0:
						j = CPX_BRANCH_GLOBAL;
					  }
					direction[k++] = j;
					}
			}
		else if (p) {
			priority = colindex + listsize;
			for(; i < nnames; i++)
				if (p[i]) {
					colindex[k] = i;
					if ((j = p[i]) < 0) {
						badpri++;
						j = 0;
						}
					priority[k++] = j;
					}
			}
		else {
			direction =  colindex + listsize;
			for(; i < nnames; i++)
				if (dir[i]) {
					switch(dir[i]) {
					  case -1:
						j = CPX_BRANCH_DOWN;
						break;
					  case 1:
						j= CPX_BRANCH_UP;
						break;
					  default:
						baddir++;
						/* no break */
					  case 0:
						j = CPX_BRANCH_GLOBAL;
					  }
					direction[k] = j;
					colindex[k++] = i;
					}
			}
		listsize = k;
		if (baddir)
			fprintf(Stderr,
		 "Treating %d .direction values outside [-1, 1] as 0.\n",
				baddir);
		if (badpri)
			fprintf(Stderr,
				"Treating %d negative .priority values as 0\n",
				badpri);
		}
	else if ((nnames = mip_pri(&start, &num, &pri, 2147483647))) {
		for(; i < nnames; i++)
			listsize += num[i];
		ci = colindex = (int *)Malloc(listsize*(2*sizeof(int)));
		priority = p = ci + listsize;
		for(k = 0; k < nnames; k++) {
			i = num[k];
			pvk = pri[k];
			sk = start[k];
			while(--i >= 0) {
				*ci++ = sk++;
				*p++ = pvk;
				}
			}
		listsize = ci - colindex;
		}
	else
		return;
	if (CPXcopyorder(Env,cpx,listsize,colindex,priority,direction))
		printf("error in CPXcopyorder!\n");
	free(colindex);
	}

 static void
mipinit_loop(ASL *asl, int **rp, int *np, real **xp, int j, int k)
{
	int i, *r;
	real *L, *U, t, *v, *x;

	v = X0 + j;
	L = LUv + j;
	U = Uvx + j;
	r = *rp;
	x = *xp;

	for(i = 0; i < k; ++i) {
		*r++ = i + j;
		if ((t = v[i]) < L[i])
			t = L[i];
		else if (t > U[i])
			t = U[i];
		*x++ = t;
		}
	*rp = r;
	*xp = x;
	*np -= k;
	}

 static void
set_mipinit(ASL *asl, cpxlp *cpx, int nint)
{
	int i, n1, *r, *r1;
	real *x, *x1;
#ifndef NO_DEPRECATED
	int k, m, m0;
#endif
#if CPX_VERSION_VERSION >= 12
	int beg[2];

	if (mipstval == 1)
		nint = n_var;
#endif
	x = x1 = (real*)Malloc(nint*(sizeof(real) + sizeof(int)));
	r = r1 = (int*)(x + nint);
	n1 = nint;
#if CPX_VERSION_VERSION >= 12
	if (
#ifndef NO_DEPRECATED
		mipstval ==
#endif
			    1)
		{
		mipinit_loop(asl, &r1, &n1, &x1, 0, nint);
		beg[0] = 0;
		beg[1] = nint;
		if ((i = CPXaddmipstarts(Env, cpx, 1, nint, beg, r, x, 0, 0)))
			badret("CPXaddmipstarts", i, 0);
		}
	else {
#endif
#ifndef NO_DEPRECATED
	if ((k = nlvbi))
		mipinit_loop(asl, &r1, &n1, &x1, nlvb - k, k);
	m0 = nlvb;
	m = nlvc - nlvb;
	if ((k = nlvci))
		mipinit_loop(asl, &r1, &n1, &x1, m0 + m - k, k);
	m0 += m;
	if ((k = nlvoi))
		mipinit_loop(asl, &r1, &n1, &x1, m0 + nlvo - nlvc - k, k);
	if (n1)
		mipinit_loop(asl, &r1, &n1, &x1, n_var - (niv + nbv), n1);
	if ((i = CPXcopymipstart(Env, cpx, nint, r, x)))
		badret("CPXcopymipstart", i, 0);
#ifdef CPX_PARAM_MIPSTART
	else if (i =  CPXsetintparam(Env, CPX_PARAM_MIPSTART, 1))
		badret("CPXsetintparam(CPX_PARAM_MIPSTART)", i, 0);
#endif
#endif
#if CPX_VERSION_VERSION >= 12
		}
#endif
	free(x);
	}
#endif /*}*/

 static int
parval(int par)
{
	int rv;
	CPXgetintparam(Env, par, &rv);
	return rv;
	}

 static int CPXPUBLIC
lpcbf(CPXCENVptr env, void *lp, int wf, void *cbh)
{
	int i;
	static int n[4];
	static int wi[8] = {
		CPX_CALLBACK_INFO_PRESOLVE_ROWSGONE,
		CPX_CALLBACK_INFO_PRESOLVE_COLSGONE,
		CPX_CALLBACK_INFO_PRESOLVE_AGGSUBST,
		CPX_CALLBACK_INFO_PRESOLVE_COEFFS,
		CPX_CALLBACK_INFO_CROSSOVER_PPUSH,
		CPX_CALLBACK_INFO_CROSSOVER_PEXCH,
		CPX_CALLBACK_INFO_CROSSOVER_DPUSH,
		CPX_CALLBACK_INFO_CROSSOVER_DEXCH
		};
	static int wf0 = -1;
	static void *lp0;
	Not_Used(cbh);

	if ((wf != wf0 || lp != lp0) && wf0 == CPX_CALLBACK_PRESOLVE && prestats) {
		if (n[0] + n[1] + n[2])
			for(i = 0; i < 4; ++i) {
				cbi.np[i] += n[i];
				n[i] = 0;
				}
		n[3] = 0;
		}
	switch(wf) {
	  case CPX_CALLBACK_PRESOLVE:
		if (prestats)
			for(i = 0; i < 4; i++)
				CPXgetcallbackinfo(env, lp, wf, wi[i], &n[i]);
		break;
	  case CPX_CALLBACK_PRIMAL_CROSSOVER:
		cbi.xkind = "primal";
		for(i = 0; i < 2; i++)
			CPXgetcallbackinfo(env, lp, wf, wi[i+4], &cbi.nx[i]);
		break;
	  case CPX_CALLBACK_DUAL_CROSSOVER:
		cbi.xkind = "dual";
		for(i = 0; i < 2; i++)
			CPXgetcallbackinfo(env, lp, wf, wi[i+6], &cbi.nx[i]);
		break;
	  }
	lp0 = lp;
	wf0 = wf;
	return breaking > 1;
	}

 static void
stat_map(int *stat, int n, int *map, int mx, char *what)
{
	int bad, i, i1, j, j1;
	static char badfmt[] = "cplex driver: %s[%d] = %d\n";

	for(i = bad = 0; i < n; i++) {
		if ((j = stat[i]) >= 0 && j <= mx)
			stat[i] = map[j];
		else {
			stat[i] = 0;
			i1 = i;
			j1 = j;
			if (!bad++)
				fprintf(Stderr, badfmt, what, i, j);
			}
		}
	if (bad > 1) {
		if (bad == 2)
			fprintf(Stderr, badfmt, what, i1, j1);
		else
			fprintf(Stderr,
		"cplex driver: %d messages about bad %s values suppressed.\n",
				bad-1, what);
		}
	}

 static void
get_statuses(ASL *asl, cpxlp *cpx, dims *d)
{
	SufDesc *sd;
	int *cs, *cstat, havestats, i, j, m, n, *rs, *rstat;
	real *L, *U;
	static int map[] = {0, 1, 3, 0, 2, 3, 3};

	havestats = 0;
	sd = d->csd;
	n = n_var;
	if (sd->kind & ASL_Sufkind_input)
		havestats = 1;
	sd = d->rsd;
	m = n_con;
	if (sd->kind & ASL_Sufkind_input)
		havestats = 1;
	d->cs = d->rs = 0;
	if (!havestats || (!mipststat && niv + nbv))
		return;
	if (method > 0) {
		stat_map(d->cs = d->cstat, n, map, 6, "incoming cstat");
		stat_map(d->rs = d->rstat, m, map, 6, "incoming rstat");
		if (objadj) {
			d->cstat[n] = CPX_BASIC;
			d->rstat[m] = CPX_FREE_SUPER;
			}
		if ((i = CPXcopybase(Env, cpx, d->cstat, d->rstat))) {
 bad_copybase:
			badret("CPXcopybase", i, 0);
			}
		return;
		}
	cs = (int*)(use_netopt ? M1alloc((mbas + nbas)*sizeof(int))
				: Malloc((mbas + nbas)*sizeof(int)));
	rs = cs + nbas;
	m = n_con;
	n = n_var;
	L = LUrhs;
	U = Urhsx;
	cstat = d->cstat;
	rstat = d->rstat;
	for(i = 0; i < m; i++) {
		if (L[i] <= negInfinity)
			cs[i] = rstat[i] == 1 ? CPX_AT_UPPER : CPX_BASIC;
		else
			cs[i] = rstat[i] == 1 ? CPX_AT_LOWER : CPX_BASIC;
		}
	if (nranges) {
		for(j = 0; j < m; j++)
			if (L[j] < U[j]
			 && L[j] > negInfinity
			 && U[j] < Infinity)
				cs[i++] = rstat[j] == 4
					? CPX_BASIC : CPX_AT_UPPER;
		}
	L = LUv;
	U = Uvx;
	for(j = 0; j < n; j++) {
		rs[j] = CPX_AT_LOWER;
		if (L[j] > negInfinity) {
			if (L[j])
				cs[i++] = cstat[j] == 3
					? CPX_BASIC : CPX_AT_LOWER;
			else
				rs[j] = cstat[j] == 3 ?
					CPX_BASIC : CPX_AT_UPPER;
			if (U[j] < Infinity)
				goto finite_ub;
			}
		else if (U[j] < Infinity) {
 finite_ub:
			if (U[j])
				cs[i++] = cstat[j] == 4
					? CPX_BASIC : CPX_AT_LOWER;
			else
				rs[j] = cstat[j] == 4
					? CPX_BASIC : CPX_AT_LOWER;
			}
		else if (cstat[j] != 1)
			rs[j] = CPX_BASIC;
		}
	if (objadj)
		cs[i] = 1;
	i = CPXcopybase(Env, cpx, cs, rs);
	if (use_netopt) {
		d->cs = cs;
		d->rs = rs;
		}
	else
		free(cs);
	if (i)
		goto bad_copybase;
	}

 static int
qmatadj(int k, int nr, int os, int *colq, int *colqcnt, double **qmatp)
{
	double badd[3], *oqmat, osd, *qmat, t;
	int badk[3], i, nbad = 0, nb, w, w1;

	oqmat = *qmatp;
	*qmatp = qmat = (double *)Malloc(nr*sizeof(double));
	while(nr > k)
		qmat[--nr] = 0.;
	osd = os;
	while(--k >= 0) {
		t = 0.;
		if (colqcnt[k]) {
			t = oqmat[colq[k]];
			if (t*osd < 0.) {
				if (nbad < 3) {
					badd[nbad] = t;
					badk[nbad] = k;
					}
				nbad++;
				}
			}
		qmat[k] = t;
		}
	free(oqmat);
	if (!nbad)
		return 0;
	fprintf(Stderr, "%d diagonal QP coefficients of the wrong sign:\n",
		nbad);
	if ((nb = nbad) > 3)
		nb = 3;
	w = 8;
	for(i = 0; i < nb; i++) {
		w1 = strlen(var_name(badk[i]));
		if (w < w1)
			w = w1;
		}
	fprintf(Stderr, "Variable%.*s Diagonal\n", w-8, "");
	for(i = 0; i < nb; i++)
		fprintf(Stderr, "%-*s %g\n", w, var_name(badk[i]),
			0.5*badd[i]);
	solve_result_num = 540;
	asl->i.uinfo = "Diagonal QP Hessian has elements of the wrong sign.";
	return 1;
	}

 static int
refcomp(const void *a, const void *b, void *c)
{
	double d, *x;

	x = (double*)c;
	d = x[*(int*)a] - x[*(int*)b];
	if (d < 0.)
		return -1;
	if (d > 0.)
		return 1.;
	return 0;
	}

 static void
sos_kludge(int nsos, int *sosbeg, double *sosref)
{
	/* Adjust sosref if necessary to accommodate CPLEX's */
	/* undocumented requirement that sosref values differ */
	/* by at least 1e-10. */
	int i, i0, i1, j, k, m, m1, *z, z0[16];
	double t, t1;

	for(i = j = 0; i++ < nsos; ) {
		k = sosbeg[i];
		t = sosref[j];
		while(++j < k) {
			t1 = sosref[j];
			t += 1e-10;
			if (t1 <= t) {
				if (t1 < t)
					goto trysort;
				sosref[j] = t1 = t + 1e-10;
				}
			t = t1;
			}
		}
	return;
 trysort:
	j = sosbeg[i0 = i - 1];
	m = 0;
	for(i = i0; i++ < nsos; j = k) {
		k = sosbeg[i];
		if ((m1 = k - j) > m)
			m = m1;
		}
	z = z0;
	if (m > sizeof(z0)/sizeof(z0[0]))
		z = (int *)Malloc(m*sizeof(int));
	for(j = sosbeg[i = i0]; i++ < nsos; j = k) {
		k = sosbeg[i];
		m1 = k - j;
		for(i1 = 0; i1 < m1; ++i1)
			z[i1] = i1 + j;
		qsortv(z, m1, sizeof(int), refcomp, sosref);
		t = sosref[z[0] + j];
		for(i1 = 1; i1 < m1; ++i1, t = t1) {
			t1 = sosref[z[i1]];
			t += 1e-10;
			if (t1 <= t)
				sosref[z[i1]] = t1 = t + 1e-10;
			}
		}
	if (z != z0)
		free(z);
	}

#ifdef CPXERR_QCP_SENSE

 static void
Surprise(ASL *asl, int n, char *who)
{
	char msgbuf[64];

	sprintf(msgbuf, "Surprise return %d from %s.", n, who);
	nonlin(1, 561, msgbuf);
	}

 static Char**
linadj(ASL *asl, int *ka, int *ia, double *a, int **kaqp, int **kalqp, int **iaqp, double **aqp,
	int **nelqcp, int ***rowqcp, int ***colqcp, double ***qcmatp)
{
	/* Move linear parts of nqc quadratic constraints from ka, ia, a */
	/* to kaq, kalq, iaq, aq and corresponding Cgrad pointers and return a */
	/* pointer to the relevant allocation. */

	cgrad **cgp, *cg, *cga, *cg0;
	double *aq, **qcmat, *x, *y;
	fint *colqf, *cqf, *rowqf;
	int i, i1, i2, j, k, nl, nqc, m, m0, n, nelq, nqcnl;
	int **colqc, *iaq, *kalq, *kaq, *kaq0, *nelqc, **rowqc, *z, *z1;
	size_t L;

	n = n_var;
	nqc = nlc;

	j = nzc;
	for(i = nl = 0; i < j; i++)
		if (ia[i] < nqc)
			nl++;

	cg = cga = (cgrad*)Malloc(nqc*(sizeof(int)+sizeof(cgrad*)) + nl*sizeof(cgrad));
	Cgrad = cgp = (cgrad**)(cga + nl);
	kaq0 = (int*)(cgp + nqc);
	L = (char*)kaq0 - (char*)cga;
	memset(kaq0, 0, nqc*sizeof(int));
	memset(cgp, 0, nqc*sizeof(cgrad*));

	for(i = j = 0; i < n; i++) {
		for(k = ka[i+1]; j < k; j++)
			if ((i1 = ia[j]) < nqc) {
				kaq0[i1]++;
				cg->varno = i;
				cg->coef = a[j];
				cg->next = cgp[i1];
				cgp[i1] = cg++;
				}
		}
	x = LUrhs;
	y = Urhsx;
	for(i = nl = nqcnl = 0; i < nqc; i++) {
		if (x[i] > negInfinity && y[i] < Infinity) {
			badretfmt(562,
			 "Constraint %s is not convex quadratic since it is %s constraint.",
				con_name(i), x[i] == y[i] ? "an equality"
						: "a two-sided");
			exit(4);
			}
		k = nl;
		nl += kaq0[i];
		kaq0[i] = k;
		j = mqpcheck(-(i+1), 0, 0, 0);
		if (j < 0) {
			nonlin(j == -2, 558,
			 "a quadratic constraint involving division by 0");
			nonlin(1, 550, "a nonquadratic nonlinear constraint");
			}
		nonlin(j == 0, 559,
		 "CPLEX driver bug: no quadratic terms in \"nonlinear\" constraint");
		nqcnl += j;
		}

	*aqp = aq = (double*)Realloc(cga, nl*sizeof(cgrad) + nqc*sizeof(cgrad*)
				+ nqc*(sizeof(double*) + 2*sizeof(int*))
				+ (3*nqc + 2*nqcnl + nl)*sizeof(int)
				+ (nl+nqcnl)*sizeof(double));
	x = aq + nl;
	cga = (cgrad*)(x + nqcnl);
	Cgrad = cgp = (cgrad**)(cga + nl);
	*qcmatp = qcmat = (double**)(cgp + nqc);
	*rowqcp = rowqc = (int**)(qcmat + nqc);
	*colqcp = colqc = rowqc + nqc;
	*kaqp = kaq = (int*)(colqc + nqc);
	*kalqp = kalq = kaq + nqc;
	*nelqcp = nelqc = kalq + nqc;
	z = nelqc + nqc;
	*iaqp = iaq = z + 2*nqcnl;

	kaq0 = (int*)((char*)aq + L);
	if (kaq0 < kaq)
		for(i = nqc; i-- > 0; )
			kaq[i] = kaq0[i];
	else
		for(i = 0; i < nqc; i++)
			kaq[i] = kaq0[i];

	memset(cgp, 0, nqc*sizeof(cgrad*));

	for(i = j = m = 0; j < n; j++) {
		ka[j] = m0 = m;
		for(k = ka[j+1]; i < k; i++)
			if ((i1 = ia[i]) < nqc) {
				iaq[i2 = kaq[i1]++] = j;
				aq[i2] = a[i];
				}
			else {
				ia[m] = i1 - nqc;
				a[m++] = a[i];
				}
		}
	nzc = ka[n] = m;
	for(i = nqc; --i > 0; )
		kaq[i] = kaq[i-1];
	kaq[0] = 0;
	j = nl;
	for(i = nqc; i > 0; ) {
		k = kaq[--i];
		cg0 = 0;
		while(j > k) {
			cg = cga++;
			cg->next = cg0;
			cg0 = cg;
			cg->coef = aq[--j];
			cg->varno = iaq[j];
			}
		cgp[i] = cg0;
		}
	for(i = m = 0; i < nqc; i++) {
		nelqc[i] = nelq = mqpcheck(-(i+1), &rowqf, &colqf, &y);
		qcmat[i] = x;
		for(j = 0; j < nelq; j++)
			*x++ = 0.5*y[j];
		rowqc[i] = z1 = z;
		colqc[i] = z += nelq;
		cqf = colqf;
		i1 = *++cqf;
		for(j = k = 0; j < nelq; j++) {
			while (j >= i1) {
				k++;
				i1 = *++cqf;
				}
			*z1++ = rowqf[j];
			*z++  = k;
			}
		free(colqf);
		free(rowqf);
		free(y);
		kaq[i] = m0 = m;
		for(cg = cgp[i]; cg; cg = cg->next)
			if (cg->coef) {
				iaq[m] = cg->varno;
				aq[m++] = cg->coef;
				}
		kalq[i] = m - m0;
		}

	return M1record((Char*)aq);
	}

#endif /*CPXERR_QCP_SENSE*/

#ifdef Uselazy
 typedef struct
LazyInfo {
	char *sense[2], **rowname[2];
	int *colno[2], nz[2], *rowbeg[2], rows[2];
	double *matval[2], *rhs[2];
	} LazyInfo;

 Char **
lazyadj(ASL *asl, LazyInfo *LI, int n, int nint, int nqc, int *mp, double *b, char *senx,
	int *ka, int *kal, int *ia, double *a, double *rngvec, char **rname)
{
	Char *rv;
	SufDesc *lp;
	char *nn, *nn1, *sense, **rowname;
	const char *s1, *s2;
	double *matval, *rhs;
	int i, i1, i2, j, j0, j1, k, m, m0, nq, nrt, nz;
	int *colno, nr[3], nt[2], *rl, *rowbeg, *z, *z0;
	size_t Lnn, Lrl;

	nr[0] = nr[1] = nr[2] = 0;
	rv = 0;
	if (!(uselazy & 3) || !nint)
		goto done;
	lp = suf_get("lazy", ASL_Sufkind_con);
	if (!(z = lp->u.i))
		goto done;

	/* Expect .lazy values, i.e., values in z, to be 0, 1, or 2 */
	/* with 1 ==> lazy constraint, 2 ==> user cut. */
	/* For simplicity, we treat 3 as 1 and ignore bits beyond the first two. */

	m0 = n_con;
	for(i = nq = 0; i < m0; i++) {
		if ((j = z[i] & uselazy) && j <= 2) {
			if (i < nqc) {
				++nq;
				continue;
				}
			++nr[j = 1 - (j&1)];
			if (senx[i] == 'R') {
				++nr[j];
				++nr[2];
				}
			}
		}
	if (nq) {
		if (nq == 1)
			s1 = s2 = "";
		else {
			s1 = "es";
			s2 = "s";
			}
		fprintf(Stderr, "Ignoring .lazy suffix%s on %d quadratic constraint%s.\n",
			s1, nq, s2);
		}

	/* nr[0] = number of lazy constraints */
	/* nr[1] = number of user cuts */
	/* nr[2] = number of range constraints among above */

	if (!(nrt = nr[0] + nr[1]))
		goto done;
	z0 = z += nqc;
	m = *mp;
	if (m > m0) {
		z = (int*)Malloc(m*sizeof(int));
		memcpy(z, z0, m0*sizeof(int));
		memset(&z[m0], 0, (m-m0)*sizeof(int));
		}
	*mp = m - (nrt - nr[2]);
	nt[1] = 0;
	nt[0] = nr[0];
	/* nt[0] = for user cuts, nt[1] for lazy constraints */
	if (nr[2])
		for(i = 0; i < m; i++) {
			if ((j = (z[i]&3))) {
				z[i] = ++nt[j &= 1];
				if (senx[i] == 'R')
					++nt[j];
				}
			else
				z[i] = 0;
			}
	else
		for(i = 0; i < m; i++)
			z[i] = (j = (z[i]&3)) ? ++nt[j&1] : 0;
	for(i = nz = 0; i < n; i++) {
		for(j = ka[i], k = j + kal[i]; j < k; j++)
			if (z[ia[j]])
				nz++;
		}
	Lnn = Lrl = 0;
	if (rname)
		Lnn = nrt*sizeof(char*);
	if (nr[2]) { /* ranges */
		Lrl = nrt;
		for(i = 0; i < n; i++) {
			for(j = ka[i], k = j + kal[i]; j < k; j++)
				if (z[j1 = ia[j]] && senx[j1] == 'R')
					nz++;
			}
		if (rname)
			for(i = 0; i < m; i++)
				if (senx[i] == 'R')
					Lnn += strlen(rname[i]) + 3;
		}
	matval = (double*)Malloc((nz+nrt)*(sizeof(double)+sizeof(int))
				+ (Lrl+2)*sizeof(int) + nrt + Lnn);
	rv = (Char*)matval;
	rhs = matval + nz;
	if (rname) {
		rowname = (char**)(rhs + nrt);
		rowbeg = (int*)(rowname + nrt);
		}
	else {
		rowname = 0;
		rowbeg = (int*)(rhs + nrt);
		}
	colno = rowbeg + nrt + 2;
	rl = colno + nz;
	sense = (char*)(rl + Lrl);
	nn = sense + nrt;
	memset(rowbeg, 0, nrt*sizeof(int));
	if (Lrl)
		memset(rl, 0, Lrl*sizeof(int));
	for(i = 0; i < n; i++) {
		for(j = ka[i], k = j + kal[i]; j < k; j++)
			if ((j1 = z[ia[j]]))
				++rowbeg[j1-1];
		}
	if (nr[2]) {
		for(i = 0; i < m; i++)
			if ((j = z[i]) && senx[i] == 'R')
				rl[j] = rowbeg[j] = rowbeg[j-1];
		}
	for(i = k = 0; i < nrt; i++) {
		j = rowbeg[i];
		rowbeg[i] = k;
		k += j;
		}
	for(i = k = 0; i < m; i++)
		if (!z[i])
			z[i] = k--;
	for(i = j1 = 0; i < n; i++) {
		j = ka[i];
		ka[i] = j0 = j1;
		for(k = j + kal[i]; j < k; j++)
			if ((i1 = z[ia[j]]) > 0) {
				i2 = rowbeg[i1-1]++;
				colno[i2] = i;
				matval[i2] = a[j];
				}
			else {
				ia[j1] = -i1;
				a[j1++] = a[j];
				}
		kal[i] = j1 - j0;
		}
	if (nr[2]) {
		for(i = 0; i < nrt; i++)
			if ((j = rl[i])) {
				k = rowbeg[i];
				rowbeg[i] = k + j;
				memcpy(colno + k, colno - j + k, j*sizeof(int));
				memcpy(matval + k, matval - j + k, j*sizeof(double));
				}
		}
	for(i = k = 0; i < nr[0]; i++) {
		j = rowbeg[i];
		rowbeg[i] = k;
		k = j;
		}
	LI->nz[0] = k;
	LI->nz[1] = nz - k;
	for(j1 = 0; i < nrt; i++) {
		j = rowbeg[i];
		rowbeg[i] = j1;
		j1 = j - k;
		}
	rowbeg[i+1] = j1;
	for(; i > nr[0]; --i)
		rowbeg[i] = rowbeg[i-1];
	rowbeg[i] = k;
	LI->sense[0] = sense;
	LI->sense[1] = sense + nr[0];
	for(i = 0; i < m; i++) {
		if ((j = z[i]) > 0) {
			rhs[--j] = b[i];
			if ((sense[j] = senx[i]) == 'R') {
				sense[j] = 'G';
				sense[++j] = 'L';
				rhs[j] = b[i] + rngvec[i];
				if (rname) {
					nn1 = strcpy1(rowname[j] = nn, rname[i]);
					nn = strcpy(nn1, ".u");
					}
				}
			}
		else {
			b[j = -j] = b[i];
			rngvec[j] = rngvec[i];
			senx[j] = senx[i];
			}
		}

	if (rname) {
		LI->rowname[0] = rowname;
		LI->rowname[1] = rowname + nr[0];
		for(i = 0; i < m; i++)
			if ((j = z[i]) > 0)
				rowname[j-1] = rname[i];
			else
				rname[-j] = rname[i];
		}
	else
		LI->rowname[0] = LI->rowname[1]  = 0;
	LI->colno[0] = colno;
	LI->colno[1] = colno + k;
	LI->rowbeg[0] = rowbeg;
	LI->rowbeg[1] = rowbeg + k + 1;
	LI->matval[0] = matval;
	LI->matval[1] = matval + k;
	LI->rhs[0] = rhs;
	LI->rhs[1] = rhs + nr[0];
	if (z != z0)
		free(z);

 done:
	LI->rows[0] = nr[0];
	LI->rows[1] = nr[1];
	return M1record(rv);
	}
#endif /*Uselazy*/

 static void
amplin(ASL *asl, cpxlp **pcpx, FILE **nl, dims *d, int *nelqp, int *nintp, char **av)
{
	int *ia, *ia1, *ia2, *ja, *ka, *ka1, *ka2, *kal, *kal1, *kal2;
	double *a, *a1, *b, *b1, *c, *c1, *l, *l1,
		*le, *lrhs, *lx, *rngvec, *rv1, *u, *u1, *urhs, *ux;
	double nos, os;
	ograd *og, *og1;
	int i, j, k, m, m0, mq, n, n0, nbv1, nextra, nint, nint1, nqc,
		nr, nsos, nz, nz1, nzc0, nzextra, nzr;
	char **cname, **rname, *s, *senx;
	Char **zap[6];
	fint nelqf;
#ifdef BARRIER
	fint *colqf, *rowqf;
	int havestart, nelq, *colq, *colqcnt, *rowq;
	double *cdual, *cprim, *qmat, *rdual, *rprim;
#undef Want_mqpcheck
#ifdef BARRIER_FOR_QP
#ifdef BARRIER
#define Want_mqpcheck
#endif
#else /*!BARRIER_FOR_QP, i.e., >= 8.0*/
#define Want_mqpcheck
#endif
#ifndef Want_mqpcheck
#undef CPXERR_QCP_SENSE
#endif
#ifdef CPXERR_QCP_SENSE	/* CPLEX version >= 9.0 */
	double *aq, **qcmat;
	int **colqc, *iaq, *kaq, *kalq, *nelqc, **rowqc;
#endif /*CPXERR_QCP_SENSE*/
#endif /*BARRIER*/
	char probname[48];
	cpxlp *cpx;
#ifdef CPLEX_MIP /*{*/
	char *ctype;
	char *sostype;
	int nsosnz, *sosbeg, *sosind;
	double *sosref;
#ifdef CPX_PARAM_IISIND /* version < 9.2b */
	int copri[2], **p_sospri, *sospri;
#else
#define copri 0
#define p_sospri 0
#endif
#endif /*}*/
#ifdef CPX_CON_INDICATOR
	int errinfo[2], nlogc;
#else
#define nlogc 0
#endif /*CPX_CON_INDICATOR*/
#ifdef Uselazy
	LazyInfo LI;
#endif /*Uselazy*/

	/* change default display to 0 */
	CPXsetintparam(Env, CPX_PARAM_SIMDISPLAY, 0);
#ifdef CPLEX_MIP
	/* change default mipdisplay to 0 */
	CPXsetintparam(Env, CPX_PARAM_MIPDISPLAY, 0);
#endif
#ifdef BARRIER
	/* change default bardisplay to 0 */
	CPXsetintparam(Env, CPX_PARAM_BARDISPLAY, 0);
#endif
	aggtries = 0;

	if (!(mint_val[set_objno].U = n_obj))
		objno = 0;
	nqc = nlc;
#ifndef CPXERR_QCP_SENSE /* if CPLEX version >= 9.0 */
	nonlin(nqc, 550, "nonlinear constraints");
#endif
	nonlin(plterms, 554, "piecewise-linear terms");
	nonlin(n_cc, 567, "complementarity constraints");

	if (getopts(av, &Oinfo)) {
		if (amplflag)
			badretfmt(560, "Error in $cplex_options.");
		exit(1);
		}
	obj_no = objno - 1;

#if CPX_VERSION >= 12050000
	if (time_flag & 0x30)
		CPXgetnumcores(Env, &num_cores);
#endif
#ifndef CPX_PARAM_FEASOPTMODE /* < 9.2b */
	if (want_iis
#ifdef CPXERR_QCP_SENSE /* if CPLEX version >= 9.0 */
			 && !nqc
#endif
			) {
		if (want_iis == 3)
			want_iis = 2;
		/* The IIS finder requires use of primal or dual simplex. */
		if (Optimize && Optimize != CPXprimopt
#ifdef CPXERR_QCP_SENSE /* if CPLEX version >= 9.0 */
							&& Optimize != CPXdualopt
#endif
			) {
			printf("Assuming primalopt so iisfind will work.\n");
			Optimize = CPXprimopt;
			algname = baralgname = "";
			}
		}
#endif /*CPX_PARAM_FEASOPTMODE*/
	if (!Optimize || nqc) {	/* check nqc in case baropt was specified */
				/* with CPLEX 12.2, at least, we cannot use */
				/* CPXyybbaropt on problems with quadratic constraints. */
		Optimize =
#ifdef CPXERR_QCP_SENSE
			   nqc | nlo ? CPXqpopt :
#endif
						  CPXlpopt;
		lpoptalg = CPX_ALG_AUTOMATIC;
		algname = baralgname = "";
		}

	m = n_con + 1;	/* allow elbow room for obj_adj */
	n = n_var + 1;
	nzr = nzc + 1;
#ifdef CPX_CON_INDICATOR /* >= 9.2b */
	nlogc = 2*n_lcon;
#endif

	d->cstat = (int*)M1zapalloc((m+n+2+nlogc)*sizeof(int));
	d->rstat = d->cstat + n + 1;
	d->csd = suf_iput("sstatus", ASL_Sufkind_var, d->cstat);
	d->rsd = suf_iput("sstatus", ASL_Sufkind_con, d->rstat);

	b = LUrhs = (double *)M1alloc(2*Sizeof(double)*(m+n));
	Urhsx = b + m;
	LUv = Urhsx + m;
	Uvx = LUv + n;
	a = A_vals = (double *)Malloc(nzr*Sizeof(double)
				+ ((fint)nzr + n + 1)*sizeof(int));
	zap[0] = M1record((Char*)a);
	zap[1] = zap[2] = zap[3] = zap[4] = zap[5] = 0;
	ia = A_rownos = (int *)(a + nzr);
	ka = A_colstarts = ia + nzr;
	nelqf = 0;
#ifdef Want_mqpcheck
	nint = nbv + niv + nlvbi + nlvci + nlvoi;
	if (Optimize == Optimize2 || Optimize == Optimizebar || nlo || nqc
	    || (nint && mipstval))
		want_xpi0 = 3;
	want_deriv = 0;
	qp_read(*nl,ALLOW_CLP);
	*nl = 0;	/* Indicate closed file -- in case we're */
			/* subsequently interrupted. */
	if (obj_no >= 0 && obj_no < n_obj
	 && (nelqf = mqpcheck(obj_no, &rowqf, &colqf, &qmat))) {
		if (nelqf < 0) {
			nonlin(nelqf == -2, 555,
			 "a quadratic objective involving division by 0");
			nonlin(1, 551, "a nonlinear objective");
			}
#ifdef BARRIER_FOR_QP
		nonlin(nint, 553,
			"integer variables and a quadratic objective");
#endif
		crossover = 0;
#ifndef BARRIER_FOR_QP
		if (Optimize == Optimize2)
#endif
			set_baropt();
		method = 1;
		}
#ifdef BARRIER_FOR_QP
	else
		nonlin(nlvoi, 552, "nonlinear integer variables");
#endif
	if (mipbasis == -1)
		mipbasis = nelqf == 0;
#ifdef CPXERR_QCP_SENSE /* if CPLEX version >= 9.0 */

	iaq = kaq = kalq = nelqc = 0; /* silence buggy warning */
	colqc = rowqc = 0; /* ditto */
	aq = 0; /* ditto */
	qcmat = 0; /* ditto */
	if (nqc > 0)
		zap[4] = linadj(asl, ka, ia, a, &kaq, &kalq, &iaq, &aq,
				&nelqc, &rowqc, &colqc, &qcmat);

#endif /*CPXERR_QCP_SENSE*/
#else /*Want_mqpcheck*/
	nonlin(nlvoi, 552, "nonlinear integer variables");
	f_read(*nl, 0);
	*nl = 0;
#endif/*Want_mqpcheck*/

	nsos = 0;
#ifdef CPLEX_MIP
	ctype = 0; /* silence buggy warning */
	if (!relax) {
#ifdef CPX_PARAM_IISIND /* version < 9.2b */
		copri[0] = objpri;
		copri[1] = conpri;
		p_sospri = &sospri;
#endif
		i = ASL_suf_sos_explict_free;
		if (!sos)
			i |= ASL_suf_sos_ignore_sosno;
		if (!sos2)
			i |= ASL_suf_sos_ignore_amplsos;
		nsos = suf_sos(i, &nsosnz, &sostype, p_sospri, copri,
				&sosbeg, &sosind, &sosref);
		}
#endif
	m = m0 = n_con;
	n = n0 = nr = n_var;
	nzc0 = nzc;
	nint = nbv + niv + nlvoi + nlvci + nlvbi;
	if (nint
#ifdef CPLEX_MIP
		&& relax
#endif
			  ) {
		if (nint) {
			printf("Ignoring integrality of %d variable%s.\n",
				nint, nint > 1 ? "s" : "");
			need_nl = nint = 0;
			}
		}
	*nelqp = nelqf;
	*nintp = nint1 = nint + nsos;
	if (use_netopt == 1 && !lnc)
		use_netopt = 0;
	if (nint1 || m <= 0 || nqc > 0) {
		method = 1;
		use_netopt = 0;
		}
	else if (!method)
		method = use_netopt ? 1
			: (dual_thresh > 0 && m - n > dual_thresh) || m > dual_ratio*n ? -1 : 1;
	nbv1 = nbv;
	if (method < 0 && use_netopt == 1)
		use_netopt = 0;
	if (!nint)
		niv = nbv = nbv1 = 0;
	og = 0;
	if (obj_no >= 0 && obj_no < n_obj) {
		og = Ograd[obj_no];
		if (!objsen)
			objsen = objtype[obj_no] ? -1 : 1;
		obj_adj = objconst(obj_no);
		}
	else if (!objsen)
		objsen = 1;
#ifndef BARRIER
	nonlin(nl_obj(obj_no), 556, "nonlinear objective terms");
#endif
	*stub_end = 0;
	for(i = 0, og1 = og; og1; i++)
		og1 = og1->next;
	sprintf(probname, "c%dv%di%do%d", n_var, n_con, niv, i);
	if (!(cpx = CPXcreateprob(Env, &i, probname)))
		badret("CPXcreateprob", i, 531);
	*pcpx = cpx;
	n = n0;
	m = m0;
	nextra = nzextra = 0;
	if (method < 0) { /* solve dual */
		l = LUv;
		u = Uvx;
		for(le = l + n0; l < le; l++, u++) {
			if (*l > negInfinity && *l)
				nextra++;
			if (*u < Infinity && *u)
				nextra++;
			}
		if (nranges) {
			ja = (int *)Malloc(m0*sizeof(int));
			memset((char *)ja, 0, m0*sizeof(int));
			ia1 = ia + nzc0;
			while(ia1 > ia)
				ja[*--ia1]++;
			l = LUrhs;
			u = Urhsx;
			for(le = l + m0; l < le; l++, u++)
				if (*l < *u
				 && *l > negInfinity
				 && *u < Infinity)
					nzextra += ja[l - LUrhs];
			free(ja);
			}
		m = n0;
		n = m0 + nranges + nextra;
		}
#ifdef OBJ_ADJ
	if (obj_adj || m <= 0) {
		objadj = 1;
		n++;
		nextra++;
		}
#endif
	nz = nzc0 + nextra;
	if (method < 0) /* dual */ {
		nz1 = nz += nzextra;
		a = (double *)Malloc((nz1+2L*n+m)*sizeof(double)
					+ (nz1+1L+n)*sizeof(int));
		*zap[0] = (Char*)a;
		l = a + nz1;
		u = l + n;
		b = u + n;
		ia = (int *)(b + m);
		ka = ia + nz1;
		nr = n;
		}
	else {
		m += objadj;
#ifdef Student_Edition
		if (n > 300 || m > 300) {
			fflush(stdout);
			fprintf(Stderr,
 "\nSorry, the student edition is limited to 300 variables and 300\n\
constraints.  After adjustments to handle the constant term in\n\
the objective, you have %d variables and %d constraints.\n", n, m);
		exit(1);
		}
#endif
		nr += objadj;
		nz1 = nz;
		l = LUv;
		u = Uvx;
		if (nint1) {
#ifdef CPLEX_MIP
			ctype = M1alloc(nr);
			if ((i = n - nbv1 - niv - nextra))
				memset(ctype, 'C', i);
			if (nbv1)
				memset(ctype+i, 'B', nbv1);
			if (niv)
				memset(ctype+i+nbv1, 'I', niv);
			if (nextra)
				memset(ctype+i+nbv1+niv, 'C', nextra);
#ifndef BARRIER_FOR_QP
			if (nlvbi)
				for(j = nlvb, i = j - nlvbi; i < j; i++)
					ctype[i] = l[i] == 0. && u[i] == 1. ? 'B' : 'I';
			if (nlvci)
				for(j = nlvc, i = j - nlvci; i < j; i++)
					ctype[i] = l[i] == 0. && u[i] == 1. ? 'B' : 'I';
			if (nlvoi)
				for(j = nlvo, i = j - nlvoi; i < j; i++)
					ctype[i] = l[i] == 0. && u[i] == 1. ? 'B' : 'I';
#endif
#else
			badretfmt(557,
		"CPLEX's MIP option is needed to solve problems with\n%s.\n",
			"nonconvex (nonconcave) piecewise-linear terms");
			exit(1);
#endif
			}
		}

	d->x = (double *)M1zapalloc((m+nlogc+2L*nr)*sizeof(double) + m);
	d->y = d->x + nr;
	d->c = c = d->y + m + nlogc;
	d->rtype = senx = (char*)(c + nr);

	rngvec = (double *)Malloc(m*Sizeof(double) + nr*Sizeof(int));
	zap[1] = M1record((Char*)rngvec);
	kal = (int *)(rngvec + m);
#ifdef CPLEX_MIP
	if (nint1) {
		Optimize = CPXmipopt;
		baralgname = "";
		}
#endif
	if (method > 0) {
		b = (double*)Malloc(m*Sizeof(double));
		zap[2] = M1record((Char*)b);
		}
	b1 = b + m;
	s = senx + m;
	rv1 = rngvec + m;
	ia1 = ia + nz;
	a1 = a + nz;
	l1 = l + n;
	u1 = u + n;
	kal1 = kal + n;
	ka1 = ka + n;
	lrhs = LUrhs + m0;
	urhs = Urhsx + m0;
	if (method < 0) { /* dual */
		os = objsen;
		nos = -os;
		objsen = -objsen;
		memset(b, 0, m*sizeof(double));
		c1 = c + (m0 + nranges + nextra);
#ifdef OBJ_ADJ
		if (objadj) {
			*--l1 = *--u1 = 1;
			*--kal1 = 0;
			*--ka1 = nz1;
			*--c1 = obj_adj;
			}
#endif
		j = n0;
		lx = LUv + n0;
		ux = Uvx + n0;
		do {
			*--rv1 = 0;
			*--s = 'E';
			--j;
			if (*--ux < Infinity) {
				if (*ux) {
					*--kal1 = 1;
					*--ka1 = --nz1;
					*--ia1 = j;
					*--a1 = -1;
					*--c1 = nos * *ux;
					*--l1 = 0;
					*--u1 = Infinity;
					}
				else
					*s = 'G';
				}
			if (*--lx > negInfinity) {
				if (*lx) {
					*--kal1 = 1;
					*--ka1 = --nz1;
					*--ia1 = j;
					*--a1 = 1;
					*--c1 = os * *lx;
					*--l1 = 0;
					*--u1 = Infinity;
					}
				else
					*s = 'L';
				}
			}
			while(lx > LUv);
		for(; og; og = og->next)
			b[og->varno] = os * og->coef;
		memset((char *)kal, 0, m0*sizeof(int));
		ia1 = A_rownos;
		ia2 = ia1 + nzc0;
		while(ia1 < ia2)
			kal[*ia1++]++;
		for(i = j = 0; i < m0; i++)
			ka[i] = j += kal[i];
		a1 = A_vals + nzc0;
		j = n0;
		while(--j >= 0) {
			ia1 = A_rownos + A_colstarts[j];
			while(ia2 > ia1) {
				ia[i = --ka[*--ia2]] = j;
				a[i] = *--a1;
				}
			}
		if (nranges) {
			/* This could be merged with the following while
			 * loop but for cplex's undocumented assumption
			 * that elements of ka are in ascending order!
			 */
			kal2 = kal + m0;
			ka2 = ka + m0;
			for(i = m0; --i >= 0; ) {
				--ka2;
				--kal2;
				if (*--lrhs < *--urhs
				 && *lrhs > negInfinity
				 && *urhs < Infinity) {
					*--kal1 = j = *kal2;
					*--ka1 = nz1 -= j;
					memcpy((char *)(a+nz1),
						(char *)(a+*ka2),
						j*sizeof(double));
					memcpy((char *)(ia+nz1),
						(char *)(ia+*ka2),
						j*sizeof(int));
					*--c1 = os * *urhs;
					*--l1 = negInfinity;
					*--u1 = 0;
					}
				}
			lrhs = LUrhs + m0;
			urhs = Urhsx + m0;
			}
		i = m0;
		while(--i >= 0) {
			--c1;
			--l1;
			--u1;
			--urhs;
			if (*--lrhs <= negInfinity) {
				*c1 = os * *urhs;
				*l1 = negInfinity;
				*u1 = 0.;
				}
			else {
				*u1 = Infinity;
				*c1 = os * *lrhs;
				*l1 = *urhs == *lrhs ? negInfinity : 0.;
				}
			}
		free(A_vals);
		}
	else { /* primal */
		memset(c, 0, n0*sizeof(double));
#ifdef OBJ_ADJ
		if (objadj) {
			c[n0] = obj_adj; /* fixed adjustment to objective */
			*--l1 = 0;
			*--u1 = 2;
			*--ia1 = m0-nqc;
			*--a1 = *--b1 = 1;
			*--kal1 = 1;
			*--ka1 = --nz1;
			*--s = 'E';
			*--rv1 = 0;
			}
#endif
		for(; og; og = og->next)
			c[og->varno] = og->coef;
		i = m0;
		while(--i >= 0) {
			--b1;
			*--rv1 = 0;
			--s;
			--urhs;
			if ((*b1 = *--lrhs) <= negInfinity) {
				*s = 'L';
				*b1 = *urhs;
				}
			else if (*urhs >= Infinity)
				*s = 'G';
			else if (*lrhs < *urhs) {
				*s = 'R';
				*rv1 = *urhs - *lrhs;
				}
			else
				*s = 'E';
			}
		ka2 = ka + n0;
		*ka2 = nz1;
		while(ka2 > ka) {
			--ka2;
			*--kal1 = ka2[1] - ka2[0];
			}
		}
	mbas = m;
	nbas = n;
	cname = rname = 0;
	if (startbas || endbas || startsol || endsol || wrtfname
	 || startvec || endvec
	 || Optimize == CPXprimopt || Optimize == CPXdualopt
	 || Optimize == CPXlpopt || Optimize == CPXqpopt
	 || parval(CPX_PARAM_XXXIND)) {
		make_names(&cname, &rname);
		zap[3] = M1record((Char*)cname);
		}
	mq = m - nqc;
#ifdef Uselazy
	i = mq;
	zap[5] = lazyadj(asl, &LI, n, nint, nqc, &mq, b+nqc, senx+nqc, ka, kal,
				 ia, a, rngvec+nqc, rname ? rname+nqc : 0);
	if ((j = LI.rows[0] + LI.rows[1] - (i - mq)) > 0) {
		/* unlikely; allocate extra space... */
		j += n_con + mq + nqc + nlogc;
		d->y = (double*)M1zapalloc(j*(sizeof(double) + sizeof(int)));
		ia1 = (int*)(d->y + j);
		memcpy(ia1, d->rstat, n_con * sizeof(int));
		d->rsd->u.i = d->rstat = ia1;
		}
#endif /*Uselazy*/
	i = CPXcopylpwnames(Env, cpx, n, mq, objsen, c, b+nqc, senx+nqc, ka,
				kal, ia, a, l, u, rngvec+nqc, cname,
				rname ? rname+nqc : 0);
	if (i)
		badret("CPXcopylpwnames", i, 531);
#ifdef CPXERR_QCP_SENSE /* if CPLEX version >= 9.0 */
	lrhs = LUrhs;
	urhs = Urhsx;
	for(i = 0; i < nqc; i++) {
		j = CPXaddqconstr(Env, cpx, kalq[i], nelqc[i], b[i], senx[i],
			iaq+kaq[i], aq+kaq[i], rowqc[i], colqc[i], qcmat[i],
			rname ? rname[i] : 0);
		if (j)
			Surprise(asl, j, "CPXaddqconstr");
		}
	mbas += nqc;
#endif /*CPXERR_QCP_SENSE*/
#ifdef CPLEX_MIP
	if (nint1) {
		if ((i = CPXcopyctype(Env, cpx, ctype)))
			badret("CPXcopyctype", i, 531);
		if (nsos) {
			sos_kludge(nsos, sosbeg, sosref);
#ifdef CPX_PARAM_IISIND /* version < 9.2b */
			i = CPXcopysos(Env,cpx,nsos,nsosnz,sostype,sospri,
					sosbeg,sosind,sosref);
#else
			i = CPXcopysos(Env,cpx,nsos,nsosnz,sostype,
					sosbeg,sosind,sosref,0);
#endif
			free(sosref);
			if (i)
				badret("CPXcopysos", i, 531);
			}
		if (nint) {
			mip_priorities(asl, cpx);
			if (mipstval && X0 && (i = nbv1 + niv + nlvbi + nlvci + nlvoi))
				set_mipinit(asl, cpx, i);
			}
		if (workfiledir)
			CPXsetstrparam(Env, CPX_PARAM_WORKDIR, workfiledir);
		}
#endif
#ifdef Uselazy
	if (LI.rows[0]) {
		CPXgetintparam(Env, CPX_PARAM_REDUCE, &i);
		if (i & ~1) {
			i &= 1;
			fprintf(Stderr, "Assuming prereduce = %d for lazy constraints.\n", i);
			if ((j = CPXsetintparam(Env, CPX_PARAM_REDUCE, i)))
				badret("CPXsetintparam(CPX_PARAM_REDUCE)", j, 531);
			}
		if ((i = CPXaddlazyconstraints(Env, cpx, LI.rows[0],
			LI.nz[0], LI.rhs[0], LI.sense[0], LI.rowbeg[0], LI.colno[0],
			LI.matval[0], LI.rowname[0])))
				badret("CPXaddlazyconstraints", i, 531);
		}
	if (LI.rows[1]) {
		CPXgetintparam(Env, CPX_PARAM_PRELINEAR, &i);
		if (i) {
			fprintf(Stderr, "Assuming prelinear = 0 for user cuts.\n");
			if ((j = CPXsetintparam(Env, CPX_PARAM_PRELINEAR, 0)))
				badret("CPXsetintparam(CPX_PARAM_PRELINEAR)", j, 531);
			}
		if ((i = CPXaddusercuts(Env, cpx, LI.rows[1],
			LI.nz[1], LI.rhs[1], LI.sense[1], LI.rowbeg[1], LI.colno[1],
			LI.matval[1], LI.rowname[1])))
				badret("CPXaddusercuts", i, 531);
		}
#endif
#ifdef BARRIER
	havestart = 0;
	if (Optimize == Optimize2 || Optimize == Optimizebar)
		havestart = startcomp(asl, n, m, nextra, ka, kal, ia,
			a, b, c, &cdual, &cprim, &rdual, &rprim);
#endif
	if (nelqf) {
		nelq = (int)nelqf;
		colq = (int*)Malloc(2*sizeof(int)*nr);
		colqcnt = colq + nr;
		if (sizeof(fint) != sizeof(int)) {
			/* should not happen */
			rowq = (int*)Malloc(nelq*sizeof(int));
			for(i = 0; i < nelq; i++)
				rowq[i] = (int)rowqf[i];
			free(rowqf);
			}
		else
			rowq = (int*)rowqf;
		k = n_var;
		for(i = 0; i < k; i++) {
			colq[i] = (int)colqf[i];
			colqcnt[i] = (int)(colqf[i+1] - colqf[i]);
			}
		for(j = colqf[k]; i < nr; i++) {
			colq[i] = j;
			colqcnt[i] = 0;
			}
		free(colqf);
		for(i = 0; ; i++) {
			if (i == k) {
				i = qmatadj(k, nr, objsen, colq,
					colqcnt, &qmat);
				if (!i) {
					for(i = k; i < nr; i++)
						qmat[i++] = 0.;
					i = CPXcopyqpsep(Env,cpx,qmat);
					if (i)
						badret("CPXcopyqpsep",
							i,531);
					}
				free(rowq);
				free(qmat);
				free(colq);
				if (i)
					exit(1);
				baralgname = "separable QP ";
				break;
				}
			if (!(j = colqcnt[i]))
				continue;
			if (j > 1 || rowq[colq[i]] != i) {
				i = CPXcopyquad(Env,cpx,colq,
					colqcnt,rowq,qmat);
				free(rowq);
				free(qmat);
				free(colq);
				if (i)
					badret("CPXcopyquad", i, 531);
				baralgname = "QP ";
				break;
				}
			}
		}
#ifdef CPX_CON_INDICATOR
	if (nlogc && (i = indicator_constrs_ASL(asl, cpx, add_indic, errinfo))) {
		switch(i) {
		  case 1:
			badretfmt(563,"logical constraint %s is not an indicator constraint.\n",
				lcon_name(errinfo[0]));
			break;
		  case 2:
			badretfmt(563,"logical constraint %s is not an indicator constraint\n\
	due to bad comparison with %s.\n", lcon_name(errinfo[0]), var_name(errinfo[1]));
			break;
		  case 3:
			badret("CPXaddindconstr", i, 531);
			break;
		  default:
			badret("indicator_constrs_ASL", i, 531);
		  }
		exit(4);
		}
#endif /*CPX_CON_INDICATOR*/
#ifdef BARRIER
	if (havestart) {
		if ((i = CPXcopystart(Env,cpx,0,0,cprim,rprim,cdual,rdual)))
			badret("CPXcopystart", i, 0);
		}
#endif
	for(i = 6; i > 0; )
		if (zap[--i]) {
			free(*zap[i]);
			*zap[i] = 0;
			}
	get_statuses(asl, cpx, d);
	cbi.disp = parval(CPX_PARAM_SIMDISPLAY);
#ifdef CPLEX_MIP
	cbi.mipdisp = nint && parval(CPX_PARAM_MIPDISPLAY);
#endif
	if (prestats
#ifndef NO_BARRIER
	 || Optimize == Optimize2
#endif
				 ) {
		zap_lpcbf = 1;
		cbi.nint = nint;
		cbi.pres = parval(CPX_PARAM_PREIND);
		CPXsetlpcallbackfunc(Env, lpcbf, 0);
#ifndef NO_BARRIER
		if (Optimize != Optimize2)
#endif
			{
			zap_lpcbf = 2;
			CPXsetmipcallbackfunc(Env, lpcbf, 0);
			}
		}
	breaking = 1;
	}

/* The following routine was named "round", but HPUX 11 reportedly */
/* could not cope with this name. */
 static int
xround(double *x, fint n, int assign, double *w)
{
	double d, dx, *xe, y;
	int m = 0;

	dx = *w;
	for(xe = x + n; x < xe; x++) {
		y = floor(*x + 0.5);
		if ((d = *x - y) != 0.) {
			if (d < 0)
				d = -d;
			if (dx < d)
				dx = d;
			m++;
			if (assign)
				*x = y;
			}
		}
	*w = dx;
	return m;
	}

 static char *
failstat(int j)
{
	static char buf[32];
	switch(j) {	/* cases and text provided by Ed Klotz */
	  case CPXERR_NO_MEMORY:
		return "insufficient memory";
	  case CPXERR_NO_PROBLEM:
		return "no problem";
	  case CPXERR_SINGULAR:
		return "singular basis";
	  case CPXERR_PRIIND:
		return "invalid pricing setting";
#ifdef CPXERR_BOUNDS_INFEAS
	  case CPXERR_BOUNDS_INFEAS:
		return "contradictory bounds";
#endif
	  case CPXERR_PRESLV_INForUNBD:
		return "presolve: infeasible or unbounded problem";
	  case CPXERR_INDEX_RANGE:
		return "illegal list length";
#ifdef CPXERR_NET_IMBALANCE
	  case CPXERR_NET_IMBALANCE:
		return "an infeasible or unbounded network";
#endif
	  case CPXERR_BAD_METHOD:
		return "unknown method after network optimization";
#ifdef CPXERR_PRESLV_TIME_LIM /* CPLEX 11 */
	  case CPXERR_PRESLV_TIME_LIM:
		return "time limit exceeded in CPLEX's presolve";
	  case CPXERR_THREAD_FAILED:
		return "thread creation failed";
#endif
	  }
	sprintf(buf, "CPLEX error # %d.", j);
	return buf;
	}

 static int
statadjust(cpxlp *cpx, int stat)
{
#ifdef CPX_STAT_OPTIMAL	/* >= 8.0 {{*/
	static int statmap[23] = { 1, 3, 2, 31, 11, 9, -1, -1, -1, 5,
				 7, 4, 12, 57, 58, 59, 60, 61, 62, 24,
				 22, 23, 31 };
	int adj, dfeas, pfeas, solmeth, soltype;

	if (stat >= 101 && stat <= 119)
		return stat + (36 - 101);
	if (stat > 22 || stat <= 0) {
#ifndef CPXERR_ALGNOTLICENSED
#define CPXERR_ALGNOTLICENSED 32024
#endif
		if (stat == CPXERR_ALGNOTLICENSED)
			return (obj_no >= 0 && obj_no < nlo) || nlc ? 34 : 32;
#ifdef CPXMIP_POPULATESOL_LIM
		if (stat == CPXMIP_POPULATESOL_LIM
		 || stat == CPXMIP_OPTIMAL_POPULATED
		 || stat == CPXMIP_OPTIMAL_POPULATED_TOL)
			stat = 36;
		else
#endif
#ifdef CPX_STAT_FIRSTORDER
		if (stat == CPX_STAT_FIRSTORDER)
			stat = 35;
		else
#endif
		stat = 0;
		}
	else {
		adj = 0;
		pfeas = dfeas = 1;
		switch(stat) {
			case CPX_STAT_INFEASIBLE:
				CPXsolninfo(Env, cpx, &solmeth, &soltype, &pfeas, &dfeas);
				if (solmeth == CPX_ALG_BARRIER && (!pfeas || !dfeas)) {
					adj = pfeas ? 16 : 17;
					if (!dfeas)
						adj += 2;
					}
				/* no break */
			case CPX_STAT_UNBOUNDED:
				if (method < 0)
					stat = CPX_STAT_INFEASIBLE + CPX_STAT_UNBOUNDED - stat;
				break;
			case CPX_STAT_NUM_BEST:
				CPXsolninfo(Env, cpx, &solmeth, &soltype, &pfeas, &dfeas);
				if (solmeth == CPX_ALG_BARRIER) {
					adj = 16;
					if (!pfeas)
						adj++;
					if (!dfeas)
						adj += 2;
					}
				break;
			case CPX_STAT_ABORT_IT_LIM:
			case CPX_STAT_ABORT_TIME_LIM:
				CPXsolninfo(Env, cpx, &solmeth, &soltype, &pfeas, &dfeas);
				if (!dfeas)
					adj = 1;
				break;
			case CPX_STAT_ABORT_USER:
				CPXsolninfo(Env, cpx, &solmeth, &soltype, &pfeas, &dfeas);
				if (solmeth == CPX_ALG_BARRIER) {
					if (pfeas && dfeas)
						adj = 5;
					else {
						adj = dfeas ? 1 : 2;
						if (!pfeas)
							adj += 2;
						}
					}
			}
		stat = statmap[stat-1] + adj;
		}
#else /*}{*/
	if (stat > 100 && stat < 120)
		stat -= 66;
	else if (stat > 18 && stat < 44)
		stat -= 13;
	else if (stat == 1101)
		stat = 31;	/* should not happen */
	else if (stat < 0 || stat > 18)
		stat = 0;
	else if (stat >= 2 && stat <= 3 && Optimize == CPXdualopt)
		stat = 5 - stat;
#endif /*}}*/
	return stat;
	}

 static void
write_basis(cpxlp *cpx)
{
	if (CPXmbasewrite(Env, cpx, endbas)) {
		printf("Could not write ending basis to \"%s\"\n", endbas);
		need_nl = 0;
		}
	}

 static void
cud_cadjust(ASL *asl)
{
	int i, m;
	real *L, *U;

	m = n_con;
	L = LUrhs;
	U = Urhsx;
	for(i = 0; i < m; i++)
		if (L[i] <= negInfinity)
			L[i] = U[i];
	}

 static void
Kludge_getbase(dims *d, int *rstat)
{ /* Kludge around CPLEX bug in rstat values for <= constraints... */
	int i, m = mbas;
	char *rt = d->rtype;
	for(i = 0; i < m; i++)
		if (rt[i] == 'L' && rstat[i] == 0)
			rstat[i] = 2;
	}

 static int
send_statuses(ASL *asl, cpxlp *cpx, dims *d)
{
	int i, j, m, n;
	int *cs, *cstat, *rs, *rstat;
	real *L, *U;
	static int map[] = {3, 1, 4, 6, 0};

	if (!(Oinfo.wantsol & 1) && !amplflag)
		return 0;
	cstat = d->cstat;
	rstat = d->rstat;
	memset(cstat, 0, n_var*sizeof(int));
	memset(rstat, 0, n_con*sizeof(int));
	if (method > 0) {
		if ((i = CPXgetbase(Env, cpx, cstat, rstat))) {
 no_basis:
			if (i != CPXERR_NO_BASIS)
				badret("CPXgetbase", i, 0);
			d->csd->kind &= ~ASL_Sufkind_output;
			d->rsd->kind &= ~ASL_Sufkind_output;
			return 1;
			}
		Basedebug("spx.bas",cstat,rstat);
		Kludge_getbase(d, rstat);
		stat_map(cstat, n_var, map, 4, "outgoing cstat");
		stat_map(rstat, n_con, map, 4, "outgoing rstat");
		equ_adjust(cstat, rstat);
		if (costsens) {
			L = (real*)M1alloc(2*sizeof(real)*(nbas + mbas));
			U = L + nbas;
			if ((i = CPXobjsa(Env, cpx, 0, nbas-1, L, U)))
				badret("CPXobjsa", i, 0);
			else {
				suf_rput("down", ASL_Sufkind_var, L);
				suf_rput("up", ASL_Sufkind_var, U);
				suf_rput("current", ASL_Sufkind_var, d->c);
				}
			L = U + nbas;
			U = L + mbas;
			if ((i = CPXrhssa(Env, cpx, 0, mbas-1, L, U)))
				badret("CPXrhssa", i, 0);
			else {
				cud_cadjust(asl);
				suf_rput("down", ASL_Sufkind_con, L);
				suf_rput("up", ASL_Sufkind_con, U);
				suf_rput("current", ASL_Sufkind_con, LUrhs);
				}
			}
		return 0;
		}
	cs = (int*)M1alloc((mbas + nbas)*sizeof(int));
	rs = cs + nbas;
	if ((i = CPXgetbase(Env, cpx, cs, rs)))
		goto no_basis;
	m = n_con;
	n = n_var;
	Kludge_getbase(d, rs);
	L = LUrhs;
	U = Urhsx;
	for(i = 0; i < m; i++) {
		if (L[i] <= negInfinity)
			rstat[i] = cs[i] == CPX_AT_UPPER ? 1 : 4;
		else
			rstat[i] = L[i] == U[i]
				||   cs[i] == CPX_BASIC ? 3 : 1;
		}
	if (nranges) {
		for(j = 0; j < m; j++)
			if (L[j] < U[j]
			 && L[j] > negInfinity
			 && U[j] < Infinity
			 && cs[i++] == CPX_BASIC)
				rstat[j] = 4;
		}
	L = LUv;
	U = Uvx;
	for(j = 0; j < n; j++) {
		if (L[j] > negInfinity) {
			if (L[j])
				cstat[j] = cs[i++] == CPX_BASIC ? 3 : 1;
			else
				cstat[j] = rs[j] == CPX_BASIC ? 3 : 1;
			if (U[j] < Infinity)
				if ((U[j] ? cs[i++] : rs[j]) == CPX_BASIC)
					cstat[j] = 4;
			}
		else if (U[j] < Infinity) {
			if (U[j])
				cstat[j] = cs[i++] == CPX_BASIC ? 4 : 1;
			else
				cstat[j] = rs[j] == CPX_BASIC ? 4 : 1;
			}
		else
			cstat[j] = rs[j] == CPX_BASIC ? 2 : 1;
		}
	equ_adjust(cstat, rstat);
	return 0;
	}

 static int
send_dray(ASL *asl, cpxlp *cpx, int nelqf)
{
	real proof, *y, *y1, *ye;

	if (nelqf) {
		/* CPLEX 8.0 faulted in this case; */
		/* CPLEX 8.1 would issue a message. */
		return 2;
		}
	y = (real*)M1zapalloc(mbas*sizeof(real));
	if (CPXdualfarkas(Env, cpx, y, &proof))
		return 1;
	if (obj_no >= 0 && objtype[obj_no]) {
		/* circumvent CPLEX bug */
		for(y1 = y, ye = y + mbas; y1 < ye; ++y1)
			*y1 = -*y1;
		}
	suf_rput("dunbdd", ASL_Sufkind_con, y);
	return 0;
	}

 static int
send_ray(ASL *asl, cpxlp *cpx)
{
	int i, j, n, rc;
	real *r;

	if (CPXgetijdiv(Env, cpx, &i, &j))
		return 1;
	if ((n = n_var) < nbas)
		n = nbas;
	r = (real*)M1zapalloc(n*sizeof(real));
	if (!(rc = CPXgetray(Env, cpx, r)))
		suf_rput("unbdd", ASL_Sufkind_var, r);
	return rc;
	}

#ifdef CPX_PARAM_IISIND
 static void
iis_put(ASL *asl, int kind, int n, int *ind, int *stat)
{
	int i, *x;

	x = (int*)M1zapalloc((&asl->i.n_var_)[kind]*sizeof(int));
	suf_iput("iis", kind, x);
	for(i = 0; i < n; i++)
		x[ind[i]] = stat[i] + 1;
	}
#endif /*CPX_PARAM_IISIND*/

 static int
retfrom(char *hbuf, int k, const char *what)
{ return Sprintf(hbuf, "\nReturn %d from %s.", k, what); }

#ifdef CPX_PARAM_FEASOPTMODE /* >= 9.2b */
 /* use conflict finder to determine IIS */

 static int
confiis_put(ASL *asl, int kind, int n, int *ind, int *stat)
{
	int i, j, k, *x;
	static int stmap[7] = {0, 5, 6, 7, 4, 1, 3 };

	x = (int*)M1zapalloc((&asl->i.n_var_)[kind]*sizeof(int));
	suf_iput("iis", kind, x);
	for(i = k = 0; i < n; i++) {
		if ((j = stat[i] + 1) < 0 || j > 6)
			j = 8;
		else
			j = stmap[j];
		if ((x[ind[i]] = j))
			k++;
		}
	return k;
	}

 static int
send_confiis_ext(ASL *asl, cpxlp *cpx, char *hbuf)
{
	char *grptype;
	double *L, *U, *gp;
	int *ci, *grpbeg, *grpind, *gs, i, i1, j, j1, k, n, nc, nv, rv, *vi;
	static int	lbmap[4] = {0, 6, 1, 8},
			ubmap[4] = {0, 7, 3, 8},
			cmap[4]  = {0, 5, 4, 8};

	L = LUv;
	U = Uvx;
	k = n_var;
	n = n_con;
	for(i = 0; i < k; i++) {
		if (L[i] > negInfinity)
			n++;
		if (U[i] < Infinity)
			n++;
		}
	gp = (double*)Malloc(n*(sizeof(double)+3*sizeof(int)+1));
	grpbeg = (int*)(gp + n);
	grpind = grpbeg + n;
	gs = grpind + n;
	grptype = (char*)(gs + n);
	for(i = 0; i < n; i++) {
		grpbeg[i] = i;
		gp[i] = 1.;
		}
	for(i = j = 0; i < k; i++) {
		if (L[i] > negInfinity) {
			grpind[j] = i;
			grptype[j++] = CPX_CON_LOWER_BOUND;
			}
		if (U[i] < Infinity) {
			grpind[j] = i;
			grptype[j++] = CPX_CON_UPPER_BOUND;
			}
		}
	for(i = 0, k += n_con - nlc; j < k; j++) {
		grpind[j] = i++;
		grptype[j] = CPX_CON_LINEAR;
		}
	for(i = 0, k += nlc; j < k; j++) {
		grpind[j] = i++;
		grptype[j] = CPX_CON_QUADRATIC;
		}
	if ((i = CPXrefineconflictext(Env, cpx, n, n, gp, grpbeg, grpind, grptype))) {
		rv = retfrom(hbuf, i, "CPXrefineconflictext");
		goto ret;
		}
	if ((i = CPXgetconflictext(Env, cpx, gs, 0, n-1))) {
		rv = retfrom(hbuf, i, "CPXgetconflictext");
		goto ret;
		}
	for(i = 0; i < n; i++)
		if ((j = ++gs[i]) < 0 || j > 1)
			gs[i] = j == CPX_CONFLICT_MEMBER ? 2 : 3;
	vi = (int*)M1zapalloc((n_var+n_con)*sizeof(int));
	ci = vi + n_var;
	k = n_var;
	for(i = j = nv = 0; i < k; i++) {
		i1 = 0;
		if (L[i] > negInfinity) {
			vi[i] = lbmap[gs[j++]];
			i1 = 1;
			}
		if (U[i] < Infinity) {
			i1 = 1;
			j1 = ubmap[gs[j++]];
			if (vi[i] && L[i] == U[i])
				j1 = 2;
			vi[i] = j1;
			}
		nv += i1;
		}
	for(i = nc = 0, k = n_con; i < k; i++)
		if ((ci[i] = cmap[gs[j++]]))
			nc++;
	suf_iput("iis", ASL_Sufkind_var, vi);
	suf_iput("iis", ASL_Sufkind_con, ci);
	rv = Sprintf(hbuf, "\nReturning an IIS involving %d variables and %d constraints.",
		nv, nc);
 ret:
	free(gp);
	return rv;
	}

 static int
send_confiis(ASL *asl, cpxlp *cpx, char *hbuf)
{
	int cs, k, mc, mc1, nc, nc1, rv;
	int *colind, *colst, *rowind, *rowst;
	static char *abort_reason[7] = {
		"contradiction", "time limit", "iteration limit", "node limit",
		"objective limit", "memory limit", "user request" };

	if ((k = CPXrefineconflict(Env, cpx, &mc, &nc)))
		return retfrom(hbuf, k, "CPXrefineconflict");

	k = mc + nc;
	if (mc < 0 || nc < 0 || k == 0)
		return Sprintf(hbuf,
			"\nSurprise values mc = %d, nc = %d from CPXrefineconflict.",
			mc, nc);

	colind = (int*)Malloc(2*sizeof(int)*k);
	colst = colind + nc;
	rowind = colst + nc;
	rowst = rowind + mc;

	if ((k = CPXgetconflict(Env, cpx, &cs, rowind, rowst, &mc1, colind, colst, &nc1))) {
		rv = retfrom(hbuf, k, "CPXgetconflict");
		goto ret;
		}
	if (cs == CPX_STAT_CONFLICT_FEASIBLE) {
		rv = Sprintf(hbuf, "\nNo IIS after all: problem is feasible!");
		goto ret;
		}
	if (cs != CPX_STAT_CONFLICT_MINIMAL) {
		if ((cs < CPX_STAT_CONFLICT_ABORT_CONTRADICTION) < 0
		 || cs > CPX_STAT_CONFLICT_ABORT_CONTRADICTION + 6)
			rv = Sprintf(hbuf,
				"\nSurprise conflict status = %d from CPXgetconflict\n", cs);
		else
			rv = Sprintf(hbuf, "\nSearch for conflicts aborted because of %s\n",
				abort_reason[cs - CPX_STAT_CONFLICT_ABORT_CONTRADICTION]);
		goto ret;
		}
	if (mc1 > mc || nc1 > nc) {
		rv = Sprintf(hbuf, "\nSurprise confnumrows = %d (should be <= %d),"
				"\nconfnumcols = %d (should be <= %d) from CPXgetconflict.",
				mc1, mc, nc1, nc);
		goto ret;
		}
	nc = confiis_put(asl, ASL_Sufkind_var, nc1, colind, colst);
	mc = confiis_put(asl, ASL_Sufkind_con, mc1, rowind, rowst);
	return Sprintf(hbuf,
		"\nReturning an IIS of %d variables and %d constraints.", nc, mc);

 ret:
	free(colind);
	return rv;
	}
#endif /*CPX_PARAM_FEASOPTMODE*/

 static int
send_iis(ASL *asl, cpxlp *cpx, char *hbuf)
{
	char *ckind;
	int i, j, m, m1, n, n1;
	int *ci, *cs, *ri, *rs;

	if (!want_iis)
		return 0;
#ifdef CPX_PARAM_FEASOPTMODE /* >= 9.2b */
	if (nlc)
		return send_confiis_ext(asl, cpx, hbuf);
	if (want_iis > 2
	|| (nbv + niv + nlvoi + nlvci + nlvbi > 0 && !relax)
	|| (Optimize != CPXprimopt && Optimize != CPXdualopt))
		return send_confiis(asl, cpx, hbuf);
#endif
	if (nlc) {
		printf("Ignoring iisfind request because of nonlinearities.\n");
		return 0;
		}
#ifdef CPX_PARAM_IISIND
	CPXsetintparam(Env, CPX_PARAM_IISIND, want_iis-1);
	if (i = CPXfindiis(Env, cpx, &m, &n))
		return retfrom(hbuf, i, "CPXfindiis");
#else
	if ((i = CPXrefineconflict(Env, cpx, &m, &n)))
		return retfrom(hbuf, i, "CPXrefineconflict");
#endif
	ci = (int*)M1alloc((2*(m+n))*sizeof(int));
	cs = ci + n;
	ri = cs + n;
	rs = ri + m;
	ckind = "";
#ifdef CPX_PARAM_IISIND
	if (i = CPXgetiis(Env, cpx, &j, ri, rs, &m1, ci, cs, &n1))
		return retfrom(hbuf, i, "CPXgetiis");
	iis_put(asl, ASL_Sufkind_var, n, ci, cs);
	iis_put(asl, ASL_Sufkind_con, m, ri, rs);
	if (j != 1)
		ckind = "partial ";
#else
	if ((i = CPXgetconflict(Env, cpx, &j, ri, rs, &m1, ci, cs, &n1)))
		return retfrom(hbuf, i, "CPXgetconflict");
	switch(j) {
	  case CPX_STAT_CONFLICT_MINIMAL:
		confiis_put(asl, ASL_Sufkind_var, n, ci, cs);
		confiis_put(asl, ASL_Sufkind_con, m, ri, rs);
		break;
	  case CPX_STAT_CONFLICT_FEASIBLE:
		return Sprintf(hbuf, "\nCPXgetconflict claims the problem is feasible.");
	  default:
		return Sprintf(hbuf, "\nSurprise conflict status %d from CPXgetconflict.", j);
	  }
#endif
	return Sprintf(hbuf,
		"\nReturning %siis of %d variables and %d constraints.",
		ckind, n, m);
	}

 static int
wantray(int dual, cpxlp **pcpx, int *itcp, int *itcip, int nelqf)
{
	Optalg Rayalg;
	cpxlp *cpx, *cpx1;
	int i;
	static Optalg PDalg[2] = { CPXprimopt, CPXdualopt };

	*itcp = 0;
	if (want_iis || !(rays & (1 << dual)))
		return 0;
	Rayalg = PDalg[dual];
	if ((!dual && Optimize == CPXprimopt)
	 || (nelqf && (Optimize == CPXdualopt || Optimize == CPXprimopt))) {
		*itcip = 0;
		return 1;
		}
	if (Optimize == CPXmipopt || Optimize == Optimizebar)
		return 0;
	cpx = *pcpx;
	if (Resolve && dual) {
		i = 0;
		CPXgetintparam(Env, CPX_PARAM_PREIND, &i);
		if (i) {
			if (!(cpx1 = CPXcloneprob(Env, cpx, &i)))
				badret("CPXcloneprob", i, 531);
			CPXfreeprob(Env, pcpx);
			*pcpx = cpx = cpx1;
			CPXsetintparam(Env, CPX_PARAM_PREIND, 0);
			}
		}
	(*Rayalg)(Env, cpx);
	*itcp = CPXgetitcnt(Env, cpx);
	*itcip = CPXgetphase1cnt(Env, cpx);
	algname = dual ? "dual " : "";
	return 1;
	}

#ifdef CPX_PARAM_POPULATELIM /*{*/
 static void
poolwrite(ASL *asl, cpxlp *cpx, dims *d, int nelqf, char *hb, int *hbi, size_t hblen)
{
	char buf[32], buf1[1200], *fname, *s;
	const char *whence;
	int i, j, k, k1, m, n, nm1, npm, nprt, nsols, oprt, rft, stat1, status;
	int *cstat, *ocstat, *orstat, *rstat;
	real feastol, inttol, obj, *x, *y, *y1;
	size_t L;
	static char
		errfmt1[] = "\nCPLEX solution status %d with fixed integers "
			"on solution pool\nmember %d (file %s):\n\t%s",
		errfmt2[] = "\nCPLEX status %d from CPXsolution "
			"on solution pool\nmember %d (file %s).";

	i = *hbi;
	j = rft = 0;
	y = 0;
	ocstat = orstat = 0;
	if (populate == 1) {
		if ((k = CPXpopulate(Env, cpx))) {
			i += snprintf(hb+i, hblen-i, "\nReturn %d from CPXpopulate.", k);
			goto done;
			}
		}
	nsols = CPXgetsolnpoolnumsolns(Env,cpx);
	if (nsols <= 0) {
		i += snprintf(hb+i, hblen-i, "\n%d solutions in solution pool.", nsols);
		goto done;
		}
	m = n_con;
	n = npm = n_var;
	nprt = 0; /* silence buggy warning */
	if ((k = pooldual && m > 0 && !nlc)) {
		nprt = nelqf ? CPXPROB_FIXEDMIQP : CPXPROB_FIXEDMILP;
		npm += m;
		}
	L = strlen(poolstub);
	x = (real*)Malloc(npm*(sizeof(real) + k*sizeof(int)) + L + 32);
	fname = (char*)(x + n);
	if (k) {
		y = x + n;
		if (d->cstat && (rstat = d->rstat)) {
			cstat = d->cstat;
			ocstat = (int*)(y + m);
			orstat = ocstat + n;
			fname = (char*)(orstat + m);
			memcpy(ocstat, cstat, n*sizeof(int));
			memcpy(orstat, rstat, m*sizeof(int));
			}
		if ((rft = !CPXgetdblparam(Env,CPX_PARAM_EPRHS,&feastol)
		 && !CPXgetdblparam(Env,CPX_PARAM_EPINT,&inttol)
		 && inttol > feastol))
			CPXsetdblparam(Env, CPX_PARAM_EPRHS, inttol);
		}
	strcpy(fname, poolstub);
	nm1 = n - 1;
	whence = 0;
	for(j = k = 0; j < nsols; j++) {
		if ((k = CPXgetsolnpoolx(Env, cpx, j, x, 0, nm1))) {
			whence = "CPXgetsolnpoolx";
			break;
			}
		if ((k = CPXgetsolnpoolobjval(Env, cpx, j, &obj))) {
			whence = "CPXgetsolnpoolobjval";
			break;
			}
		g_fmtop(buf, obj);
		k1 = snprintf(buf1, sizeof(buf1),
			"Solution pool member %d (of %d); objective %s",
			j+1, nsols, buf);
		sprintf(fname+L, "%d.sol", j+1);
		if ((y1 = y)) {
			oprt = CPXgetprobtype(Env, cpx);
			if ((k = CPXchgprobtypesolnpool(Env, cpx, j, nprt))) {
				whence = "CPXchgprobtypesolnpool";
				break;
				}
			status = CPXprimopt(Env, cpx);
			stat1 = CPXgetstat(Env, cpx);
			if (status || stat1 != 1) {
				snprintf(buf1+k1, sizeof(buf1)-k1, errfmt1,
					stat1, j, fname,
					solinfo[statadjust(cpx,stat1)].msg);
				y1 = 0;
				}
			else if (CPXsolution(Env, cpx, &stat1, &obj, x, y, 0, 0)) {
				y1 = 0;
				status = CPXgetstat(Env, cpx);
				stat1 = statadjust(cpx, status);
				snprintf(buf1+k1, sizeof(buf1)-k1, errfmt2,
					status, j, fname);
				}
			else if (send_statuses(asl, cpx, d))
				snprintf(buf1+k1, sizeof(buf1)-k1, "\nNo basis.");
			if ((k = CPXchgprobtype(Env, cpx, oprt))) {
				whence = "CPXchgprobtype<restore>";
				break;
				}
			}
		if (write_solf_ASL(asl, buf1, x, y1, 0, fname)) {
			whence = "";
			break;
			}
		}
	i += snprintf(hb+i,hblen-i,"\nWrote %d solution%s in solution pool",
			j, "s" + (j == 1));
	s = L > 32 ? "\n" : " ";
	if (j > 0)
	  switch(j) {
		case 1:
			i += snprintf(hb+i,hblen-i,"\nto file %s1.sol", poolstub);
			break;
		case 2:
			i += snprintf(hb+i,hblen-i,"\nto files %s1.sol%sand %s2.sol",
				poolstub, s, poolstub);
			break;
		default:
			i += snprintf(hb+i,hblen-i,"\nto files %s1.sol%s... %s%d.sol",
				poolstub, s, poolstub, j);
		}
	i += snprintf(hb+i,hblen-i,".");
	if (whence) {
		if (k)
			i += snprintf(hb+i,hblen-i,
				"\nSolution pool writing stopped by return %d from %s.",
				k, whence);
		else
			i += snprintf(hb+i,hblen-i, "\nCould not open \"%s\"", fname);
		}
	if (ocstat) {
		memcpy(cstat, ocstat, n*sizeof(int));
		memcpy(rstat, orstat, m*sizeof(int));
		}
	free(x);
 done:
	if (rft) /*restore "feasibility" tolerance*/
		CPXsetdblparam(Env, CPX_PARAM_EPRHS, feastol);
	d->npool = j;
	suf_iput("npool", ASL_Sufkind_prob, &d->npool);
	*hbi = i;
	}

 struct
Element
{
	keyword *kw;
	union {
		double dval;
		int ival;
		} u;
	};

 typedef struct
ParInfo
{
	double align;
	Element *nd;
	AVL_Tree *tree;
	int n_dbl, n_int;
	} ParInfo;

 static int
int_elcmp(void *v, const Element *L, const Element *R)
{ return (Intcast L->kw->info) - (Intcast R->kw->info); }

 static ParInfo *
getparinfo(void)
{
	AVL_Tree *at;
	Element *e;
	ParInfo *PI;
	double w[4];
	int k, z[4];
	keyword *kw, *kwe;
	size_t L, n;

	/* Use of AVL stuff here is only needed if the choice between */
	/* (sf_int and sf_int2) or (sf_dbl and sf_dbl2) is wrong in keywds. */

	n = nkeywds;
	L = n*sizeof(Element) + sizeof(ParInfo); /* overkill */
	PI = (ParInfo*)Malloc(L);
	PI->nd = e = (Element*)(PI + 1);
	PI->tree = at = AVL_Tree_alloc(0, int_elcmp, mymalloc_ASL);
	PI->n_dbl = PI->n_int = 0;
	kw = keywds;
	for(kwe = kw + nkeywds; kw < kwe; ++kw) {
		k = Intcast kw->info;
		if (kw->kf == sf_int) {
			CPXinfointparam(Env, k, z, z+1, z+2);
			CPXgetintparam(Env, k, z+3);
			if (z[0] != z[3]) {
				e->kw = kw;
				if (!AVL_insert(e,at)) {
					e->u.ival = z[3];
					++e;
					PI->n_int++;
					}
				}
			}
		else if (kw->kf == sf_dbl) {
			CPXinfodblparam(Env, k, w, w+1, w+2);
			CPXgetdblparam(Env, k, w+3);
			if (w[0] != w[3]) {
				e->kw = kw;
				if (!AVL_insert(e,at)) {
					e->u.dval = w[3];
					++e;
					PI->n_dbl++;
					}
				}
			}
		}
	return PI;
	}

 static void
tunewrite(char *fname, const char *fkeywd, char *hb, int *hbi, size_t hblen)
{
	Element *e, *ee;
	FILE *f;
	ParInfo *PI;
	int i, n;
	keyword *kw;

	i = *hbi;
	if (!(f = fopen(fname, "w"))) {
		*hbi = i += snprintf(hb+i, hblen-i, "\nCould not open %s file \"%s\"",
				fkeywd, fname);;
		solve_result_num = 565;
		return;
		}
	PI = getparinfo();
	n = 0;
	e = PI->nd;
	n = AVL_Tree_size(PI->tree);
	ee = e + n;
	for(; e < ee; ++e) {
		kw = e->kw;
		if (kw->kf == sf_int)
			fprintf(f, "%s = %d\n", kw->name, e->u.ival);
		else
			fprintf(f, "%s = %.g\n", kw->name, e->u.dval);
		}
	fclose(f);
	*hbi = i += snprintf(hb+i, hblen-i, "\n%d settings written to %s file \"%s\"",
			n, fkeywd, fname);
	AVL_Tree_free(&PI->tree);
	free(PI);
	}
#endif /*}*/

#if CPX_VERSION >= 1000 /*{*/
 static void
tunewriteprm(char *fname, const char *fkeywd, char *hb, int *hbi, size_t hblen)
{
	char *fmt;
	int i = *hbi;

	if (CPXwriteparam(Env, fname)) {
		fmt = "\nError writing %s file \"%s\"";
		solve_result_num = 565;
		}
	else
		fmt = "\n%s file \"%s\" written";
	*hbi = i += snprintf(hb+i, hblen-i, fmt, fkeywd, fname);
	}
#endif /*}*/

#if CPX_VERSION >= 1100 /*{*/

 typedef struct
Tunefix {
	AVL_Tree *avi, *avr;
	Element *e;
	int *num0, *nume, *numi, *numr;
	} Tunefix;

 static void
tunefixname(Tunefix *Tf, char *s, const char *what)
{
	AVL_Tree *T;
	Element *e;
	Kwfunc *f;
	char *peq;
	int L;
	keyword *kw;

	if ((kw = (keyword *)b_search_ASL(keywds, (int)sizeof(keyword),
			nkeywds, &s, &peq))) {
		f = kw->kf;
		e = Tf->e;
		e->kw = kw;
		if (f == sf_int || f == sf_int2) {
			if (!(T = Tf->avi))
				Tf->avi = T = AVL_Tree_alloc(0, int_elcmp, mymalloc_ASL);
			if (!AVL_insert(e,T)) {
				*Tf->numi++ = Intcast kw->info;
				++Tf->e;
				}
			}
		else if (f == sf_dbl || f == sf_dbl2) {
			if (!(T = Tf->avr))
				Tf->avr = T = AVL_Tree_alloc(0, int_elcmp, mymalloc_ASL);
			if (!AVL_insert(e,T)) {
				*--Tf->numr = Intcast kw->info;
				++Tf->e;
				}
			}
		}
	else {
		for(L = 0; s[L] > ' '; ++L);
		badretfmt(566, "unrecognized keyword \"%.*s\" in %s",
			L, s, what);
		}
	}

 static void
tunefix_proc(Tunefix *Tf, char *s)
{
	char buf[64], *s1;
	size_t inc, L;

	for(;;) {
		while(*s <= ' ')
			if (!*s++)
				return;
		for(L = inc = 0; s[L] > 0; ++L) {
			if (s[L] == ',') {
				inc = 1;
				break;
				}
			}
		if (L && L < sizeof(buf)) {
			s1 = s;
			if (inc) {
				memcpy(s1 = buf, s, L);
				s1[L] = 0;
				}
			tunefixname(Tf, s1, "tunefix");
			}
		s += L + inc;
		}
	}

 static void
tunefixfile_proc(Tunefix *Tf)
{
	FILE *f;
	char buf[4096], *b, *b1, *be;

	if (!(f = fopen(tunefixfile, "rb")))
		badretfmt(566, "Could not open tunefixfile \"%s\"", tunefixfile);
 top:
	while(fgets(buf, sizeof(buf), f)) {
		b = buf;
		be = b + strlen(b);
		if (be > b && be[-1] == '\n')
			*--be = 0;
		for(;;) {
			while(*b <= ' ')
				if (!*b++)
					goto top;
			b1 = b;
			while(*++b > ' ');
			tunefixname(Tf, b1, "tunefixfile");
			}
		}
	fclose(f);
	}

 static void
tunerun(cpxlp *cpx, char *hb, int *hbi, size_t hblen)
{
	Element *E;
	Tunefix TF;
	char *what, whatbuf[64];
	double *dblval;
	int *intval;
	int i, j, k, n_dbl, n_int;

	dblval = 0;
	intval = 0;
	TF.avi = TF.avr = 0;
	TF.e = E = (Element *)Malloc((nkeywds)*(sizeof(Element) + sizeof(int)));
	TF.num0 = TF.numi = (int*)(E + nkeywds);
	TF.nume = TF.numr = TF.num0 + nkeywds;
	if (tunefix)
		tunefix_proc(&TF, tunefix);
	if (tunefixfile)
		tunefixfile_proc(&TF);
	n_dbl = TF.nume - TF.numr;
	n_int = TF.numi - TF.num0;
	if (n_dbl + n_int) {
		dblval = (double*)Malloc(n_dbl*sizeof(double) + n_int*sizeof(int));
		intval = (int*)(dblval + n_dbl);
		for(i = 0; i < n_int; ++i)
			CPXgetintparam(Env, TF.num0[i], intval + i);
		for(i = 0; i < n_dbl; ++i)
			CPXgetdblparam(Env, TF.numr[i], dblval + i);
		}
	k = CPXtuneparam(Env, cpx, n_int, TF.num0, intval, n_dbl, TF.numr, dblval,
			0, 0, 0, &j);
	i = *hbi;
	if (k)
		i += snprintf(hb + i, hblen - i,
				"\nSurprise return %d from CPXtuneparam.", k);
	else {
		switch(j) {
		 case 0:
			what = "finished";
			break;
		 case CPX_TUNE_ABORT:
			what = "aborted";
			break;
		 case CPX_TUNE_TILIM:
			what = "time limit reached";
			break;
		 default:
			what = whatbuf;
			snprintf(whatbuf, sizeof(whatbuf),
				"got surprise tunestat = %d", j);
		 }
		i += snprintf(hb + i, hblen - i, "\nTuning %s.", what);
		}
	*hbi = i;
	if (dblval) {
		free(dblval);
		if (TF.avi)
			AVL_Tree_free(&TF.avi);
		if (TF.avr)
			AVL_Tree_free(&TF.avr);
		}
	free(E);
	}
#endif /*}*/

 static void
amplout(ASL *asl, cpxlp **pcpx, dims *d, int status, int nelqf, int nint1, int *nosp)
{
	char buf[32], hbuf[4096], *wb;
	cpxlp *cpx, *cpx1;
	int bitc, i, ii, itc, itci, j, m, mipstat, n, needsol, nint;
	int nodecnt, nodelim, nos, npt, nround, opt, stat, stat0, stat1, stat10;
	real *bb, *bn, *l, *le, *u, w[4], *x, *x2, *y, *z, *z1;
	real absmipgap, bbound, bcond, bobj, bobj0, feastol, inttol;
	real relmipgap, obj, t;
	Sol_info *SI;
	static Sol_info solinfo1[] = {
	 { "QP Hessian is not positive semi-definite", 542, 0 },
	 { "primal optimal (no dual solution available)", 004, 1 },
	 { "ran out of memory", 522, 0 }
#ifdef CPXERR_Q_NEG_ZERO_COMP	/* CPLEX 7.1 */
	 ,{ "QP Hessian has diag. elements of the wrong sign", 541, 0 }
#endif
		};
	static char *netmsg[] = {
		"no network",
		/*"an optimal network solution"*/ "",
		"an infeasible network",
		"an unbounded network",
		"an infeasible or unbounded network",
		"iteration limit with a feasible network solution",
		"iteration limit with an infeasible network solution",
		"time limit with a feasible network solution",
		"time limit with an infeasible network solution",
		"a feasible aborted network solution",
		"an infeasible aborted network solution",
		"a bug"
#ifdef Student_Edition
		,"the adjusted problem too large"
#endif
		 };
	static char *dray_msg[4] = {
		"\nconstraint.dunbdd returned",
		"\nfailed to compute variable.dunbdd",
		"\nQP simplex cannot compute variable.dunbdd"
		};
#if CPX_VERSION < 1100 /*{{*/
	int cl0, cl, cov, ccpr = 0;
#else /*}{*/
	typedef struct
CutInfo {
		const char *cutname;
		int cuttype;
		} CutInfo;
	CutInfo *CI;
	static CutInfo Cut_Info[] = {
		{ "cover",			CPX_CUT_COVER },
		{ "GUB-cover",			CPX_CUT_GUBCOVER },
		{ "flow-cover",			CPX_CUT_FLOWCOVER },
		{ "clique",			CPX_CUT_CLIQUE },
		{ "Gomory",			CPX_CUT_FRAC },
		{ "flow-path",			CPX_CUT_FLOWPATH },
		{ "disjunctive",		CPX_CUT_DISJ },
		{ "implied-bound",		CPX_CUT_IMPLBD },
		{ "zero-half",			CPX_CUT_ZEROHALF },
		{ "mixed-integer rounding",	CPX_CUT_MIR },
#ifdef CPX_CUT_LOCALCOVER
		{ "local cover",		CPX_CUT_LOCALCOVER },
#endif
#ifdef CPX_CUT_MCF
		{ "multi-commodity flow",	CPX_CUT_MCF },
#endif
#ifdef CPX_CUT_LANDP
		{ "split",			CPX_CUT_LANDP },
#endif
#if 0
		{ "",	CPX_CUT_TIGHTEN },
		{ "",	CPX_CUT_OBJDISJ },
		{ "",	CPX_CUT_USER },
		{ "",	CPX_CUT_TABLE },
		{ "",	CPX_CUT_SOLNPOOL }
#endif
		{ 0, 0}};
#ifdef CPX_PARAM_MIPKAPPASTATS /*{*/
	typedef struct
QualityInfo {
	const char *desc;
	int key;
	} QualityInfo;
	enum {QinfoN = 6};
	static QualityInfo QualInfo[QinfoN] = {
		{"max. condition number of optimal bases", CPX_KAPPA_MAX},
		{"attention statistic (weighted sum of percentages of bases with\n"
		 "high condition numbers", CPX_KAPPA_ATTENTION},
		{"percentage of stable bases (condition number < 1e7)",
			CPX_KAPPA_STABLE},
		{"percentage of suspect bases (condition between 1e7 and 1e10)",
			CPX_KAPPA_SUSPICIOUS},
		{"percentage of unstable bases (condition number between 1e10 and 1e14)",
			CPX_KAPPA_UNSTABLE},
		{"percentage of ill-posed bases (condition number > 1e14)",
			CPX_KAPPA_ILLPOSED}};
	int haveqmet = 0;
	real qmet[QinfoN];
#endif /*}*/

#endif /*}}*/

	cpx = *pcpx;
	absmipgap = obj = relmipgap = Infinity;
	stat0 = stat10 = -1;
	stat1 = 1;
	ii = mipstat = nodecnt = npt = nround = opt = 0;
	m = n_con;
	n = n_var;
	x = d->x;
	y = d->y;
	if (Logf)
		CPXsetlogfile(Env, NULL); /* suppress idiotic Error 1217 */
#ifdef BARRIER
	bitc = CPXgetbaritcnt(Env, cpx);
#else
	bitc = 0;
#endif
#ifdef CPLEX_MIP
	bobj = bobj0 = objsen * Infinity;
#endif
	itc = CPXgetitcnt(Env, cpx);
	itci = CPXgetphase1cnt(Env, cpx);
	if (Logf)
		CPXsetlogfile(Env, Logf);
	nos = nosp ? *nosp : 0;
	switch (status) {
	  case 0: break;
	  case CPXERR_NO_MEMORY:
		SI = solinfo1 + 2;
		stat = -1;
		goto have_SI;
	  case CPXERR_Q_NOT_POS_DEF:
		SI = solinfo1;
		stat = -1;
		goto have_SI;
#ifdef CPXERR_Q_NEG_ZERO_COMP
	  case CPXERR_Q_NEG_ZERO_COMP:
		SI = solinfo1 + 3;
		stat = -1;
		goto have_SI;
#endif
	  case CPXERR_PRESLV_INForUNBD:
		SI = solinfo + 31;
		goto have_SI;
	  case 32024:
		stat = 33;
		goto more_default;
	  default:
		stat = 0;
 more_default:
		stat0 = status;
		x = y = 0;
		goto have_stat;
		}
	if (Optimize == CPXlpopt
#ifndef   BARRIER_FOR_QP
				 || Optimize == CPXqpopt
#endif /* BARRIER_FOR_QP */
				) {
		if (lpoptalg == CPX_ALG_AUTOMATIC)
			switch(i = CPXgetmethod(Env,cpx)) {
			  case CPX_ALG_DUAL:
				Optimize = CPXdualopt;
				algname = "dual ";
				break;
			  case CPX_ALG_PRIMAL:
				Optimize = CPXprimopt;
				break;
#ifdef CPX_STAT_OPTIMAL	/* >= 8.0 */
			  case CPX_ALG_SIFTING:
			  case CPX_ALG_BARRIER:
#ifdef CPX_ALG_FEASOPT /* >= 9.2b */
			  case CPX_ALG_FEASOPT:
#endif

				Optimize = CPXlpopt;
				lpoptalg = i;
#endif /* CPX_STAT_OPTIMAL */
			  }
		}
#ifdef CPLEX_MIP
	else if (Optimize == CPXmipopt) {
#if CPX_VERSION < 1100 /*{*/
		if (prestats && nint1
		 && (parval(CPX_PARAM_CLIQUES) || parval(CPX_PARAM_COVERS))) {
			cl0 = CPXgetgenclqcnt(Env, cpx);
			cl  =  CPXgetclqcnt(Env, cpx);
			cov =  CPXgetcovcnt(Env, cpx);
			ccpr = 1;
			}
#endif /*}*/
		if (status > 0) {
			mipstat = statadjust(cpx, stat0 = status);
			status = 0;
			}
		}
#endif
#ifdef CPLEX_MIP
	if (nint1) {
#ifdef CPX_PARAM_MIPKAPPASTATS /*{*/
		if (CPXgetintparam(Env, CPX_PARAM_MIPKAPPASTATS, &j) == 0 && j >= 0) {
		    memset(qmet, 0, sizeof(qmet));
		    j = CPXgetdblquality(Env, cpx, &qmet[0], QualInfo[0].key);
		    if (j == 0 && qmet[0] > 0.) {
			haveqmet = 1;
			for(j = 1; j < QinfoN; ++j) {
				CPXgetdblquality(Env, cpx, &qmet[j], QualInfo[j].key);
				if (qmet[j] > 0.)
					haveqmet = j;
				}
			}
		    }
#endif /*}*/
		if (mipstat) {
			stat = mipstat;
			if (mipstat < 4 || mipstat > 8)
				goto have_ii;
			switch(mipstat) {
			  /* Increase limits to permit recovering best sol'n */
			  case 5: case 6:
				CPXsetintparam(Env, CPX_PARAM_ITLIM,
						2*parval(CPX_PARAM_ITLIM));
				break;
			  case 7: case 8:
				{ double tl = 1e75;
				CPXgetdblparam(Env, CPX_PARAM_TILIM, &tl);
				CPXsetdblparam(Env, CPX_PARAM_TILIM, 2.*tl);
				}
			  }
			}
		else {
			stat0 = CPXgetstat(Env, cpx);
			if (stat0 == CPXMIP_INForUNBD && Resolve) {
				CPXgetintparam(Env, CPX_PARAM_PREIND, &i);
				if (i) {
					if (!(cpx1 = CPXcloneprob(Env, cpx, &i)))
						badret("CPXcloneprob", i, 531);
					CPXsetintparam(Env, CPX_PARAM_PREIND, 0);
					CPXfreeprob(Env, &cpx);
					CPXmipopt(Env, *pcpx = cpx = cpx1);
					stat0 = CPXgetstat(Env, cpx);
					}
				}
			stat = statadjust(cpx, stat0);
			}
		ii = CPXgetmipitcnt(Env, cpx);
		nodecnt = CPXgetnodecnt(Env, cpx);
 have_ii:
		if (CPXgetbestobjval(Env, cpx, &bobj))
			bobj = bobj0;
		if (!solinfo[stat].wantobj)
			x = y = 0;
		if (solinfo[stat].wantobj || endbas) {
			if (zap_lpcbf) {
				if (zap_lpcbf == 2)
					CPXsetmipcallbackfunc(Env, breakcallback, 0);
				CPXsetlpcallbackfunc(Env, breakcallback, 0);
				zap_lpcbf = 0;
				}
			if (!CPXgetmipobjval(Env, cpx, &obj) && bobj != bobj0) {
				if ((absmipgap = obj - bobj) < 0.)
					absmipgap = -absmipgap;
				if ((t = obj) < 0.)
					t = -t;
				relmipgap = absmipgap / (1e-10 + t);
				}
			if (endbas)
				mipbasis = 1;
			needsol = 0;
			if ((i = CPXgetmipx(Env, cpx, x, 0, n-1))) {
				fprintf(Stderr,
					"Surprise return %d from CPXgetmipx\n", i);
				x = 0;
				}
			if (!nlc && mipbasis) {
				opt = CPXgetprobtype(Env, cpx);
				npt = nelqf ? 8 : 3;
				if (npt != opt)
					CPXchgprobtype(Env, cpx, npt);
				if ((i = !CPXgetdblparam(Env,CPX_PARAM_EPRHS,&feastol)
				 && !CPXgetdblparam(Env,CPX_PARAM_EPINT,&inttol)
				 && inttol > feastol))
					CPXsetdblparam(Env, CPX_PARAM_EPRHS, inttol);
				status = CPXprimopt(Env, cpx);
				stat1 = CPXgetstat(Env, cpx);
				if (i)	/*restore "feasibility" tolerance*/
					CPXsetdblparam(Env, CPX_PARAM_EPRHS, feastol);
				if (status || stat1 != 1) {
					fprintf(Stderr,
					"CPLEX solution status %d with fixed integers:\n\t%s\n",
						stat1, solinfo[statadjust(cpx,stat1)].msg);
					y = 0;
					}
				else if (x)
					needsol = 1;
				if (endbas && !nelqf)
					write_basis(cpx);
				}
			if (!solinfo[stat].wantobj)
				goto have_stat;
			if (needsol) {
				if (CPXsolution(Env, cpx, &stat1, &obj, x, y, 0, 0))
					stat1 = CPXgetstat(Env, cpx);
				stat1 = statadjust(cpx, stat10 = stat1);
				}
			/* round to integer */
			if (!relax && Round >= 0 && x) {
				w[0] = 0;
				if ((nint = niv + nbv)) {
					x2 = x + n_var - nint;
					nround = xround(x2, nint, Round & 1, w);
					}
				if ((nint = nlvbi)) {
					x2 = x + (nlvb - nint);
					nround += xround(x2, nint, Round & 1, w);
					}
				if ((nint = nlvci)) {
					x2 = x + (nlvc - nint);
					nround += xround(x2, nint, Round & 1, w);
					}
				if ((nint = nlvoi)) {
					x2 = x + (nlvo - nint);
					nround += xround(x2, nint, Round & 1, w);
					}
				if (nround) {
					CPXinfodblparam(Env, CPX_PARAM_EPINT,
							w+1, w+2, w+3);
					if (w[0] <= intwarn_tol && !(Round & 8))
						nround = 0;
					else
						CPXgetdblparam(Env, CPX_PARAM_EPINT, w+1);
					}
				if (nround) {
					if ((solinfo[stat].code & ~1) == 2
					 && !(Round & 2))
						stat += 19;
					if (Round & 4)
						nround = 0;
					else if (!(Round & 1))
						nround = -nround;
					}
				}
			}
		}
	else
#endif
		{
		switch(CPXsolution(Env, cpx, &stat, &obj, x, y, 0, 0)) {
		 case 0: break;
#ifdef CPXERR_NO_DUAL_SOLN
		 case CPXERR_NO_DUAL_SOLN:
			y = 0;
			if (!CPXsolution(Env, cpx, &stat, &obj, x, y, 0, 0)) {
				if (method > 0 && stat == 1) {
					SI = solinfo1 + 1;
					goto have_SI;
					}
				break;
				}
#endif /*CPXERR_NO_DUAL_SOLN*/
		 default:
			stat = CPXgetstat(Env, cpx);
		 }
		stat = nos == 2 ? 2 : statadjust(cpx, stat0 = stat);
		}
	if (method < 0) {	/* dual */
		z = x;
		x = y;
		y = z;
		z1 = z + m;
		l = LUrhs;
		u = Urhsx;
		/* fix dual variables */
		for(le = l + m; l < le; l++, u++, z++)
			if (*u > *l && *l > negInfinity && *u < Infinity) {
				if (*z < -*z1)
					*z = *z1;
				z1++;
				}
		if (objsen > 0) {
			for(z = y, z1 = y + m; z < z1; z++)
				*z = -*z;
			for(z = x, z1 = x + n; z < z1; z++)
				*z = -*z;
			}
		}
 have_stat:
	SI = solinfo + stat;
 have_SI:
	i = Snprintf(hbuf, sizeof(hbuf), "%s: %s", Oinfo.bsname, SI->msg);
	if (!stat)
		i += Snprintf(hbuf+i, sizeof(hbuf)-i,
			": %s", failstat(stat0));
#ifdef CPLEX_MIP
	else if (stat == 16 || stat == 17) {
		CPXgetintparam(Env, CPX_PARAM_NODELIM, &nodelim);
		i += Snprintf(hbuf+i, sizeof(hbuf)-i,
			".\nCurrent node limit = %d", nodelim);
		}
#endif
	else if (!SI->wantobj)
		i += Snprintf(hbuf+i, sizeof(hbuf)-i, ".");
	if (SI->wantobj) {
		if (obj == Infinity)
			CPXgetobjval(Env, cpx, &obj);
		g_fmtop(buf, obj);
		i += Snprintf(hbuf+i, sizeof(hbuf)-i, "; objective %s", buf);
		}
	solve_result_num = SI->code;
	if (nosp) {
		if (nos < 0 || nos > 12)
			nos = 11;
		if (nos != 1)
			i += Snprintf(hbuf+i, sizeof(hbuf)-i,
				"\nnetopt found %s.", netmsg[nos]);
		if (net_nodes != lnc || net_arcs != nwv)
			i += Snprintf(hbuf+i, sizeof(hbuf)-i,
			"\nNetwork extractor found %d nodes and %d arcs.",
				net_nodes, net_arcs);
		i += Snprintf(hbuf+i, sizeof(hbuf)-i,
				"\n%d network simplex iterations.", netiters);
		}
	if (nint1 && !mipstat)
		i += Snprintf(hbuf+i, sizeof(hbuf)-i,
		 "\n%d MIP simplex iterations\n%d branch-and-bound nodes",
			ii, nodecnt);
	if (itc > 0 || itci > 0
	 || Optimize == CPXprimopt || Optimize == CPXdualopt) {
		i += Snprintf(hbuf+i, sizeof(hbuf)-i,
			"\n%d %s%ssimplex iterations (%d in phase I)%s",
				itc, nelqf ? "QP " : "", algname, itci,
				method > 0 ? "" : " on the dual problem");
		if (itc == 0 && itci == 0 && netiters <= 0
		 && solve_result_num == 0)
			solve_result_num = 4;
		}
#ifdef CPX_STAT_OPTIMAL	/* >= 8.0 */
	else if (Optimize == CPXlpopt && lpoptalg == CPX_ALG_SIFTING) {
		itc = CPXgetsiftitcnt(Env, cpx);
		itci = CPXgetsiftphase1cnt(Env, cpx);
		i += Snprintf(hbuf+i, sizeof(hbuf)-i,
			"\n%d sifting subproblems solved (%d in phase I)", itc, itci);
		}
#endif
	if (bitc > 0) {
		i += Snprintf(hbuf+i, sizeof(hbuf)-i, "\n%d %sbarrier iterations",
			bitc, baralgname);
		if (cbi.nx[0] | cbi.nx[1])
			i += Snprintf(hbuf+i, sizeof(hbuf)-i,
				"\n%d push, %d exchange %s crossover iterations",
				cbi.nx[0], cbi.nx[1], cbi.xkind);
		}
	if (cbi.np[0] + cbi.np[1] + cbi.np[3])
		i += Snprintf(hbuf+i, sizeof(hbuf)-i,
			"\nCPLEX's %spresolve eliminated totals of %d "
			"constraints and %d variables\n(perhaps repeatedly) "
			"and made %d coefficient changes.",
			nint1 ? "MIP " : "", cbi.np[0], cbi.np[1], cbi.np[3]);
	if (cbi.np[2])
		i += Snprintf(hbuf+i, sizeof(hbuf)-i,
			"\nCPLEX's aggregator made %d substitutions.", cbi.np[2]);
#ifdef CPLEX_MIP
	if (bestnode) {
		bbound = bobj0;
		if ((j = SI->code) < 200) {
			if (objsen*(bbound = bobj) > objsen*obj)
				bbound = obj;
			}
		else if (j < 400) {
			if (j >= 300)
				bbound = -bbound;
			}
		else if (j < 500 && objsen*bobj <= objsen*obj)
			bbound = bobj;
		if (n_obj > 0) {
			if ((bestnode & 3) == 3) {
				bb = (real*)M1zapalloc(2*sizeof(real)*n_obj);
				bn = bb + n_obj;
				}
			else
				bn = bb = (real*)M1zapalloc(n_obj*sizeof(real));
			if (bestnode & 1) {
				if (obj_no >= 0)
					bb[obj_no] = bbound;
				suf_rput("bestbound", ASL_Sufkind_obj, bb);
				}
			if (bestnode & 2) {
				if (obj_no >= 0)
					bn[obj_no] = bobj;
				suf_rput("bestnode", ASL_Sufkind_obj, bn);
				}
			}
		if (bestnode & 1)
			suf_rput("bestbound", ASL_Sufkind_prob, &bbound);
		if (bestnode & 2)
			suf_rput("bestnode", ASL_Sufkind_prob, &bobj);
		}
#if CPX_VERSION < 1100
	if (ccpr) {
		if (cl0 | cl | cov) {
			if (cl0)
				i += Snprintf(hbuf+i, sizeof(hbuf)-i,
			  "\n%d of %d clique inequalities used", cl, cl0);
			if (cov)
				i += Snprintf(hbuf+i, sizeof(hbuf)-i,
			  "\n%d cover cuts added", cov);
			}
		else
			i += Snprintf(hbuf+i, sizeof(hbuf)-i,
				"\nNo clique or cover cuts used.");
		}
#endif
	if (absmipgap > 0. && absmipgap < Infinity && !(retmipgap & 4))
			i += Snprintf(hbuf+i, sizeof(hbuf)-i,
				"\nabsmipgap = %g, relmipgap = %g",
				absmipgap, relmipgap);
	if (retmipgap & 1) {
		suf_rput("relmipgap", ASL_Sufkind_obj, &relmipgap);
		suf_rput("relmipgap", ASL_Sufkind_prob, &relmipgap);
		}
	if (retmipgap & 2) {
		suf_rput("absmipgap", ASL_Sufkind_obj, &absmipgap);
		suf_rput("absmipgap", ASL_Sufkind_prob, &absmipgap);
		}
#endif
	if (stat1 != 1) {
		SI = solinfo + stat1;
		i += Snprintf(hbuf+i, sizeof(hbuf)-i,
			"\nStatus recovering solution: %s", SI->msg);
		if (!stat1)
			i += Snprintf(hbuf+i, sizeof(hbuf)-i, ": stat = %d", stat10);
		if (solve_result_num < SI->code)
			solve_result_num = SI->code;
		}
	if (aggtries > 1)
		i += Snprintf(hbuf+i, sizeof(hbuf)-i,
			"\nTried aggregator %d times", aggtries);
	if (solve_result_num >= 200 && solve_result_num < 300) {
		i += send_iis(asl, cpx, hbuf+i);
		}
	if (send_statuses(asl, cpx, d))
		i += Snprintf(hbuf+i, sizeof(hbuf)-i, "\nNo basis.");
	else if (asl->i.flags & 1 && method > 0)
		switch(stat) {
			/* Note: stat values have been adjusted to accord with */
			/* earlier CPLEX versions -- don't use CPX_STAT_INFEASIBLE */
			/* or CPX_STAT_UNBOUNDED here... */
		  case 2: /*infeasible*/
			if (wantray(1,&cpx,&itc,&itci,nelqf)) {
				j = send_dray(asl,cpx,nelqf);
				i += Snprintf(hbuf+i, sizeof(hbuf)-i, dray_msg[j]);
				}
			goto j_check;
		  case 3: /*unbounded*/
			if (wantray(0,&cpx,&itc,&itci,0))
			  i += Snprintf(hbuf+i, sizeof(hbuf)-i, (j = send_ray(asl,cpx))
				? "\nfailed to compute variable.unbdd"
				: "\nvariable.unbdd returned");
 j_check:
			*pcpx = cpx;
			if (itc) {
				i += Snprintf(hbuf+i, sizeof(hbuf)-i,
				  "\n%d extra %ssimplex iterations for ray",
					itc, algname);
				if (itci)
				  i += Snprintf(hbuf+i, sizeof(hbuf)-i,
					" (%d in phase I)",itci);
				}
			if (!j)
				solve_result_num += 10;
		  }
	if (basis_cond) {
		if (CPXgetdblquality(Env, cpx, &bcond, CPX_KAPPA)) {
			i += Snprintf(hbuf+i, sizeof(hbuf)-i,
				"\nBasis condition is unavailable.");
			bcond = 0;
			}
		else
			i += Snprintf(hbuf+i, sizeof(hbuf)-i,
				"\nBasis condition = %g", bcond);
		if (n_obj > 0)
			suf_rput("basis_cond", ASL_Sufkind_obj, &bcond);
		suf_rput("basis_cond", ASL_Sufkind_prob, &bcond);
		}
#ifdef CPX_PARAM_MIPKAPPASTATS
	if (haveqmet) {
		i += Snprintf(hbuf+i, sizeof(hbuf)-i,
			"\n%s = %.3g", QualInfo[0].desc, qmet[0]);
		for(j = 1; j < haveqmet; ++j)
			if (qmet[j] > 0.)
				i += Snprintf(hbuf+i, sizeof(hbuf)-i,
					"\n%s = %.3g", QualInfo[j].desc, qmet[j]);
		}
#endif
	if (nround) {
		wb = "";
		if (nround < 0) {
			nround = -nround;
			wb = "would be ";
			}
		i += Snprintf(hbuf+i, sizeof(hbuf)-i,
			"\n%d integer variables %srounded (maxerr = %g).",
			nround, wb, w[0]);
		if (w[0] > intwarn_tol && 0.5*w[0] > w[2]) {
			i += Snprintf(hbuf+i, sizeof(hbuf)-i,
			  "\nAssigning integrality = %.1g might help.", 0.5*w[0]);
			i += Snprintf(hbuf+i, sizeof(hbuf)-i,
			  "\nCurrently integrality = %g.", w[1]);
			}
		}
#if CPX_VERSION >= 1100
	if (npt != opt)
		CPXchgprobtype(Env, cpx, opt);
	if (cutstats && nint1)
	    for(CI = Cut_Info; CI->cutname; ++CI)
		if (!CPXgetnumcuts(Env, cpx, CI->cuttype, &j) && j > 0)
			i += Snprintf(hbuf+i, sizeof(hbuf)-i,
				"\n%d %s cut%s", j, CI->cutname, "s" + (j == 1));
#endif
#ifdef CPX_PARAM_POPULATELIM
	if (Optimize == CPXmipopt && poolstub && x)
		poolwrite(asl, cpx, d, nelqf, hbuf, &i, sizeof(hbuf));
	if (pretunefile)
		tunewrite(pretunefile, "pretunefile", hbuf, &i, sizeof(hbuf));
#endif
#if CPX_VERSION >= 1000 /*{*/
	if (pretunefileprm)
		tunewriteprm(pretunefileprm, "pretunefileprm", hbuf, &i, sizeof(hbuf));
#ifdef CPX_TUNE_TILIM /* CPLEX 11 */
	if (tunefile || tunefileprm) {
		tunerun(cpx, hbuf, &i, sizeof(hbuf));
		if (tunefile)
			tunewrite(tunefile, "tunefile", hbuf, &i, sizeof(hbuf));
		if (tunefileprm)
			tunewriteprm(tunefileprm, "tunefileprm", hbuf, &i, sizeof(hbuf));
		}
#endif
#endif /*}*/
	write_sol(hbuf, x, y, &Oinfo);
	if (Logf) {
		CPXsetlogfile(Env, NULL);
		CPXfclose(Logf);
		Logf = 0;
		}
	}

 static void
show_times(void)
{
	int i;

	Times[3] = xectim_();
	for(i = 1; i <= 2; ++i)
	    if (time_flag & i) {
		fprintf(i == 1 ? stdout : Stderr,
		"\nTimes (seconds):\nInput =  %g\nSolve =  %g\nOutput = %g\n",
			Times[1] - Times[0], Times[2] - Times[1],
			Times[3] - Times[2]);
		}
#if CPX_VERSION >= 12050000
	if (!DTimes_failed)
		for(i = 4; i <= 8; i += 4)
		    if (time_flag & i) {
		fprintf(i == 4 ? stdout : Stderr,
		"\nTimes (ticks):\nInput =  %g\nSolve =  %g\nOutput = %g\n",
			DTimes[1] - DTimes[0], DTimes[2] - DTimes[1],
			DTimes[3] - DTimes[2]);
		}
	if (time_flag & 0x30 && num_cores) {
		for(i = 16; i <= 32; i += 16)
			if (time_flag & i)
				fprintf(i == 16 ? stdout : stderr,
					"\n%d logical cores are available.\n", num_cores);
		}
#endif
	}

 static void
basread(ASL *asl, cpxlp *cpx)
{
	if (CPXreadcopybase(Env, cpx, startbas)) {
		printf("Could not read starting basis file \"%s\"\n",
			startbas);
		need_nl = 0;
		}
	}

 static void CPXPUBLIC
mymsgfunc(void *handle, Const char *msg)
{
	static int wantnl;
	Const char *msg1;

	handle = handle; /* shut up non-use warning */
	if (Logf)
		CPXfputs(msg, Logf);
	if (!strncmp(msg, "Tried aggregator ", 17))
		aggtries = atoi(msg+17);
	else if (prestats) {
		if (!netopting || !strncmp(msg,"Extract",7))
			goto print_msg;
		}
	else if (!netopting) {
		msg1 = msg;
		switch(*msg) {
		  case 'M':
			if (msg[1] != 'I')
				break;
			msg1++;
			/* no break */
		  case 'L':
		  case 'Q':
			if (msg1[1] == 'P' && msg1[2] == ' ')
				msg1 += 3;
		  }
		switch(*msg1) {
		  case 'A':
			if (strncmp(msg,"Aggregator",10))
				goto accept;
			break;
		  case 'N':
			if (strncmp(msg,"No presolve",11))
				goto accept;
			break;
		  case 'P':
			if (!strncmp(msg,"Presolve",8)) {
				wantnl = !strchr(msg+8,'\n');
				break;
				}
			/* no break */
		  default:
		  accept:
			if (wantnl)
				wantnl = !strchr(msg,'\n');
			else {
 print_msg:
				printf("%s", msg);
				fflush(stdout); /* in case AMPL is */
						/* writing a log file */
				}
		  }
		}
	}

 static char iis_table[] = "\n\
0	non	not in the iis\n\
1	low	at lower bound\n\
2	fix	fixed\n\
3	upp	at upper bound\n\
4	mem	member\n\
5	pmem	possible member\n\
6	plow	possibly at lower bound\n\
7	pupp	possibly at upper bound\n\
8	bug\n";

 static SufDecl
suftab[] = {
	{ "absmipgap", 0, ASL_Sufkind_obj  | ASL_Sufkind_outonly },
	{ "absmipgap", 0, ASL_Sufkind_prob  | ASL_Sufkind_outonly },
	{ "basis_cond", 0, ASL_Sufkind_obj  | ASL_Sufkind_outonly },
	{ "basis_cond", 0, ASL_Sufkind_prob  | ASL_Sufkind_outonly },
	{ "bestbound", 0, ASL_Sufkind_obj  | ASL_Sufkind_outonly },
	{ "bestbound", 0, ASL_Sufkind_prob | ASL_Sufkind_outonly },
	{ "bestnode", 0, ASL_Sufkind_obj  | ASL_Sufkind_outonly },
	{ "bestnode", 0, ASL_Sufkind_prob | ASL_Sufkind_outonly },
	{ "current", 0, ASL_Sufkind_con | ASL_Sufkind_outonly },
	{ "current", 0, ASL_Sufkind_var | ASL_Sufkind_outonly },
	{ "direction", 0, ASL_Sufkind_var },
	{ "down", 0, ASL_Sufkind_con | ASL_Sufkind_outonly },
	{ "down", 0, ASL_Sufkind_var | ASL_Sufkind_outonly },
	{ "dunbdd", 0, ASL_Sufkind_con | ASL_Sufkind_outonly },
	{ "iis", iis_table, ASL_Sufkind_var | ASL_Sufkind_outonly },
	{ "iis", 0, ASL_Sufkind_con | ASL_Sufkind_outonly },
#ifdef CPX_PARAM_FEASOPTMODE /* >= 9.2b */
	{ "lazy", 0, ASL_Sufkind_con },
#endif
#ifdef CPX_PARAM_POPULATELIM
	{ "npool", 0, ASL_Sufkind_prob | ASL_Sufkind_outonly},
#endif
	{ "priority", 0, ASL_Sufkind_var },
	{ "ref", 0, ASL_Sufkind_var | ASL_Sufkind_real },
	{ "relmipgap", 0, ASL_Sufkind_obj  | ASL_Sufkind_outonly },
	{ "relmipgap", 0, ASL_Sufkind_prob  | ASL_Sufkind_outonly },
	{ "sos", 0, ASL_Sufkind_var },
	{ "sos", 0, ASL_Sufkind_con },
	{ "sosno", 0, ASL_Sufkind_var | ASL_Sufkind_real },
	{ "sosref", 0, ASL_Sufkind_var | ASL_Sufkind_real },
	{ "sstatus", 0, ASL_Sufkind_var, 1 },
	{ "sstatus", 0, ASL_Sufkind_con, 1 },
	{ "unbdd", 0, ASL_Sufkind_var | ASL_Sufkind_outonly},
	{ "up", 0, ASL_Sufkind_con | ASL_Sufkind_outonly },
	{ "up", 0, ASL_Sufkind_var | ASL_Sufkind_outonly }
	};

 static void
adjust_version(ASL *asl)
{
	char *t;
	const char *s;
	int n;

	s = CPXversion(Env);
	n = 2*strlen(s) + 19;
	Oinfo.bsname = M1alloc(n);
	t = strcpy1(Oinfo.bsname, "CPLEX ");
	t = strcpy1(t, s);
	Oinfo.version = ++t;
	t = strcpy1(t, "AMPL/");
	t = strcpy1(t, Oinfo.bsname);
	}

 static int
Optimize1(CPXENVptr e,  cpxlp *c)
{
	int k, rc;

#ifdef CPX_PARAM_POPULATELIM
	if (Optimize == CPXmipopt) {
		if (populate == 2 && poolstub)
			rc = CPXpopulate(e,c);
		else {
			if (!poolstub)
				CPXsetintparam(Env, CPX_PARAM_SOLNPOOLCAPACITY, 0);
			rc = CPXmipopt(e,c);
			}
		}
	else
#endif
		rc = Optimize(e, c);
	if (rc == CPXERR_PRESLV_INForUNBD
#ifdef CPX_STAT_INForUNBD	/* >= 8.0 */
		|| CPXgetstat(e,c) == CPX_STAT_INForUNBD
#endif
		) {
		CPXsetintparam(e, CPX_PARAM_AGGIND, 0);
		CPXsetintparam(e, CPX_PARAM_PREIND, 0);
		rc = Optimize(e, c);
		}
#ifdef CPX_PARAM_FEASOPTMODE /* >= 9.2b */
	if (want_feasopt)
		switch(k = CPXgetstat(e,c)) {
		  case CPX_STAT_INFEASIBLE:
		  case CPX_INTEGER_INFEASIBLE:
		  case CPX_STAT_INForUNBD:
		  case CPXMIP_INForUNBD:
			rc = feasopt(e, c);
		  }
#endif /*CPX_PARAM_FEASOPTMODE*/
	return rc;
	}

 int
main(int argc, char **argv)
{
	int	nelqf, nint1, nos, *nosp, rc, status, z;
	dims	d;
	char	*solmsg, *stub;
	CPXCHANNELptr cpxresults;
	size_t L;
	 /* static values for use with longjmp */
	static cpxlp *cpx;
	static sig_func_type *oic;
	static FILE *nl;

	Times[0] = xectim_();
	Oinfo.bsname = basename(argv[0]); /* in case CPXopenCPLEX fails */

#ifdef BARRIER
	asl = ASL_alloc(ASL_read_fg);
#else
	asl = ASL_alloc(ASL_read_f);
#endif

	mdbl_val[0].U = Infinity;
	want_deriv = 0;
	Env = 0;
	cpx = 0;
	breaking = 0;
	oic = 0;
	nl = 0;
	rc = setjmp(Jb);
	if (rc) {
		if (nl)
			fclose(nl);
		--rc;
		if (solve_result_num > 0){
			if (amplflag | (Oinfo.wantsol & 1))
				rc = 0;
			goto ws_now;
			}
		goto done;
		}
	if (argc < 2)
#ifdef Stand_alone
		return cpxmain(argc,argv);
#else
		usage_ASL(&Oinfo, 1);
#endif
#ifndef KEEP_BANNER
	/* keep banner only if first arg is -v */
	if (argc > 1 && *(stub = argv[1]) == '-' && stub[1] == 'v')
		nos = -1;
	else {
		nos = dup(2);
		close(2);
		z = open("/dev/null",O_WRONLY); /* should return 2 */
		}
#endif
#ifdef LICENSE_FILE
	if (!(stub = getenv("ILOG_LICENSE_FILE")) || !*stub)
		putenv("ILOG_LICENSE_FILE=" LICENSE_FILE);
#endif
	if ((Env = CPXopenCPLEX(&z))) {
#if CPX_VERSION >= 12050000
		DTimes_failed = CPXgetdettime(Env, &DTimes[0]);
#endif
		adjust_version(asl);
		}
#ifndef KEEP_BANNER
	if (nos >= 0) {
		close(2);
		dup(nos);	/* should dup onto 2 */
		close(nos);
		}
#endif

	/* prepare to return error messages in .sol file */
	if (!(stub = getstub(&argv, &Oinfo)))
		usage_ASL(&Oinfo, 1);
	nl = jac0dim(stub, (fint)0);
#ifdef CPXERR_IN_INFOCALLBACK /* CPLEX 11 */
	CPXsetterminate(Env, &breaknow);
#endif
	oic = signal(SIGINT, intcatch);

	if (!Env) {
		if (amplflag) {
			badlic(570, z);
			goto ws_now;
			}
		badlic(0, z);
		goto done;
		}

	if ((z = CPXgetchannels(Env, &cpxresults, 0,0,0))) {
		badret("CPXgetchannels", z, 531);
		goto ws_now;
		}

	suf_declare(suftab, sizeof(suftab)/sizeof(SufDecl));
	amplin(asl, &cpx, &nl, &d, &nelqf, &nint1, argv);
	if (!cpx) {
		if (solve_result_num > 0) {
 ws_now:
			L = strlen(Oinfo.bsname) + strlen((char*)asl->i.uinfo);
			sprintf(solmsg = (char*)M1alloc(L+3), "%s: %s",
				Oinfo.bsname, (char*)asl->i.uinfo);
			write_sol(solmsg, 0, 0, &Oinfo);
			}
		else {
			printf("Failed to load the problem!\n");
			rc = 1;
			}
		goto done;
		}
	if (startsol && CPXreadcopysol(Env, cpx, startsol))
		printf("Failed to read startsol file \"%s\".\n", startsol);
	if (wrtfname)
		CPXwriteprob(Env, cpx, wrtfname, NULL);
	if (startbas)
		basread(asl, cpx);
#ifdef BARRIER
#ifndef NO_DEPRECATED
	if (startvec && (status = CPXreadcopyvec(Env, cpx, startvec)))
		printf("\n*** return %d from CPXreadcopyvec.\n", status);
#endif
#endif
#ifdef CPX_PARAM_IISIND /* version < 9.2b */
#ifdef CPLEX_MIP
	if (nint1 && starttree)
		treeio(cpx, starttree, "read", CPXreadcopytree);
#endif
#endif
	fflush(stdout);
	Times[1] = xectim_();
#if CPX_VERSION >= 12050000
	if (!DTimes_failed)
		DTimes_failed = CPXgetdettime(Env, &DTimes[1]);
#endif
	nosp = 0;

	disconnectchannel(cpxresults);
	addfuncdest(cpxresults, 0, mymsgfunc);
	if (use_netopt) {
		nos = netopt(cpx,&d,&net_status,&net_nodes,&net_arcs,&netiters);
		if (nos)
			net_status = 0;
		nosp = &net_status;
		}
	if (prestats && Optimize == CPXmipopt) {
		CPXpresolve(Env, cpx, CPX_ALG_NONE);
		/* Collect presolve statistics. */
		lpcbf(Env, 0, -1, 0);
		/* Turn off call-backs to allow use of multiple threads in CPXmipopt. */
		CPXsetmipcallbackfunc(Env, 0, 0);
		CPXsetlpcallbackfunc(Env, 0, 0);
		}
	status = Optimize1(Env,cpx);
	Times[2] = xectim_();
#if CPX_VERSION >= 12050000
	if (!DTimes_failed)
		DTimes_failed = CPXgetdettime(Env, &DTimes[2]);
#endif
	breaking = 1;	/* reset in case, e.g., amplout runs CPXprimopt */
	if (endbas && !nint1)
		write_basis(cpx);
	if (endsol && CPXsolwrite(Env, cpx, endsol))
		printf("Failed to write endsol file \"%s\"\n", endsol);
#ifdef CPX_PARAM_IISIND /* version < 9.2b */
#ifdef CPLEX_MIP
	if (nint1 && endtree)
		treeio(cpx, endtree, "write", (Treeio)CPXtreewrite);
#endif
#endif
	amplout(asl, &cpx, &d, status, nelqf, nint1, nosp);
 done:
	if (oic)
		signal(SIGINT, oic);
#if CPX_VERSION >= 12050000
	if (!DTimes_failed)
		DTimes_failed = CPXgetdettime(Env, &DTimes[3]);
#endif
	if (cpx)
		CPXfreeprob(Env, &cpx);
	if (Env)
		CPXcloseCPLEX(&Env);
	ASL_free(&asl);
	show_times();
	return rc;
	}

 void
mainexit_ASL(int rc)
{ longjmp(Jb, rc+1); }
