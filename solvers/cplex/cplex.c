/*******************************************************************
Copyright (C) 2019 AMPL Optimization, Inc.; written by David M. Gay.

Permission to use, copy, modify, and distribute this software and its
documentation for any purpose and without fee is hereby granted,
provided that the above copyright notice appear in all copies and that
both that the copyright notice and this permission notice and warranty
disclaimer appear in supporting documentation.

The author and AMPL Optimization, Inc. disclaim all warranties with
regard to this software, including all implied warranties of
merchantability and fitness.  In no event shall the author be liable
for any special, indirect or consequential damages or any damages
whatsoever resulting from loss of use, data or profits, whether in an
action of contract, negligence or other tortious action, arising out
of or in connection with the use or performance of this software.
*******************************************************************/

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

#if CPX_APIMODEL == CPX_APIMODEL_LARGE
#include "cplexx.h"
#endif
#include "cplex.h"

#if CPX_VERSION < 12060200
#define OBJ_ADJ
#endif

#if CPX_VERSION < 12060000
#undef Want_Distmipopt
#else
#undef filename
#ifdef Want_Distmipopt
#include "cplexdistmip.h"
#endif
#endif

#if CPX_VERSION < 12090000
#undef NO_MOkwf
#define NO_MOkwf
#endif

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
#ifndef CPX_PARAM_QPMETHOD
#define CPX_PARAM_PDSWITCH 1063
#endif
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

#if CPX_APIMODEL == CPX_APIMODEL_LARGE
#define ForceZ | ASL_use_Z
	typedef size_t CStype;
#undef nzc
#undef A_colstarts
#define nzc asl->i.nZc_
#define A_colstarts asl->i.A_colstartsZ_
#define CPXcopylpwnames(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17) CPXLcopylpwnames(a1,a2,a3,a4,a5,a6,a7,a8,(CPXLONG*)a9,a10,a11,a12,a13,a14,a15,(char const *const *)(a16),(char const *const *)(a17))
#define CPXaddlazyconstraints(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10) CPXLaddlazyconstraints(a1,a2,a3,a4,a5,a6,(CPXLONG*)a7,a8,a9,(char const *const *)(a10))
#define CPXaddusercuts(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10) CPXLaddusercuts(a1,a2,a3,a4,a5,a6,(CPXLONG*)a7,a8,a9,(char const *const *)(a10))
#else
	typedef int CStype;
#define ForceZ /*nothing*/
#endif

static ASL *asl;
 static double Times[4];
#if CPX_VERSION >= 12050000
 static double DTimes[4];
 static int DTimes_failed, num_cores;
#endif

 typedef struct
PBuf { char *s, *se; } PBuf;

 static void
Bpf(PBuf *b, const char *fmt, ...)
{
	va_list ap;

	if (b->s < b->se) {
		va_start(ap, fmt);
		b->s += Vsnprintf(b->s, b->se - b->s, fmt, ap);
		va_end(ap);
		}
	}


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
	set_incompat	= 27,
	set_objrep	= 28,
	set_qcdual	= 29,
	set_multiobj	= 30
	};
#ifdef CPX_PARAM_FEASOPTMODE /* >= 9.2b */
#define Uselazy
#endif

 enum { /* sf_mdbl f values */
	set_dual_ratio = 0,
	set_droptol = 1
	};

 static mint_values
mint_val[31] = {
	/* set_crossover */	{0, 2, 1},
	/* set_dualthresh */	{-0x7fffffff, 0x7fffffff, 0},
	/* set_netopt */	{0, 3, 1},
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
	/* set_mipbasis */	{-1, 1, -1},
	/* set_basis_cond */	{0, 1, 0},
	/* set_retmipgap */	{0, 7, 0},
	/* set_feasopt */	{0, 2, 0},
	/* set_feasoptobj */	{1, 3, 1},
	/* set_lazy */		{0, 3, 3},
	/* set_populate */	{0, 2, 0},
	/* set_pooldual */	{0, 1, 0},
	/* set_resolve */	{0, 1, 1},
	/* set_cutstats */	{0, 1, 0},
	/* set_incompat */	{0, 2, 1},
	/* set_objrep */	{0, 3, 2},
	/* set_qcdual */	{0, 1, 1},
	/* set_multiobj */	{0, 1, 0}
	};

 static mdbl_values
mdbl_val[] = {
	/* set_dual_ratio */	{1., 1e30, 3.},
	/* set_droptol */	{0., 1e30, 0.}
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
#define objrep		mint_val[28].val
#define want_qcdual	mint_val[29].val
#define multiobj	mint_val[30].val
#define dual_ratio	mdbl_val[0].val
#define droptol		mdbl_val[1].val

 static int hybmethod = CPX_ALG_PRIMAL;
 static int netiters = -1;
#if CPX_VERSION >= 12080000
#define USE_CHANNELS
#else
#undef USE_CHANNELS
#endif
#ifdef USE_CHANNELS
 static CPXCHANNELptr cpxlog;
#else
 static CPXFILEptr Logf;
#endif
 static char cplex_version[] = "AMPL/CPLEX with bad license\0\nAMPL/CPLEX Driver Version 20190501\n";
 static char *baralgname, *endbas, *endsol, *endtree, *endvec, *logfname;
 static char *paramfile, *poolstub, *pretunefile, *pretunefileprm;
 static char *startbas, *startsol, *starttree, *startvec, *tunefile, *tunefileprm;
 static char *tunefix, *tunefixfile, *vmconfig, *workfiledir, *wrtfname, *wrtmipstart;
 static int bestnode, breaking, costsens, lpoptalg, mbas, method;
 static int nbas, netopting, nosolve, objadj, objsen, relax, zap_lpcbf;
 static int aggtries, net_status, net_nodes, net_arcs;
 static double obj_adj;
 static char *algname = "";
 static real intwarn_tol = 1e-9;
 static real qcdmax = 1e9;
#ifdef CPXERR_QCP_SENSE
 static real qctol1 = 1e-5, qctol2 = 1e-5, qctol3 = 1e-5;
#endif

 typedef struct { char *msg; int code; int wantobj; }
Sol_info;
 static Sol_info solinfo[] = {
	 { "unrecoverable failure", 500, 0 },
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
	 { "optimal relaxed quadratic penalty in feasopt", 126, 1}, /* --> 62 */
	 { "dettimelim reached with feasible solution", 440, 1},
	 { "dettimelim reached with infeasible solution", 441, 0}
#ifdef CPX_PARAM_FEASOPTMODE /* >= 9.2b {*/
	,{ "minimal sum of constraint relaxations", 230, 1 }, /* --> 65 */
	 { "optimal among minimal sum of constraint relaxations", 231, 1 },
	 { "minimal number of constraint relaxations", 232, 1 },
	 { "optimal among minimal number of constraint relaxations", 233, 1 },
	 { "minimal sum of squares of constraint relaxations", 234, 1 },
	 { "optimal among minimal sum of squares of constraint relaxations", 235, 1 },
	 { "integer minimal sum of constraint relaxations", 240, 1 },
	 { "integer optimal among minimal sum of constraint relaxations", 241, 1 },
	 { "integer minimal number of constraint relaxations", 242, 1 },
	 { "integer optimal among minimal number of constraint relaxations", 243, 1 },
	 { "integer minimal sum of squares of constraint relaxations", 244, 1 },
	 { "integer optimal among minimal sum of squares of constraint relaxations", 245, 1 }
#ifdef CPX_STAT_MULTIOBJ_OPTIMAL /*{*/
	,{ "multiobjective optimal", 6, 1 }, /*77*/
	 { "multiobjective infeasible", 211,  0 },
	 { "multiobjective infeasible or unbounded", 212,  0 },
	 { "multiobjective unbounded", 311, 0 },
	 { "multiobjective nonoptimal", 140, 1 },
	 { "multiobjective stopped", 450, 1 }
#endif /*}*/
#endif /*}*/
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
	/* for linadj and qcduals */
	cde	*consave;
	cde	*objsave;
	ograd	**ogsave;
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
badretfmt(int rc, const char *fmt, ...)
{
	va_list ap;
	char buf[4200], *s;
	int k;

	va_start(ap, fmt);
	k = Vsnprintf(buf, sizeof(buf), fmt, ap) + 1;
	if (rc) {
		solve_result_num = rc;
		memcpy(s = (char*)M1alloc(k), buf, k);
		asl->i.uinfo = s;
		}
	else
		fprintf(Stderr, "%s\n", buf);
	}

 static void
badret(const char *what, int i, int rc)
{
	badretfmt(rc, "%s failed; error code %d.", what, i);
	if (rc)
		exit(1);
	}

 static void
badwrite(int rc, const char *fmt, const char *fname, int status)
{
	badretfmt(rc, fmt, fname, status);
	exit(1);
	}

#ifdef CPX_CON_INDICATOR /* >= 9.2b */
#define ALLOW_CLP ASL_allow_CLP

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

	rv = CPXbaropt(e, c);
#if !defined(NO_DEPRECATED) && CPX_VERSION < 12070000
	if (!rv) {
		int s;
		if (endvec && (s = CPXvecwrite(e, c, endvec)))
			printf("\n*** return %d from CPXvecwrite.\n", s);
		}
#endif
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
startcomp(ASL *asl, int n0, int m, int nextra, CStype *ka, int *kal, int *ia,
	double *a, double *b, double *c,
	double **cdualp, double **cprimp, double **rdualp, double **rprimp)
{
	CStype i, j, k, n;
	double *r, t, *x, *y;

	/* Why can't CPXcopystart do this calculation? */

	if (!X0 && !pi0)
		return 0;
	n = n0 - objadj;
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
	set_bestbound	= 13,
	set_nosolve	= 14,
	set_benders	= 15
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

#ifdef CPXPARAM_Benders_Strategy
	 case set_benders:
		Optimize = CPXbendersopt;
		break;
#endif

	 case set_nosolve:
		nosolve = 1;
	 }
	return v;
	}

 static void
baddval(Option_Info *oi, keyword *kw, double t, double L, double U)
{
	printf("Rejecting %s %g; must be between %g and %g\n",
		kw->name, t, L, U);
	badopt_ASL(oi);
	}

 static void
badival(Option_Info *oi, keyword *kw, int t, int L, int U)
{
	printf("Rejecting %s %d; must be between %d and %d\n",
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
	switch(i) {
#ifdef BARRIER
	  case set_crossover:
		set_baropt();
		break;
#endif
	  case set_netopt:
		if (t == 3) {
			m->val = 0;
			CPXsetintparam(Env, CPX_PARAM_LPMETHOD, CPX_ALG_NET);
			}
		break;
#ifdef CPLEX_MIP
	  case set_mipcuts:
		{static int op[] = {
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
		break;
#endif
		}
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
	const char *fmt, *what = kw->name;

	f = (int)strtol(v, &rv, 10);
	if (rv == v || (*rv != '=' && *rv != ' ')) {
		printf("Expected an integer parameter number after %s%s, not \"%s\"\n",
			what, oi->eqsign, v);
		badopt_ASL(oi);
		while(*v > ' ' && *v++ != ',');
		return v;
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
		rv = v + 1;
		if (CPXgetintparam(Env, f, &t))
			goto badopt;
		printf("%s=%d=%d\n", what, f, t);
		oi->option_echo &= ~ASL_OI_echothis;
		return rv;
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
		if (CPXinfointparam(Env, f, z, z+1, z+2)) {
 badopt:
			fmt = "Rejecting %s=%d=%d; bad iparam number %d.\n";
			z[1] = f;
			}
		else
			fmt = "Rejecting %s=%d=%d; assigned value must be between %d and %d\n";
		printf(fmt, what, f, t, z[1], z[2]);
		badopt_ASL(oi);
		}
	return rv;
	}

 static char *
sf_spar(Option_Info *oi, keyword *kw, char *v)
{
	int c, f, q;
	char *b, *be, buf[CPX_STR_PARAM_MAX], buf1[CPX_STR_PARAM_MAX], *rv, *v0 = v;
	const char *what = kw->name;

	f = (int)strtol(v, &rv, 10);
	if (rv == v || (*rv != '=' && *rv != ' ')) {
		printf("Expected an integer parameter number after %s%s, not \"%s\".\n",
			what, oi->eqsign, v);
 bad:
		badopt_ASL(oi);
		while(*v > ' ' && *v++ != ',');
		return v;
		}
	for(v = rv; *v <= ' '; v++)
		if (!*v) {
			printf("Expected an integer value after %s%s.\n", what, oi->eqsign);
			goto bad;
			}
	if (*v == '=') {
		while(*++v <= ' ')
			if (!*v) {
				printf("Expected an integer value after %s%s%s.\n",
					what, oi->eqsign, v0);
				goto bad;
				}
		}
	if (*v == '?' && v[1] <= ' ') {
		++v;
		if (!CPXgetstrparam(Env, f, buf))
			printf("%s=%d=\"%s\"\n", what, f, buf);
		else {
 badparnum:
			printf("Bad string parameter number %d.\n", f);
 badpar:
			badopt_ASL(oi);
			}
		oi->option_echo &= ~ASL_OI_echothis;
		return v;
		}
	b = buf;
	be = buf + sizeof(buf);
	if ((q = *v) == '"' || q == '\'') {
		while((c = *++v)) {
			if (c == q) {
				++v;
				break;
				}
			if (b < be)
				*b++ = c;
			}
		}
	else do {
		if (b < be)
			*b++ = q;
		} while((q = *++v) > ' ');
	if (b < be) {
		*b = 0;
		if (CPXsetstrparam(Env, f, buf)) {
			if (CPXgetstrparam(Env, f, buf1))
				goto badparnum;
			printf("Bad value \"%s\" for %s=%d.\n", buf, what, f);
			goto badpar;
			}
		}
	else {
		printf("Oversize value for %s=%d.\n", what, f);
		badopt_ASL(oi);
		}
	return v;
	}

 static char *
sf_str(Option_Info *oi, keyword *kw, char *v)
{
	int c, f, q;
	char *b, *be, buf[CPX_STR_PARAM_MAX];
	const char *what = kw->name;

	f = Intcast kw->info;
	if (*v == '?' && v[1] <= ' ') {
		++v;
		buf[0] = 0;
		CPXgetstrparam(Env, f, buf);
		printf("%s=\"%s\"\n", what, buf);
		oi->option_echo &= ~ASL_OI_echothis;
		return v;
		}
	b = buf;
	be = buf + sizeof(buf);
	if ((q = *v) == '"' || q == '\'') {
		while((c = *++v)) {
			if (c == q) {
				++v;
				break;
				}
			if (b < be)
				*b++ = c;
			}
		}
	else do {
		if (b < be)
			*b++ = q;
		} while((q = *++v) > ' ');
	if (b < be) {
		*b = 0;
		if (CPXsetstrparam(Env, f, buf)) {
			printf("Bad value \"%s\" for %s.\n", buf, what);
			badopt_ASL(oi);
			}
		}
	else {
		printf("Oversize value for %s.\n", what);
		badopt_ASL(oi);
		}
	return v;
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
		printf("Rejecting %s %g; must be between %g and %g\n",
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
		while(*v > ' ' && *v++ != ',');
		return v;
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
		printf("Rejecting %s=%d=%g; assigned value must be between %g and %g\n",
			what, f, t, z[1], z[2]);
		badopt_ASL(oi);
		}
	return rv;
	}


 static char **file_name[22] = { &endbas, &endtree, &startbas, &starttree,
				&startsol, &endsol, &logfname, &wrtfname,
				&workfiledir, &poolstub, &paramfile,
				&pretunefile, &pretunefileprm, &tunefile,
				&tunefileprm, &tunefix, &tunefixfile, &startvec,
				&endvec, &vmconfig, &wrtmipstart };

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
	set_endvector	= 18,
	set_vmconfig	= 19,
	set_wrtmipstart = 20
	};

 static char *
sf_char(Option_Info *oi, keyword *kw, char *v)
{
	char *rv, *s, *t;
	const char *fmt;
	int c, f, q;
	size_t n;

	if (!*v) {
		printf("Rejecting %s: no following file name\n", kw->name);
		badopt_ASL(oi);
		return v;
		}
	f = Intcast kw->info;
	if (*v == '?' && v[1] <= ' ') {
		fmt = "%s=%s\n";
		if ((t = *file_name[f])) {
			for(s = t; (c = *s); ++s)
				if (c <= ' ')
					goto qfmt;
			}
		else {
			t = "";
 qfmt:
			fmt = "%s=\"%s\"\n";
			}
		printf(fmt, kw->name, t);
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
#ifdef USE_CHANNELS
		if (CPXsetlogfilename(Env, logfname, "w"))
#else
		if (!(Logf = CPXfopen(t,"w")))
#endif
			{
			printf("Cannot open logfile \"%s\"\n", t);
			logfname = 0;
			badopt_ASL(oi);
			}
		else
#ifdef USE_CHANNELS
			if (CPXgetchannels(Env, NULL, NULL, NULL, &cpxlog)) {
				printf("Cannot get cpxlog channel; ignoring logfile.\n");
				CPXsetlogfilename(Env, NULL, NULL);
				logfname = 0;
				}
#else
			CPXsetlogfile(Env, Logf);
#endif
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

#ifdef NO_MOkwf
#define MOkwf 0
#else /*{*/
 typedef struct
MOkwval {
	struct MOkwval *next;
	keyword *kw;
	char *val;
	int Objno;
	} MOkwval;

 static MOkwval *first_MOkw, **pnext_MOkw = &first_MOkw;
 static int *indobj, *mopriority, nmo, nmopri, *objind, *pri0;
 static real *moabstol, *momem, *moreltol, *moweight;
 static struct paramset **mopars;

 static int CPXPUBLIC
Optimizemo(CPXCENVptr env, cpxlp *cpx)
{
	int i, j, k, no, nz, ono, *vn;
	ograd *og;
	real *c;

	if ((i = CPXsetnumobjs(Env, cpx, no = nmo)))
		badret("CPXsetnumobjs", i, 1);
	ono = obj_no;
	if ((i = objind[ono]) > 0) {
		j = indobj[i];
		k = indobj[0];
		indobj[0] = j;
		indobj[i] = k;
		objind[j] = 0;
		objind[k] = i;
		}
	c = (real*)M1alloc(n_var*(sizeof(real) + sizeof(int)));
	vn = (int*)(c + n_var);
	i = 0;
	if (indobj[0] == ono) {
#ifndef NAN
		/* CPX_NO_OFFSET_CHANGE is #defined as NAN, so we construct a NaN */
		/* using constants in arith.h. */
		union { double d; unsigned int x[2]; } u;
		u.x[0] = QNaN0;
		u.x[1] = QNaN1;
#undef CPX_NO_OFFSET_CHANGE
#define CPX_NO_OFFSET_CHANGE u.d
#endif
		if ((k = CPXmultiobjchgattribs(env, cpx, 0, CPX_NO_OFFSET_CHANGE,
				moweight[ono], mopriority[ono], moabstol[ono],
				moreltol[ono], obj_name(ono))))
			badret("CPXmultiobjchgattribs", k, 1);
		i = 1;
		}
	for(; i < no; ++i) {
		j = indobj[i];
		for(nz = 0, og = Ograd[j]; og; og = og->next) {
			c[nz] = og->coef;
			vn[nz++] = og->varno;
			}
		if ((k = CPXmultiobjsetobj(env, cpx, i, nz, vn, c, objconst(j), moweight[j],
				mopriority[i], moabstol[j], moreltol[j], obj_name(j))))
			badret("CPXmultiobjsetobj", k, 1);
		}
	return CPXmultiobjopt(env, cpx, (const struct paramset * const*)mopars);
	}

 static fint
MOkwf(char *key, fint klen)
{
	ASL *asl;
	MOkwval *v;
	char *key0;
	int c, on;
	keyword *kw;

	key0 = key;
	if (klen < 7 || strncmp(key, "obj_", 4)) {
 badkey:
		printf("Unrecognized keyword \"%s\"\n", key0);
		return 1;
		}
	on = 0;
	for(key += 4;;) {
		c = *key++;
		if (c == '_')
			break;
		if (c < '0' || c > '9')
			return 1;
		on = 10*on + (c - '0');
		}
	asl = Oinfo.asl;
	if (on <= 0 || on > n_obj) {
		printf("Rejecting obj_%d; obj_n must have 1 <= n <= %d\n", on, n_obj);
		return 1;
		}
	if (!(kw = (keyword*)b_search_ASL(Oinfo.keywds, (int)sizeof(keyword),
					Oinfo.n_keywds, &key, &Oinfo.eqsign)))
		goto badkey;
	*pnext_MOkw = v = (MOkwval*)mem(sizeof(MOkwval) + strlen(key) + 1);
	pnext_MOkw = &v->next;
	v->next = 0;
	v->kw = kw;
	strcpy(v->val = (char*)(v+1), key);
	v->Objno = on - 1;
	return 0;
	}
 static char modisplay_desc[] = "how much to report during multiobjective optimization:\n\
			0 = nothing\n\
			1 = summary after each subproblem (default)\n\
			2 = subproblem logs as well as summaries.";
 static char multiobj_desc[] = "whether to do multi-objective optimization:\n\
			0 = no (default)\n\
			1 = yes\n\
		When multiobj = 1 and several linear objectives are\n\
		present, suffixes .objpriority, .objweight, .objreltol,\n\
		and .objabstol on the objectives are relevant.\n\
		Objectives with greater (integer) .objpriority values\n\
		have higher priority.  Objectives with the same\n\
		.objpriority are weighted by .objweight.  Objectives\n\
		with positive .objabstol or .objreltol are allowed to\n\
		be degraded by lower priority objectives by amounts not\n\
		exceeding the .objabstol (absolute) and .objreltol\n\
		(relative) limits.  The objective must all be linear.\n\
		Objective-specific values may be assigned via keywords\n\
		of the form obj_n_name, such as obj_1_pricing to\n\
		specify \"pricing\" for the first objective.  If no\n\
		.objweight values are provided, 1. is assumed for all.\n\
		Similarly, if no .objpriority values are given, 1 is\n\
		assumed for all.  For .objreltol and .objabstol, if\n\
		no values are given, all are assumed to be 0."
		;

 static int
pricomp(const void *a, const void *b, void *c)
{
	int d, *x;

	x = (int*)c;
	d = x[*(int*)a] - x[*(int*)b];
	if (d < 0)
		return -1;
	if (d > 0)
		return 1.;
	return 0;
	}

 static void
MOcheck(ASL *asl)
{
	Kwfunc *kf;
	MOkwval *kwv;
	SufDesc *oat, *opr, *ort, *owt;
	char *s, *x;
	int Nlo, *apri, f, i, j, k, nbad, no, nzw, *operm, *opinv, *pri, z[3];
	keyword *kw;
	real *at, d, *rt, *wt, zd[3];
	size_t L;
	struct paramset *ps;
	static char BadOption[] = "bad multiobjective keyword assignment";
	static char NonlinOjb[] = "multiobjective problem involving a nonlinear objective";
	static char Toofew[] = "fewer than 2 nonzero multiobjective weights";

	pri = 0;
	at = rt = wt = 0;
	no = n_obj;
	L = no*(4*sizeof(int) + sizeof(struct parmset*));
	if ((oat = suf_get("objabstol", ASL_Sufkind_obj)))
		at = oat->u.r;
	if ((opr = suf_get("objpriority", ASL_Sufkind_obj)))
		pri = opr->u.i;
	if ((ort = suf_get("objreltol", ASL_Sufkind_obj)))
		rt = ort->u.r;
	else
		L +=  no*sizeof(real);
	if ((owt = suf_get("objweight", ASL_Sufkind_obj)))
		wt = owt->u.r;
	else
		L += no*sizeof(real);
	if (!at)
		L +=  no*sizeof(real);
	if (!rt)
		L +=  no*sizeof(real);
	if (!wt)
		L +=  no*sizeof(real);
	momem = (real*)M1alloc(L);
	x = (char*)momem;
	if (!at) {
		memset(at = (real*)x, 0, no*sizeof(real));
		x = (char*)(at + no);
		}
	if (!(rt)) {
		memset(rt = (real*)x, 0, no*sizeof(real));
		x = (char*)(rt + no);
		}
	if (!(wt)) {
		wt = (real*)x;
		for(i = 0; i < no; ++i)
			wt[i] = 1.;
		x = (char*)(wt + no);
		}
	memset(mopars = (struct paramset **)x, 0, no*sizeof(struct paramset*));
	x = (char*)(mopars + no);
	operm = (int*)x;
	opinv = operm + no;
	pri0 = opinv + no;
	apri = pri0 + no;
	x = (char*)(apri + no);
	for(i = 0; i < no; ++i)
		operm[i] = i;
	nzw = 0;
	if (owt) {
		for(i = 0; i < no; ++i) {
			if (wt[i] == 0.) {
				operm[i] = -1;
				++nzw;
				}
			}
		}
	if ((Nlo = nlo)) {
		for(i = 0; i < Nlo; ++i) {
			if (operm[i] >= 0) {
				solve_result_num = 532;
				asl->i.uinfo = NonlinOjb;
				longjmp(Jb, 1);
				}
			}
		}
	for(i = 0; i < no; ++i)
		opinv[i] = -1;
	if (nzw) {
		if (nzw >= no - 1) {
			solve_result_num = 534;
			asl->i.uinfo = Toofew;
			longjmp(Jb, 1);
			}
		for(i = j = 0; i < no; ++i) {
			if (operm[i] >= 0)
				operm[j++] = i;
			}
		no = j;
		}
	nmo = no;
	if (pri) {
		qsortv(operm, no, sizeof(int), pricomp, pri);
		apri[0] = j = 1;
		pri0[0] = pri[operm[0]];
		for(i = 1; i < no; ++i) {
			if ((k = pri[operm[i]]) > pri[operm[i-1]])
				pri0[j++] = k;
			apri[i] = j;
			}
		nmopri = j;
		}
	else {
		nmopri = 1;
		pri0[0] = 0;
		for(i = 0; i < no; ++i)
			apri[i] = 1;
		}
	for(i = 0; i < no; ++i)
		opinv[operm[i]] = i;
	nbad = 0;
	for(kwv = first_MOkw; kwv; kwv = kwv->next) {
		if ((j = operm[kwv->Objno]) >= 0) {
			j = apri[j];
			if (!(ps = mopars[j])) {
				if (!(mopars[j] = ps = CPXparamsetcreate(Env, &i)))
					badret("CPXparamsetcreate", i, 1);
				}
			kw = kwv->kw;
			kf = kw->kf;
			if (kf == sf_int) {
				j = (int)strtol(s = kwv->val, &s, 10);
				if (s == kwv->val) {
					printf("Expected an integer value for %s, not \"%s\"\n",
						kw->name, s);
					badopt_ASL(&Oinfo);
					}
				else if ((i = CPXparamsetaddint(Env, ps, f = Intcast kw->info, j))) {
					++nbad;
					if (i == CPXERR_PARAM_INCOMPATIBLE)
						incompatible(&Oinfo, kw, s, kwv->val);
					else {
						z[2] = 0;
						z[1] = 1;
						CPXinfointparam(Env, f, z, z+1, z+2);
						badival(&Oinfo, kw, j, z[1], z[2]);
						}
					}
				}
			else if (kf == sf_dbl) {
				d = strtod(s = kwv->val, &s);
				if (s == kwv->val) {
					printf("Expected a numeric value for %s, not \"%s\"\n",
						kw->name, s);
					badopt_ASL(&Oinfo);
					}
				else if ((i = CPXparamsetadddbl(Env, ps, f = Intcast kw->info, d))) {
					++nbad;
					if (i == CPXERR_PARAM_INCOMPATIBLE)
						incompatible(&Oinfo, kw, s, kwv->val);
					else {
						zd[2] = 0.;
						zd[1] = 1.;
						CPXinfodblparam(Env, f, zd, zd+1, zd+2);
						printf("Rejecting obj_%d_%s %g; must be "
							"between %g and %g\n", kwv->Objno+1,
							kw->name, d, z[1], z[2]);
						badopt_ASL(&Oinfo);
						}
					}
				}
			else {
				++nbad;
				printf("Rejecting inappropriate obj_%d_%s assignment.\n",
					kwv->Objno + 1, kw->name);
				badopt_ASL(&Oinfo);
				}
			}
		}
	if (nbad) {
		solve_result_num = 533;
		asl->i.uinfo = BadOption;
		longjmp(Jb, 1);
		}
	mopriority = apri;
	moabstol = at;
	objind = operm;
	indobj = opinv;
	moreltol = rt;
	moweight = wt;
	}
#endif /*}*/

#if 0
 char undoc[] = "Undocumented.";
#endif
#ifdef CPLEX_MIP /*{*/
 static char
	absmipap_desc[] = "Absolute mixed-integer optimality gap tolerance\n\
		(for difference between current best integer solution\n\
		and optimal value of LP relaxation).  Default 0.",
	aggcutlim_desc[] = "Bound on the number of constraints aggregated to\n\
		generate flow-cover and mixed-integer-rounding cuts;\n\
		default = 3.",
#ifdef CPX_PARAM_AUXROOTTHREADS
	auxrootthreads_desc[] = "Controls the number of threads used for auxiliary\n\
		chores when solving the root node of a MIP problem.\n\
		When N threads are available (possibly limited by\n\
		\"threads\"), auxrootthreads must be less than N.\n\
		Possible values:\n\
		    0 = automatic choice (default)\n\
		    n < N:  use N-n threads for the root node and\n\
		         n threads for auxiliary chores.",
#endif
	bbinterval_desc[] = "For nodeselect = 2, select the best-bound node,\n\
		rather than the best-estimate node, every\n\
		bbinterval iterations (default 7); 0 means\n\
		always use the best-estimate node.",
#ifdef CPXPARAM_Benders_Strategy /*{*/
	benders_feascut_desc[] = "Tolerance for violations of feasibility cuts in Benders\n\
		algorithm.  Default = 1e-6.",
	benders_optcut_desc[] = "Tolerance for violations of optimality cuts in Benders\n\
		algorithm.  Default = 1e-6.",
	benders_strategy_desc[] = "How to decompose the problem for Benders algorithm:\n\
		   -1 = do not apply Benders algorithm\n\
		    0 = automatic choice (default): if suffix benders\n\
				is present on variables, variables that have\n\
				.benders = 0 go into the master and CPLEX\n\
				assigns other variables to workers; otherwise\n\
				integer variables go into the master and\n\
				continuous variables into workers\n\
		    1 = use suffix benders to determine which variables\n\
				are for the master (.benders = 0) and which\n\
				for workers (.benders = n > 0 ==> worker n\n\
		    2 = similar to 0, but suffix benders is required\n\
		    3 = similar to 0, but ignore suffix benders (if\n\
				present).",
	bendersopt_desc[] = "Single-word phrase:  use Benders algorithm.\n\
		Both integer and continuous variables must be present.",
#endif /*}*/
	bestbound_desc[] = "Single-word phrase requesting return of suffix\n\
		.bestbound on the objective and problem for the\n\
		best known bound on the objective value.  For MIP\n\
		problems with .bestnode value from a feasible node\n\
		(see below), .bestbound = .bestnode.",
	bestnode_desc[] = "Single-word phrase requesting return of suffix\n\
		.bestnode on the objective and problem for the\n\
		objective value at the best feasible MIP node.\n\
		For non-MIP problems and for MIP problems for which\n\
		a feasible node has not yet been found, this value\n\
		is +Infinity for minimization problems and -Infinity\n\
		for maximization problems.",
	boundstr_desc[] = "Whether to use bound strengthening in solving MIPs:\n\
		   -1 (default) = automatic choice\n\
		    0 = never\n\
		    1 = always.",
#ifdef CPX_PARAM_BQPCUTS
	bqpcuts_desc[] = "Whether to generate boolean quadratic polytope (BQP)\n\
		cuts for nonconvex QP amd MIQP problems when solved\n\
		to optimality:\n\
		   -1 = do not generate BQP cuts\n\
		    0 = automatic choice (default)\n\
		    1 = generate BQP cuts moderateely\n\
		    2 = generate BQP cuts agressively\n\
		    3 = generate BQP cuts very agressively.",
#endif
	branch_desc[] = "Branching direction for integer variables:\n\
		-1 = down, 0 = algorithm decides, 1 = up; default = 0.",
	branchdir_desc[] = "Synonym for \"branch\".",
	cliquecuts_desc[] = "Synonym for \"cliques\".",
	cliques_desc[] = "Whether to use clique cuts in solving MIPs:\n\
		   -1 = never\n\
		   0 = automatic choice (default)\n\
		   1, 2, 3 = ever more aggressive generation.",
	coeffreduce_desc[] = "Whether to use coefficient reduction when\n\
		preprocessing MIPS:\n\
		  -1 = automatic choice (default)\n\
		   0 = no\n\
		   1 = reduce only integral coefficients\n\
		   2 = reduce all potential coefficients\n\
		   3 = reduce aggressively with tiling.",
	covercuts_desc[] = "Synonym for \"covers\".",
	covers_desc[] = "Whether to use cover cuts in solving MIPs:\n\
		   -1 = never\n\
		   0 = automatic choice (default)\n\
		   1, 2, 3 = ever more aggressive generation.",
#ifdef CPX_PARAM_CPUMASK
	cpumask_desc[] = "Whether and how to bind threads to cores\n\
		on systems where this is possible:\n\
		    off = no CPU binding\n\
		    auto = automatic binding (default).\n\
		Values other than \"off\" and \"auto\" must be a\n\
		hexadecimal string (digits 0-9 and a-f, ignoring\n\
		case, so values A-F and a-f are treated alike).\n\
		The lowest order bit is for the first logical CPU.\n\
		For example, \"a5\" and \"A5\" indicate that CPUs 0, 2,\n\
		5, and 7 are available for binding to threads, since\n\
		hex value a5 = 2^7 + 2^5 + 2^2 + 2^0.",
#endif
	cutpass_desc[] = "Number of passes permitted when generating MIP\n\
		cutting plane:\n\
		   -1 = none\n\
		   0 = automatic choice (default)\n\
		   positive = at most that many passes.",
cutsfactor_desc[] =
#if CPX_VERSION >= 12060200
		"Limit on MIP cuts added:\n\
		   > 1 ==> (cutsfactor-1)*m, where m\n\
		is the original number of rows (after presolve);\n\
		   < 0 ==> no limit;\n\
		   0 <= cutsfactor <= 1 ==> no MIP cuts\n\
		Default = -1 (no limit).",
#else
		"Limit MIP cuts added to (cutsfactor-1)*m, where m\n\
		is the original number of rows (after presolve).\n\
		Default = 4.",
#endif
#if CPX_VERSION >= 1100
	cutstats_desc[] = "0 or 1 (default 0):  Whether the solve_message report\n\
		the numbers and kinds of cuts used.",
#endif
	disjcuts_desc[] = "Whether to generate MIP disjunctive cuts:\n\
		   -1 = no\n\
		   0 = automatic choice (default)\n\
		   1, 2, 3 = ever more aggressive generation.",
#ifdef CPX_PARAM_IISIND /* version < 9.2b */
	endtree_desc[] = "File for writing final branch & bound search tree.\n\
		Withdrawn from CPLEX 10.",
#endif
	flowcuts_desc[] = "Whether to use flow cuts in solving MIPs:\n\
		   -1 = never\n\
		   0 = automatic choice (default)\n\
		   1, 2 = ever more aggressive use.",
	flowpathcuts_desc[] = "Whether to generate MIP flow-path cuts:\n\
		   -1 = no\n\
		   0 = automatic choice (default)\n\
		   1, 2 = ever more aggressive generation.",
#ifdef CPX_PARAM_FPHEUR
	fpheur_desc[] = "Whether to use the feasibility pump heuristic on MIP\n\
		problems:\n\
		   -1 = no\n\
		    0 = automatic choice (default)\n\
		    1 = yes, focus on finding a feasible solution\n\
		    2 = yes, focus on finding a good objective\n\
		value at a feasible solution.",
#endif
#ifndef NO_CPLEX66 /* for versions prior to CPLEX 6.6 */
	fraccand_desc[] = "Limit on number of candidate variables when\n\
		generating Gomory cuts for MIP problems:\n\
		default = 200.",
	fraccuts_desc[] = "Whether to generate MIP fractional Gomory\n\
		cuts:\n\
		   -1 = no\n\
		    0 = decide automatically (default)\n\
		    1 = generate moderately\n\
		    2 = generate aggressively.",
	fracpass_desc[] = "Limit on number of passes to generate MIP\n\
		fractional Gomory cuts:\n\
		   0 = automatic choice (default)\n\
		   positive = at most that many passes.",
	fractionalcuts_desc[] = "Synonym for \"fracpass\".",
#endif
	gubcuts_desc[] = "Whether to use GUB cuts in solving MIPs:\n\
		   -1 = never\n\
		   0 = automatic choice (default)\n\
		   1, 2 = ever more aggressive generation.",
#ifdef CPX_PARAM_HEURFREQ
	heurfreq_desc[] = "How often to apply \"node heuristics\" for MIPS:\n\
		   -1 = never\n\
		   0 = automatic choice (default)\n\
		   n > 0 = every n nodes.",
#endif
#ifdef CPX_PARAM_HEURISTIC
	heuristic_desc[] = "Deprecated:  withdrawn from CPLEX 9.0.\n\
		What heuristic to use on MIPS for finding an\n\
		initial integer solution:\n\
		   -1 = none\n\
		   0 = decide automatically (default)\n\
		   1 = use a rounding heuristic at node 0.",
#endif
#ifdef CPX_PARAM_HEURFREQ
	heuristicfreq_desc[] = "Synonym for \"heurfreq\".",
#endif
	iisfind_desc[] = "Whether to find and return an IIS (irreducible\n\
		infeasible set of variables and constraints) if\n\
		the problem is infeasible:\n\
		0 = no (default)\n\
		1 = find an IIS.\n\
		IIS details are returned in suffix .iis, which\n\
		assumes one of the values \"non\" for variables\n\
		and constraints not in the IIS; \"low\" for\n\
		variables or inequality constraint bodies whose lower\n\
		bounds are in the IIS; \"upp\" for variables and\n\
		inequality constraint bodies whose upper bounds are\n\
		in the IIS; and \"fix\" for equality constraints that\n\
		are in the IIS.",
	impliedcuts_desc[] = "Whether to use implied cuts in solving MIPs:\n\
		   -1 = never\n\
		    0 = automatic choice (default)\n\
		    1, 2 = ever more aggressive use.",
	integrality_desc[] = "Amount by which an integer variable can differ\n\
		from the nearest integer and still be considered\n\
		feasible.  Default = 1e-5; must be in [1e-9, 0.5].\n\
		(The upper bound was not enforced prior to CPLEX 11.)",
	intwarntol_desc[] = "Do not warn about perturbations to \"integer\"\n\
		variables to make them integers when the maximum\n\
		perturbation is at most intwarntol (default 1e-9);\n\
		see \"round\".",
#ifdef CPX_PARAM_LBHEUR
	lbheur_desc[] = "Whether to use a local branching heuristic in an\n\
		attempt to improve new incumbents found during a\n\
		MIP search.  (Default = 0 = no; 1 = yes.)",
#endif
#ifdef CPXPARAM_MIP_Cuts_LocalImplied
	localimpliedcuts_desc[] = "Whether to generate locally valid implied bound\n\
		cuts for MIP problems:\n\
		   -1 ==> no\n\
		    0 ==> automatic choice (default)\n\
		    1 ==> yes, moderately\n\
		    2 ==> yes, aggressively\n\
		    3 ==> yes, very aggressively.",
#endif
	lowercutoff_desc[] = "For maximization problems involving integer\n\
		variables, skip any branch whose LP relaxation's\n\
		optimal value is less than lowercutoff.  Warning:\n\
		if lowercutoff is too large, the problem will\n\
		appear infeasible.  Default = -1e75.",
#ifdef CPX_PARAM_MCFCUTS
	mcfcuts_desc[] = "Whether to use multi-commodity flow (MCF) cuts:\n\
		   -1 = no\n\
		    0 = let CPLEX decide (default)\n\
		    1 = generate a modest number of MCS cuts\n\
		    2 = generate MCS cuts aggressively.",
#endif
#if 0 /* old mipalg_desc; not sure when it change to the one for 12.6.1, given below */
	mipalg_desc[] = "Algorithm used on mixed-integer subproblems\n\
		(default 2):\n\
		   1 = primal simplex\n\
		   2 = dual simplex (plus primal, if dual fails)\n\
		   3 = network alg., then dual simplex\n\
		   4 = barrier with crossover\n\
		The next two are not available in CPLEX versions >= 8;\n\
		specify mipcrossover instead.\n\
		   5 = dual to iteration limit, then barrier\n\
		   6 = barrier without crossover.",
#else /* valid at least for CPLEX 12.6.1 */
	mipalg_desc[] = "Algorithm used on mixed-integer subproblems:\n\
		   0 = automatic choice (default)\n\
		   1 = primal simplex\n\
		   2 = dual simplex\n\
		   3 = network simplex\n\
		   4 = barrier\n\
		   5 = sifting.\n\
		For MIQP problems (quadratic objective, linear\n\
		constraints), settings other than 3 and 5 are treated\n\
		as 0.  For MIQCP problems (quadratic objective and\n\
		constraints), all settings are treated as 4.",
#endif
	mipalgorithm_desc[] = "Synonym for \"mipalg\".",
	mipbasis_desc[] = "Whether to compute a basis and dual variables for MIP\n\
		problems when endbasis is not specified:\n\
		  -1 = default (described below)\n\
		   0 = no\n\
		   1 = yes\n\
		This keyword is new with driver version 20040716.\n\
		When endbasis is specified, mipbasis=1 is assumed.\n\
		Otherwise, when mipbasis=0 is specified for a MIP\n\
		problem, no solver-status values for variables are\n\
		returned to AMPL.  The default is to assume 1 unless\n\
		a quadratic objective or constraint is present, in\n\
		which case qcdual is assumed if quadratic constraints\n\
		are present and 0 is assumed otherwise (as finding a\n\
		basis can be time consuming).",
	mipcrossover_desc[] = "Crossover method used when using the barrier\n\
		method for MIP subproblems:\n\
		   -1 = no crossover\n\
		    0 (default) = automatic choice\n\
		    1 = primal\n\
		    2 = dual.",
	mipcuts_desc[] = "Sets all ten of cliques, covers, disjcuts,\n\
		flowcuts, flowpathcuts, fraccuts, gubcuts,\n\
		impliedcuts, mircuts and zerohalfcuts to the\n\
		specified value.",
	mipdisplay_desc[] = "Frequency of displaying branch-and-bound\n\
		information (for optimizing integer variables):\n\
		   0 (default) = never\n\
		   1 = each integer feasible solution\n\
		   2 = every \"mipinterval\" nodes\n\
		   3 = every \"mipinterval\" nodes plus\n\
		       information on LP relaxations\n\
		       (as controlled by \"display\")\n\
		   4 = same as 2, plus LP relaxation info\n\
		   5 = same as 2, plus LP subproblem info.",
	mipemphasis_desc[] = "Whether to emphasize seeking optimality\n\
		(0 = default) or finding feasible solutions (1).\n\
		For CPLEX versions >= 8, two other values are\n\
		possible:  emphasizing optimality over\n\
		feasibility (2) and emphasizing best bound (3).",
	mipgap_desc[] = "Relative tolerance for optimizing integer\n\
		variables: stop if\n\
		   abs((best bound) - (best integer))\n\
		       < mipgap * (1 + abs(best bound)).\n\
		Default = 1e-4; must be between 1e-9 and 1.",
	mipinterval_desc[] = "Frequency of node logging for mipdisplay >= 2.\n\
		Default = 0 ==> automatic choice.  Values n > 0 ==>\n\
		every n nodes and every new incumbent; n < 0 ==> less\n\
		frequently the more negative n is.",
#ifdef CPX_PARAM_MIPKAPPASTATS
	mipkappa_desc[] = "For MIP problems, whether to compute the \"MIP kappa\",\n\
		which summarizes the condition numbers of the optimal\n\
		bases seen while solving the problem:\n\
		 -1 = no\n\
		  0 = automatic choice (default)\n\
		  1 = compute for a sample of subproblems\n\
		  2 = compute for all subproblems (possibly expensive).",
#endif
	mipordertype_desc[] = "Synonym for \"ordertype\".",
#ifdef CPX_PARAM_MIPSEARCH
	mipsearch_desc[] = "Search strategy for mixed-integer problems, new\n\
		in CPLEX 11:\n\
		   0 = automatic choice (default)\n\
		   1 = traditional branch and cut\n\
		   2 = dynamic search.",
#endif
	mipsolutions_desc[] = "Stop branch-and-bound for integer variables\n\
		after finding \"mipsolutions\" feasible solutions.\n\
		Default = 2^31 - 1.",
	mipstart_desc[] = "Synonym for \"mipstartvalue\".",
#if 0 /* old mipstartalg_desc; not sure when it change to the one for 12.6.1, given below */
	mipstartalg_desc[] = "Which LP algorithm to use in solving the initial\n\
		MIP subproblem:\n\
		   1 = primal simplex\n\
		   2 = dual simplex (default)\n\
		   3 = netopt\n\
		   4 = barrier with crossover\n\
		The next two are not available in CPLEX versions >= 8;\n\
		specify mipcrossover instead.\n\
		   5 = dual simplex to iteration limit,\n\
		       then barrier\n\
		   6 = barrier without crossover.",
#else /* valid at least for CPLEX 12.6.1 */
	mipstartalg_desc[] = "For problems with integer variables, which algorithm\n\
		to use in solving the initial MIP subproblem:\n\
		   0 = automatic choice (default)\n\
		   1 = primal simplex\n\
		   2 = dual simplex\n\
		   3 = network simplex\n\
		   4 = barrier\n\
		   5 = sifting\n\
		   6 = concurrent (several at once, if possible).\n\
		For MIQP problems (quadratic objective, linear\n\
		constraints), setting 5 is treated as 0 and 6 as 4.\n\
		For MIQCP problems (quadratic objective & constraints),\n\
		all settings are treated as 4.",
#endif
	mipstartstatus_desc[] = "Whether to use incoming variable and constraint\n\
		statuses if the problem has integer variables:\n\
		   0 = no\n\
		   1 = yes (default).",
	mipstartvalue_desc[] = "Whether to use initial guesses in problems with\n\
		integer variables:\n\
		   0 = no\n\
		   1 = yes (default)"
#if CPX_VERSION < 12070000
			"\n\
		   2 = yes, using the deprecated CPXcopymipstart\n\
		       function (as formerly done), rather than\n\
		       the now preferred CPXaddmipstarts (starting\n\
		       with driver version 20110607); this is for\n\
		       debugging and is to be withdrawn.",
#else
			".",
#endif
	mipsubalg_desc[] = "Synonym for \"mipalg\".",
#ifdef CPX_PARAM_MIPTHREADS
	mipthreads_desc[] = "Maximum threads for the MIP algorithm\n\
		(1 unless you have a CPLEX license for multiple\n\
		threads).  Withdrawn from CPLEX 11.",
#endif
#ifdef CPX_PARAM_MIQCPSTRAT
	miqcpstrat_desc[] = "Strategy for solving quadratically-constrained MIPs\n\
		(MIQCP problems):\n\
		   0 = automatic choice (default)\n\
		   1 = solve a quadratically-constrained node\n\
		       relaxation (QCP) at each node\n\
		   2 = solve an LP node relaxation at each node.",
#endif
	mircuts_desc[] = "Whether to generate MIP rounding cuts:\n\
		   -1 = no\n\
		    0 = automatic choice (default)\n\
		    1 = moderate generation\n\
		    2 = aggressive generation.",
	node_desc[] = "Synonym for \"nodes\".",
	nodefile_desc[] = "Whether to save node information in a temporary file:\n\
		   0 = no\n\
		   1 (default) = compressed node file in memory\n\
		   2 = node file on disk\n\
		   3 = compressed node file on disk.",
	nodefiledir_desc[] = "Synonym for workfiledir.  Prior to CPLEX 7.1,\n\
		this directory is just for node information files.",
#ifdef CPX_PARAM_NODEFILELIM
	nodefilelim_desc[] = "Maximum size for the node file.  Default 1e75.\n\
		Removed from CPLEX 7.1.  See workfilelim.",
	nodefilesize_desc[] = "Synonym for \"nodefilelim\".",
#endif
	nodelim_desc[] = "Synonym for \"nodes\".",
	nodes_desc[] = "Stop branch-and-bound for integer variables\n\
		after \"nodes\" LP relaxations.  Default = 2^31 - 1;\n\
		nodes = 0 ==> complete processing at the root (create\n\
		cuts, apply heuristics at root);\n\
		1 ==> allow branching from root:  create but do not\n\
		solve nodes.",
	nodesel_desc[] = "Strategy for choosing next node while optimizing\n\
		integer variables:\n\
		   0 = depth-first search;\n\
		   1 = breadth-first search (default);\n\
		   2 = best-estimate search;\n\
		   3 = alternate best-estimate search.",
	nodeselect_desc[] = "Synonym for \"nodesel\".",
	nosolve_desc[] = "Stop after loading the problem and honoring any\n\
		\"writeprob\" or \"writemipstart\" directives.",
	ordertype_desc[] = "How to generate default priorities for integer\n\
		variables when no .priority suffix is specified:\n\
		   0 = do not generate priorities (default)\n\
		   1 = use decreasing costs\n\
		   2 = use increasing bound range\n\
		   3 = use coefficient count.",
	plconpri_desc[] = "Priority (default 1) for SOS2 constraints for nonconvex\n\
		piecewise-linear terms in constraints.",
	plobjpri_desc[] = "Priority (default 2) for SOS2 constraints for nonconvex\n\
		piecewise-linear terms in objectives.",
#ifdef CPX_PARAM_POLISHAFTEREPAGAP
	polishafter_absmipgap_desc[] = "Start polishing integer solutions after the\n\
		absolute mixed-integer optimality gap is at most\n\
		polishafter_absmipgap.  Default 1e-6.",
	polishafter_intsol_desc[] = "Start polishing integer solutions after the\n\
		finding polishafter_intsol integer-feasible\n\
		solutions.  Default 2^31 - 1.",
	polishafter_mipgap_desc[] = "Start polishing integer solutions after the\n\
		relative mixed-integer optimality gap is at most\n\
		polishafter_mipgap.  Default 0.",
	polishafter_nodes_desc[] = "Start polishing integer solutions after the\n\
		processing polishafter_nodes nodes.\n\
		Default 2^31 - 1.",
	polishafter_time_desc[] = "Start polishing integer solutions after finding\n\
		at least one integer feasible solution and spending\n\
		polishafter_time CPU seconds seeking integer\n\
		solutions.  Default 1e75.",
#endif
#ifdef  CPX_PARAM_POLISHAFTERDETTIME
	polishafter_timedet_desc[] = "Start polishing integer solutions after finding\n\
		at least one integer feasible solution and spending\n\
		polishafter_time \"ticks\" seeking integer solutions.\n\
		Default 1e75.",
#endif
#ifdef CPX_PARAM_POLISHTIME
	polishtime_desc[] = "New in CPLEX 10.0:  seconds to spend \"polishing\"\n\
		integer solutions.  Default 0.\n\
		Deprecated in CPLEX 11.2 (in favor of the polishafter\n\
		keywords above).",
#endif
#ifdef CPX_PARAM_POPULATELIM
	poolagap_desc[] = "Solutions with objective worse in absolute value by\n\
		poolgap than the best solution are not kept in the\n\
		solution pool; see poolstub.  Default 1e75.",
	poolcapacity_desc[] = "Number of solutions to keep in solution pool;\n\
		see poolstub.  Default 2100000000.",
	pooldual_desc[] = "Whether to return dual values (corresponding to fixed\n\
		integer variables) for MIP and MIQP problems in the\n\
		solution pool:\n\
		   0 = no (default)\n\
		   1 = yes (which takes extra computation)",
	poolgap_desc[] = "Solutions  with objective worse in a relative sense by\n\
		poolgap than the best solution are not kept in the\n\
		solution pool; see poolstub.  Default 1e75.",
	poolintensity_desc[] = "How hard to try adding MIP solutions to the solution\n\
		pool.  Useful only if poolstub is specified.  Default 0\n\
		is treated as 1 if poolstub is specified without\n\
		populate, or 2 if populate is specified.  Larger values\n\
		(3 or 4) cause more additions to the solution pool,\n\
		possibly consuming considerable time; poolintensity 4\n\
		tries to generate all MIP solutions, which could be a\n\
		very large number.",
	poolreplace_desc[] = "Policy for replacing solutions in the solution pool if\n\
		more than poolcapacity solutions are generated:\n\
		   0 = FIFO (first-in, first-out); default\n\
		   1 = Keep best solutions\n\
		   2 = Keep most diverse solutions.",
	poolstub_desc[] = "Stub for solution files in the MIP solution pool.\n\
		New in CPLEX 11 and meaningful only if some variables\n\
		are integer or binary.  A pool of alternate MIP\n\
		solutions is computed if poolstub is specified, and the\n\
		solutions that remain in the solution pool (after some\n\
		are replaced if more than poolcapacity solutions are\n\
		found) are written to files\n\
		  (poolstub & '1') ... (poolstub & |solution pool|),\n\
		where |solution pool| is the number of solutions in the\n\
		solution pool.  That is, file names are obtained by\n\
		appending 1, 2, ... |solution pool| to poolstub.  The\n\
		value of |solution pool| is returned in suffix npool\n\
		on the objective and problem.",
	populate_desc[] = "Whether to run CPLEX's \"populate\" algorithm in an\n\
		attempt to add more solutions to the MIP solution pool.\n\
		   0 = no; just keep solutions found during the\n\
		       initial solve\n\
		   1 = run \"populate\" after finding a MIP solution\n\
		   2 = run \"populate\" instead of seeking a single\n\
		       best solution.\n\
		See poolstub.",
	populatelim_desc[] = "Limit on number of solutions added to the solution pool\n\
		by the populate algorithm.  See poolstub and populate.\n\
		Default 20.",
	prerelax_desc[] = "Whether to use CPLEX's presolve on the initial LP\n\
		relaxation of a MIP:\n\
		   -1 = automatic choice (default)\n\
		    0 = no\n\
		    1 = yes.",
#endif
	priorities_desc[] = "Whether to consider priorities for MIP branching:\n\
		0 = no\n\
		1 = yes (default).",
	probe_desc[] = "Whether to do variable probing when solving MIPs\n\
		(which sometimes dramatically affects performance,\n\
		for better or worse):\n\
		   -1 = no\n\
		    0 = automatic choice (default)\n\
		    1, 2, or 3 = ever more probing.",
#ifdef CPX_PARAM_PROBETIME
	probetime_desc[] = "Limit in seconds on time spent probing.\n\
		Default = 1e75.",
#endif
#ifdef CPX_PARAM_PROBEDETTIME
	probetimedet_desc[] = "Limit in \"ticks\" on time spent probing.\n\
		Default = 1e75.",
#endif
#ifdef Want_Distmipopt
	rampup_duration_desc[] = "How to ramp up distributed parallel optimization:\n\
		   -1 = no ramp up\n\
		    0 = automatic choice (default)\n\
		    1 = dynamic rampup, limited by rampup_timelim and\n\
		        rampup_walltimelim\n\
		    2 = infinite ramp up:  concurrent MIP optimization.",
	rampup_timelim_desc[] = "Time limit in deterministic \"ticks\" to spend on the\n\
		\"ramp up\" phase of distributed parallel optimization.\n\
		This only matters when rampup_duration = 0 or 1.\n\
		Default = 1e75.",
	rampup_walltimelim_desc[] = "Time limit limit in wall-clock seconds to spend on the\n\
		\"ramp up\" phase of distributed parallel optimization.\n\
		This only matters when rampup_duration = 0 or 1.\n\
		Default = 1e75.",
#endif
	relax_desc[] = "Single-word phrase:  ignore integrality.",
	relaxpresolve_desc[] = "Synonym for \"prerelax\".",
	relobjdif_desc[] = "Synonym for \"relobjdiff\".",
	relobjdiff_desc[] = "If the objdifference parameter is 0,\n\
		relobjdiff times the absolute value of the objective\n\
		value is added to (for maximizing) or subtracted\n\
		from (for minimizing) the best (so far) feasible\n\
		objective value while optimizing integer variables.\n\
		Subsequent nodes will be ignored if their LP\n\
		relaxations have optimal values worse than this\n\
		sum.  Default = 0.  Positive values may speed\n\
		the search -- and may cause the optimal solution\n\
		to be missed.",
	relpresolve_desc[] = "Synonym for \"prerelax\".",
#ifdef CPX_PARAM_REPAIRTRIES
	repairtries_desc[] = "How many times to try to repair in infeasible\n\
		MIP starting guess:\n\
		   -1 = none\n\
		    0 = automatic choice (default)\n\
		    > 0 = that many times.",
#endif
#ifdef CPX_PARAM_REPEATPRESOLVE
	repeatpresolve_desc[] = "Whether to repeat CPLEX's presolve at MIP nodes:\n\
		   -1 = automatic choice (default)\n\
		    0 = no\n\
		    1 = presolve again ignoring cuts\n\
		    2 = presolve again considering cuts\n\
		    3 = presolve again considering cuts and\n\
		        allowing new root cuts.",
#endif
#undef REQCONVEX
#ifdef CPX_PARAM_OPTIMALITYTARGET
#define REQCONVEX CPX_PARAM_OPTIMALITYTARGET
#elif defined(CPX_PARAM_SOLUTIONTARGET)
#define REQCONVEX CPX_PARAM_SOLUTIONTARGET
#endif
#ifdef REQCONVEX
	reqconvex_desc[] = "Whether to require a quadratic model to be convex:\n\
		   0 = automatic choice (default)\n\
		   1 = require convexity\n\
		   2 = do not require convexity; just look\n\
		       for a local solution"
#if CPX_VERSION >= 12060000
			"\n\
		   3 = globally solve if noncovex.",
#else
			".",
#endif
#endif
	resolve_desc[] = "Whether to re-solve the problem with CPLEX's\n\
		presolve turned off when it reports the problem\n\
		to be infeasible or unbounded.  Re-solving may\n\
		take extra time but should determine whether the\n\
		problem is infeasible or unbounded.\n\
		   0 = no\n\
		   1 = yes (default).",
	return_mipgap_desc[] = "Whether to return mipgap suffixes or include\n\
		mipgap values in the solve_message:  sum of\n\n\
		   1 = return relmipgap suffix\n\
		   2 = return absmipgap suffix\n\
		   4 = suppress mipgap values in solve_message\n\n\
		The suffixes are on the objective and problem;\n\
		returned suffix values are +Infinity if no integer-\n\
		feasible solution has been found, in which case no\n\
		mipgap values are reported in the solve_message.\n\
		Default = 0.  See also bestbound and bestnode above.",
#ifdef CPX_PARAM_RINSHEUR
	rinsheur_desc[] = "Relaxation INduced neighborhood Search HEURistic\n\
		for MIP problems:\n\
		   -1 = none\n\
		    0 = automatic choice of interval (default)\n\
		    n (for n > 0) = every n nodes.",
#endif
#ifdef CPXPARAM_MIP_Cuts_RLT
	rltcuts_desc[] = "Whether to use RLT (Reformulation Linearization\n\
		Technique) cuts:\n\
		   -1 = no\n\
		    0 = automatic choice (default)\n\
		    1 = generate RLT cuts moderately\n\
		    2 = generate RLT cuts aggressively\n\
		    3 = generate RLT cuts very aggressively",
#endif
#ifdef CPX_PARAM_HEURISTIC
	rootheuristic_desc[] = "Withdrawn from CPLEX 9.0.  Synonym for \"heuristic\".",
#endif
	round_desc[] = "Whether to round integer variables to integral\n\
		values before returning the solution, and whether\n\
		to report that CPLEX returned noninteger values\n\
		for integer values (default 1):  sum of\n\
		   1 ==> round nonintegral integer variables\n\
		   2 ==> do not modify solve_result\n\
		   4 ==> do not modify solve_message\n\
		   8 ==> modify solve_result and solve_message\n\
		even if maxerr < intwarntol (default 1e-9).\n\
		Modifications take place only if CPLEX assigned\n\
		nonintegral values to one or more integer variables.",
	solutionlim_desc[] = "Synonym for \"mipsolutions\".",
	sos_desc[] = "0 or 1 (default 1):  Whether to honor declared\n\
		suffixes .sosno and .ref describing SOS sets.\n\
		Each distinct nonzero .sosno value designates an SOS\n\
		set, of type 1 for positive .sosno values and of type\n\
		2 for negative values.  The .ref suffix contains\n\
		corresponding reference values.",
	sos2_desc[] = "0 or 1 (default 1): Whether to tell CPLEX about SOS2\n\
		constraints for nonconvex piecewise-linear terms.",
#ifdef CPX_PARAM_LANDPCUTS
	splitcuts_desc[] = "Whether to use lift-and-project cuts on MIP problems\n\
		(new for CPLEX 12.5.1):\n\
		   -1 = no\n\
		    0 = automatic choice (default)\n\
		    1 = moderate use\n\
		    2 = aggressive use\n\
		    3 = very aggressive use.",
#endif
	startalg_desc[] = "Synonym for \"mipstartalg\".",
	startalgorithm_desc[] = "Synonym for \"mipstartalg\".",
	startbasis_desc[] = "\"startbasis foo\" reads the initial basis\n\
		(in BAS format) from file \"foo\".",
	startsol_desc[] = "Synonym for \"readsol\".",
#ifdef CPX_PARAM_IISIND /* version < 9.2b */
	starttree_desc[] = "File for reading initial branch-and-bound search tree.",
#endif
	strongcand_desc[] = "Length of the candidate list for \"strong branching\"\n\
		when solving MIPs: default 10.",
	strongit_desc[] = "Number of simplex iterations on each variable in\n\
		the candidate list during strong branching.\n\
		Default = 0 = automatic choice.",
#ifdef CPX_PARAM_STRONGTHREADLIM
	strongthreads_desc[] = "Maximum threads during strong branching.\n\
		Default = 1.",
#endif
	subalg_desc[] = "Synonym for \"mipalg\".",
	subalgorithm_desc[] = "Synonym for \"mipalg\".",
#ifdef CPX_PARAM_SUBMIPSUBALG
	submipalg_desc[] = "Rarely used choice of algorithm for the initial\n\
		relaxation of a subMIP: not a subproblem, but an\n\
		auxiliary MIP that CPLEX sometimes forms and solves,\n\
		e.g., when\n\
			* dealing with a partial MIP start\n\
			* repairing an infeasible MIP start\n\
			* using the RINS heuristic\n\
			* branching locally\n\
			* polishing a solution.\n\
		Possible values (when appropriate):\n\
			0 = automatic choice (default)\n\
			1 = primal simplex\n\
			2 = dual simplex\n\
			3 = network simplex (not for MIQPs)\n\
			4 = barrier\n\
			5 = sifting (0 is used for MIQPs).\n\
		Only 0 is allowed for MIQCPs.",
#endif
#ifdef CPXPARAM_MIP_SubMIP_NodeLimit
	submipnodelim_desc[] = "Limit on nodes searched by relaxation induced\n\
		neighborhood search (RINS) heuristic for MIP\n\
		problems and for processing of MIP starting values.\n\
		Default = 500.",
#endif
#ifdef CPX_PARAM_SUBMIPSCAIND
	submipscale_desc[] = "Rarely used choice of scaling for auxiliary subMIPs\n\
		(described with \"submipalg\"):\n\
			-1 = no scaling\n\
			 0 = equilibration scaling (default)\n\
			 1 = more aggressive scaling.",
#endif
#ifdef CPX_PARAM_SUBMIPSTARTALG
	submipstart_desc[] = "Rarely used choice of algorithm for the initial\n\
		relaxation of a subMIP (described with \"submipalg\"):\n\
			0 = automatic choice (default)\n\
			1 = primal simplex\n\
			2 = dual simplex\n\
			3 = network simplex\n\
			4 = barrier\n\
			5 = sifting (0 is used for MIQPs)\n\
			6 = concurrent (dual, barrier and primal in\n\
				opportunistic mode; dual and barrier in\n\
				deterministic mode; 4 is used for MIPQs).\n\
		Only 0 is allowed for MIQCPs.",
#endif
#ifdef CPX_PARAM_SYMMETRY
	symmetry_desc[] = "Whether to break symmetry during\n\
		    preprocessing of MIP problems:\n\
		   -1 = automatic choice (default)\n\
		    0 = no\n\
		    1 = moderate effort\n\
		    2 = more effort\n\
		    3 = still more effort"
#if CPX_VERSION >= 11000000
			"\n\
		    4 = even more effort (new in CPLEX 11)\n\
		    5 = more effort than 4 (new in CPLEX 11).",
#else
		".",
#endif
#endif
	treelimit_desc[] = "Synonym for \"treememory\".",
	treememlim_desc[] = "Synonym for \"treememory\".",
	treememory_desc[] = "Max. megabytes of memory (default 1e75) to use for\n\
		branch-and-bound tree.",
	uppercutoff_desc[] = "For minimization problems involving integer\n\
		variables, skip any branch whose LP relaxation's\n\
		optimal value is more than uppercutoff.  Warning:\n\
		if uppercutoff is too small, the problem will\n\
		appear infeasible.  Default = 1e75.",
	varsel_desc[] = "Strategy for selecting the next branching\n\
		variable during integer branch-and-bound:\n\
		   -1 = branch on variable with smallest\n\
		        integer infeasibility\n\
		    0 = algorithm decides (default)\n\
		    1 = branch on variable with largest\n\
		        integer infeasibility"
#if CPX_VERSION >= 6050000
			"\n\
		    2 = branch based on pseudo costs\n\
		    3 = strong branching\n\
		    4 = branch based on pseudo reduced costs.",
#else
			"\n\
		    2 = branch based on pseudo reduced costs.",
#endif
	varselect_desc[] = "Synonym for \"varsel\"."
#ifdef Want_Distmipopt
	,vmconf_desc[] = "For distributed parallel MIP optimization, if vmconf\n\
		starts with @, then the remainder is the name of a file\n\
		containing a parallel MIP configuration; otherwise\n\
		vmconf itself is a parallel MIP configuration string,\n\
		which must be quoted if it contains white space."
#endif
#ifdef CPX_PARAM_ZEROHALFCUTS
	,zerohalfcuts_desc[] = "Whether to generate zero-half cuts for MIP problems:\n\
		   -1 = no\n\
		    0 = automatic choice (default):  continue\n\
		        generating until new cuts are not helpful\n\
		    1 = generate zero-half cuts moderately\n\
		    2 = generate zero-half cuts aggressively."
#endif
	;
#endif /*}*/

 static char
	advance_desc[] = "Whether to use advance basis information (initial\n\
		primal and dual variable values and basis indicators).\n\
		Default 1 (yes).",
	aggfill_desc[] = "Synonym for \"agglim\".",
	agglim_desc[] = "Variables that appear in more than agglim rows\n\
		(default 10) will not be substituted away by the\n\
		\"aggregate\" algorithm.",
	aggregate_desc[] = "Whether to make substitutions to reduce the number of\n\
		rows:  0 ==> no; n > 0 ==> apply aggregator n times.\n\
		Default -1 ==> automatic choice.",
	aggtol_desc[] = "Pivot tolerance for aggregating.  It seldom needs\n\
		fiddling.  Default = .05; must be in [1e-10, .99].",
	autoopt_desc[] = "Single-word phrase:  use CPLEX's automatic choice of\n\
		optimizer (currently dualopt for LPs).",
	autopt_desc[] = "Synonym for \"autoopt\".",
	backtrack_desc[] = "Tolerance (> 0, default 0.9999) for when to backtrack\n\
		during branch & bound.  Low values tend to pure\n\
		best-bound search.  High values (> 1) tend to pure\n\
		depth-first search.  Values less than the default\n\
		are often good when subproblems are expensive.",
#ifdef BARRIER
	baralg_desc[] = "How to start the barrier algorithm:\n\
			   0 (default) = 1 for MIP subproblems, else 3\n\
			   1 = infeasibility-estimate start\n\
			   2 = infeasibility-constant start\n\
			   3 = standard start.",
	barcorr_desc[] = "Limit on centering corrections in each iteration\n\
		of the barrier algorithm:\n\
			   -1 = decide automatically (default)\n\
			   nonnegative = at most that many.",
	bardisplay_desc[] = "Specifies how much the barrier algorithm chatters:\n\
			   0 = no output (default)\n\
			   1 = one line per iteration\n\
			   2 = more output.",
	bargrowth_desc[] = "Tolerance for detecting unbounded faces in the\n\
		barrier algorithm: higher values make the test\n\
		for unbounded faces harder to satisfy.\n\
		Default = 1e12.",
	bariterlim_desc[] = "Maximum barrier iterations allowed (default 2^31 - 1).",
	barobjrange_desc[] = "Limit on the absolute objective value before the\n\
		barrier algorithm considers the problem unbounded.\n\
		Default = 1e20.",
	baropt_desc[] = "Single-word phrase:  use the barrier algorithm\n\
		(unless there are discrete variables).",
#ifdef CPX_PARAM_BAROOC
	baroutofcore_desc[] = "Whether the barrier solver should use disk\n\
		(out-of-core) storage for Cholesky factors:\n\
		    0 = no (default)\n\
		    1 = yes.\n\
		Withdrawn in CPLEX 10.0.",
#endif
	barstart_desc[] = "Barrier starting-point algorithm:\n\
			   1 = assume dual is 0 (default)\n\
			   2 = estimate dual\n\
			   3 = average of primal estimate, 0 dual\n\
			   4 = average of primal and dual estimates.",
	barstartalg_desc[] = "Synonym for \"barstart\".",
#ifdef CPX_PARAM_BARTHREADS
	barthreads_desc[] = "Maximum threads for the barrier algorithm\n\
		(1 unless you have a CPLEX license for multiple\n\
		threads).  Withdrawn from CPLEX 11.",
#endif
#ifdef CPX_PARAM_BARVARUP
	barvarup_desc[] = "Upper bound imposed by barrier algorithm on\n\
		variables with infinite upper bound; used to\n\
		prevent trouble with problems having unbounded\n\
		optimal faces.  Default = 1e20.\n\
		Withdrawn from CPLEX 9.0.",
#endif
#endif /* BARRIER */
#ifdef BASDEBUG
	basdebug_desc[] = "Filename for output to debug basis transmission.",
#endif
	basis_cond_desc[] = "Whether to show the condition number of the simplex\n\
		basis in the solve_message and to return its value\n\
		in the problem.basis_cond and objective.basis_cond\n\
		suffixes.  (Default = 0 = no; 1 = yes).",
#ifdef CPX_PARAM_BASINTERVAL
	basisinterval_desc[] = "Number of interations between savings of current\n\
		simplex basis to a file.  Default = 50000.\n\
		Deprecated in 10.0, removed in 12.6.1.",
#endif
	clocktype_desc[] = "Kind of times CPLEX reports during the solution\n\
		process:\n\
		   0 = automatic choice\n\
		   1 = CPU time\n\
		   2 = wall clock time (total elapsed time = default)",
#ifdef BARRIER
	comptol_desc[] = "Convergence tolerance for barrier algorithm:\n\
		the algorithm stops when the relative\n\
		complementarity is < bartol (default 1e-8).",
#endif
	concurrent_desc[] = "Single-word phrase:  with CPLEX versions >= 8\n\
		and when hardware and licensing permit, try\n\
		several methods in parallel.",
	concurrentopt_desc[] = "Synonym for \"concurrent\".",
#ifdef CPXPARAM_Conflict_Algorithm
	conflictalg_desc[] = "Choice of algorithm used by the CPLEX's conflict\n\
		refiner:\n\
		   0 = automatic choice (default)\n\
		   1 = fast\n\
		   2 = propagate\n\
		   3 = presolve\n\
		   4 = IIS\n\
		   5 = limited solve\n\
		   6 = full solve.\n\
		Settings 1, 2, and 3 are fast but may not discard\n\
		many constraints; 5 and 6 work harder at this.\n\
		Setting 4 searches for an Irreducible Infeasible\n\
		Set of linear constraints (e.g., ignoring quadratic\n\
		constraints).",
#endif
#ifdef CPX_PARAM_CONFLICTDISPLAY
	conflictdisplay_desc[] = "What to report when the conflict finder is working:\n\
		   0 = nothing\n\
		   1 = summary (default)\n\
		   2 = detailed display.",
#endif
	crash_desc[] = "Crash strategy (used to obtain starting basis);\n\
		possible values = -1, 0, 1; default = 1.\n\
		The best setting is problem-dependent and\n\
		can only be found by experimentation.\n\
		0 completely ignores the objective.",
	crossover_desc[] = "Causes the barrier algorithm to be run (in the\n\
		absence of discrete variables) and specifies\n\
		whether to \"crossover\" to an optimal simplex\n\
		basis afterwards:\n\
		   0 = no crossover\n\
		   1 = crossover with primal simplex\n\
		       (default for baropt)\n\
		   2 = crossover with dual simplex.",
#if CPX_VERSION >= 12070000
	datacheck_desc[] = "debug option; possible values:\n\
		   0 = no data checking (default)\n\
		   1 = issue warnings\n\
		   2 = try to \"assist\" by warning about bad scaling.",
#elif defined(WANT_DATACHECK)
	datacheck_desc[] = "debug option; gives inappropriate complaints\n\
		about Infinity"
		"\n\
		Possible values:\n\
		   0 = no data checking (default)\n\
		   1 = issue warnings.",
#endif
#ifdef BARRIER
	dense_desc[] = "Synonym for \"densecol\".",
	densecol_desc[] = "If positive, minimum nonzeros in a column for\n\
		the barrier algorithm to consider the column dense.\n\
		If 0 (default), this tolerance is selected\n\
		automatically.",
#endif
	dependency_desc[] = "Whether to use CPLEX's presolve dependency checker:\n\
		  -1 = automatic choice (default)\n\
		   0 = no\n\
		   1 = turn on only at start of preprocessing\n\
		   2 = turn on only at end of preprocessing\n\
		   3 = turn on at both start and end of\n\
		       preprocessing.",
#ifdef CPX_PARAM_DETTILIM
	dettimelim_desc[] = "Time limit in platform-dependent \"ticks\".\n\
		Default = 1e75.  See timing.",
#endif
	dgradient_desc[] = "Pricing algorithm for dual simplex (default 0):\n\
		   0 = choose automatically\n\
		   1 = standard dual pricing\n\
		   2 = steepest-edge pricing\n\
		   3 = steepest-edge pricing in slack space\n\
		   4 = steepest-edge with unit initial norms\n\
		   5 = devex pricing.",
	display_desc[] = "Frequency of displaying LP progress information:\n\
		   0 (default) = never\n\
		   1 = each factorization\n\
		   2 = each iteration.",
	doperturb_desc[] = "1 means initially perturb the problem (by an\n\
		amount governed by \"perturbation\", which is\n\
		described below).  0 (default) means let the\n\
		algorithm decide.  Setting doperturb to 1\n\
		is occasionally helpful for highly degenerate\n\
		problems.",
	dparam_desc[] = "Used with syntax \"dparam=n=d\" (no spaces), where n\n\
		is a decimal integer, the number of a CPLEX \"double\"\n\
		(i.e., floating-point valued) parameter.  If d is a\n\
		decimal floating-point value, assign d to \"double\"\n\
		parameter n.  If d is ?, report the current value of\n\
		\"double\" parameter n.  This facility provides a way\n\
		to modify \"double\" parameters that have not (yet)\n\
		been assigned a keyword.",
#if CPX_VERSION >= 1100
	droptol_desc[] = "If droptol > 0 is specified, linear constraint\n\
		and objective coefficients less than droptol in\n\
		magnitude are treated as zero.",
#endif
	dual_desc[] = "Single-word phrase:  solve the dual problem.",
	dualopt_desc[] = "Single-word phrase:  use a dual simplex algorithm.",
	dualratio_desc[] = "If neither \"primal\" nor \"dual\" was specified and\n\
		\"dual\" is possible (e.g., no integer variables and\n\
		no node and arc declarations), choose between\n\
		\"primal\" and \"dual\" as follows.\n\
		Let m = number of rows, n = number of columns;\n\
		if m - n > dualthresh > 0 or m > dualratio*n,\n\
		solve the dual; otherwise solve the primal.\n\
		Defaults:  dualthresh = 0, dualratio = 3.",
	dualthresh_desc[] = "See dualratio.",
#ifdef CPX_PARAM_EACHCUTLIM
	eachcutlim_desc[] = "Limit on the number of cuts of each time.\n\
		Default = 2100000000.",
#endif
	endbasis_desc[] = "\"endbasis foo\" writes the final basis to\n\
		file \"foo\" (in BAS format).",
	endsol_desc[] = "File for writing the final solution as an XML file.",
#if defined(BARRIER) && !defined(NO_DEPRECATED) && CPX_VERSION < 12070000
	endvector_desc[] = "File for writing solution from barrier algorithm\n\
		without crossover:  only meaningful with\n\
		\"baropt crossover=0\".  Deprecated; use endsol instead.",
#endif
	feasibility_desc[] = "Amount by which basic variables can violate\n\
		their bounds.  Default = 1e-6; possible\n\
		values are between 1e-9 and 0.1.",
#ifdef CPX_PARAM_FEASOPTMODE /* >= 9.2b */
	feasopt_desc[] = "For infeasible problems, whether to find a feasible\n\
		point for a relaxed problem (see feasoptobj):\n\
		   0 = no (default)\n\
		   1 = find a relaxed feasible point\n\
		   2 = find a \"best\" solution among the relaxed\n\
		       feasible points.",
	feasoptobj_desc[] = "Objective for \"feasopt\":\n\
		   1 = minimize the sum of constraint and variable\n\
		       bound relaxations (default)\n\
		   2 = minimize the number of constraint and variable\n\
		       bounds relaxed (a MIP problem, generally\n\
		       harder than feasoptobj = 1)\n\
		   3 = minimize the sum of squares of constraint and\n\
		       variable bound relaxations.",
#endif
	file_desc[] = "Synonym for \"writeprob\".",
#ifdef CPX_PARAM_FINALFACTOR
	finalfactor_desc[] = "Whether to factor the basis after \"uncrushing\"\n\
		the problem at the end.\n\
		   0 = no\n\
		   1 = yes (default)\n\
		Withdrawn in CPLEX 10.0.",
#endif
#ifdef BARRIER
	growth_desc[] = "Synonym for \"bargrowth\".",
#endif
#ifdef CPXERR_PARAM_INCOMPATIBLE
	incompat_desc[] = "How to treat parameter settings that CPLEX finds\n\
		incompatible:\n\
		   0 = quietly ignore incompatibilities\n\
		   1 = report and ignore them (default)\n\
		   2 = reject them, refusing to solve.\n\
		For example, CPLEX regards the polishafter_* parameters\n\
		introduced in CPLEX 11.2 as incompatible with the older\n\
		polishtime parameter.",
#endif
	iparam_desc[] = "Used with syntax \"iparam=n=i\" (no spaces), where n\n\
		is a decimal integer, the number of a CPLEX integer\n\
		parameter.  If i is a decimal integer, assign i to\n\
		integer parameter n.  If i is ?, report the current\n\
		value of integer parameter n.  This facility provides\n\
		a way to modify integer parameters that have not (yet)\n\
		been assigned a keyword.",
	iterations_desc[] = "Limit on total LP iterations; default 2^31 - 1.",
	iterlim_desc[] = "Synonym for \"iterations\".",
#ifdef Uselazy
	lazy_desc[] = "Whether to recognize suffix .lazy on constraints\n\
		(new for CPLEX 10): sum of\n\
		   1 ==> treat .lazy = 1 as lazy constraint\n\
		   2 ==> treat .lazy = 2 as user cut\n\
		Default lazy = 3 ==> treat both.  (Suffix .lazy on\n\
		constraints is ignored if not 0, 1, or 2 modulo 3.)",
#endif
	limitperturb_desc[] = "Synonym for \"perturblimit\".",
	logfile_desc[] = "Name of file to receive all CPLEX messages.",
	lowerobj_desc[] = "Stop minimizing when the objective value\n\
		goes below lowerobj.  Default = -1e75.",
	lowerobjlim_desc[] = "Synonym for \"lowerobj\".",
	lpdisplay_desc[] = "Synonym for \"display\".",
	lpiterlim_desc[] = "Synonym for \"iterations\".",
	lptimelim_desc[] = "Synonym for \"time\".",
	markowitz_desc[] = "Pivot tolerance; default = 0.01; must be between\n\
		0.0001 and 0.99999.  Bigger values may improve\n\
		numerical properties of the solution (and may\n\
		take more time).",
	maximize_desc[] = "Single-word phrase:  maximize the objective,\n\
		regardless of model specifications.",
#ifdef CPX_PARAM_MEMORYEMPHASIS
	memoryemphasis_desc[] = "Whether to compress data to reduce the memory used,\n\
		which may make some information (e.g., basis condition)\n\
		unavailable:\n\
		   0 = no (default)\n\
		   1 = yes.",
#endif
	minimize_desc[] = "Single-word phrase:  minimize the objective,\n\
		regardless of model specifications.",
	nameround_desc[] = "Whether to mangle variable and constraint names\n\
		by turning [ and ] into ( and ), respectively:\n\
		   0 = no (default)\n\
		   1 = yes.\n\
		This only matters if you specify endbasis=...\n\
		or startbasis=... or perhaps writeprob=something.lp\n\
		and have instructed AMPL to write .row and .col files.\n\
		(It is usually better to let AMPL's status facilities\n\
		convey basis information.)  An alternative under Unix\n\
		is to use the \"tr\" command to make the above changes\n\
		if they are needed.",
	netdisplay_desc[] = "Which objective value to show when using the\n\
		network simplex algorithm with display > 0\n\
		or netopt=3:\n\
		   0 = none\n\
		   1 = true objective\n\
		   2 = penalized objective (default).",
	netfeasibility_desc[] = "Feasibility tolerance for the network simplex\n\
		algorithm.  Default = 1e-6; possible values are\n\
		between 1e-11 and 1e-1.",
	netfind_desc[] = "Algorithm for finding embedded networks:\n\
		   1 = extract only the natural network\n\
		   2 = use reflection scaling (default)\n\
		   3 = use general scaling.",
	netfinder_desc[] = "Synonym for \"netfind\".",
	netiterations_desc[] = "Limit on network simplex iterations.\n\
		Default = large (e.g., 2^31 - 1).",
	netopt_desc[] = "0 means never invoke the network optimizer.\n\
		1 (default) means invoke the network optimizer\n\
		  only if the model had node and arc declarations.\n\
		2 means always invoke the network optimizer\n\
		  (unless there are integer variables); the network\n\
		  optimizer may be able to find and exploit an\n\
		  embedded network.\n\
		3 is similar to 2, but sets CPLEX's LPMethod\n\
		  to CPX_ALG_NET rather than explicitly invoking\n\
		  the network optimizer.  This might make a\n\
		  difference if CPLEX's presolve makes relevant\n\
		  reductions.",
	netoptimality_desc[] = "Tolerance for optimality of reduced costs in the\n\
		network simplex algorithm.  Default 1e-6; must be\n\
		between 1e-11 and 1e-1.",
	netpricing_desc[] = "How to price in the network simplex algorithm:\n\
		   0 = automatic choice (default)\n\
		   1 = partial pricing\n\
		   2 = multiple partial pricing\n\
		   3 = multiple partial pricing with sorting.",
#ifdef CPX_PARAM_NUMERICALEMPHASIS
	numericalemphasis_desc[] = "Whether to try to improve numerical accuracy (at a\n\
		possible cost of time or memory):\n\
		   0 = no (default)\n\
		   1 = yes.",
#endif
	objdifference_desc[] = "Amount added to (for maximizing) or subtracted\n\
		from (for minimizing) the best (so far) feasible\n\
		objective value while optimizing integer variables.\n\
		Subsequent nodes will be ignored if their LP\n\
		relaxations have optimal values worse than this\n\
		sum.  Default = 0.  Positive values may speed\n\
		the search -- and may cause the optimal solution\n\
		to be missed.",
	objno_desc[] = "1 (default) = first objective, 2 = second, etc.;\n\
		0 ==> no objective:  just find a feasible point.",
	objrep_desc[] = "Whether to replace\n\
		        minimize obj: v;\n\
		with\n\
		        minimize obj: f(x)\n\
		when variable v appears linearly in exactly one\n\
		constraint of the form\n\
		        s.t. c: v >= f(x);\n\
		or\n\
		        s.t. c: v == f(x);\n\
		Possible objrep values:\n\
		   0 = no\n\
		   1 = yes for v >= f(x)\n\
		   2 = yes for v == f(x) (default)\n\
		   3 = yes in both cases\n\
		For maximization problems, \">= f(x)\" is changed to\n\
		\"<= f(x)\" in the description above.  This is new\n\
		with driver version 20130622.",
	optimality_desc[] = "Tolerance for optimality of reduced costs.\n\
		Default 1e-6; must be between 1e-9 and 1e-1.",
	optimize_desc[] = "Synonym for \"primal\".",
#ifdef BARRIER
	ordering_desc[] = "Ordering algorithm used by the barrier algorithm\n\
		   0 = automatic choice (default)\n\
		   1 = approximate minimum degree\n\
		   2 = approximate minimum fill\n\
		   3 = nested dissection.",
#endif
	outlev_desc[] = "Synonym for \"display\".",
#ifdef CPX_PARAM_PARALLELMODE
	parallelmode_desc[] = "Parallel optimization mode:\n\
		   -1 = opportunistic mode\n\
		    0 = automatic: let CPLEX decide (default)\n\
		    1 = deterministic mode.",
#endif
	paramfile_desc[] = "File containing param settings to import.  The file\n\
		is read and settings in it echoed when the keyword\n\
		is processed.",
#if CPX_VERSION >= 1000
	paramfileprm_desc[] = "File containing param settings in CPLEX PRM format\n\
		to import.  The file is read without echoing settings\n\
		in it when the keyword is processed.",
#endif
#ifdef CPX_PARAM_PDSWITCH
	pdswitch_desc[] = "Whether to switch algorithms (from primal to dual\n\
		or vice versa) when undoing perturbations or shifts\n\
		made during the primal or dual simplex algorithm:\n\
		   -1 = no\n\
		    0 = automatic choice (default)\n\
		    1 = yes\n\
		pdswitch=1 may help if there are many cycles\n\
		of undoing perturbations.  New in CPLEX 7.0,\n\
		hidden in CPLEX 7.1.",
#endif
	perturb_desc[] = "Synonym for \"doperturb\".",
	perturbation_desc[] = "Amount by which to perturb variable bounds\n\
		when perturbing problems (see \"doperturb\").\n\
		Default 1e-6; must be positive.",
	perturbconst_desc[] = "Synonym for \"perturbation\".",
	perturblim_desc[] = "Number of stalled simplex iterations before the\n\
		problem is perturbed.  Default = 0 = automatic.",
	perturblimit_desc[] = "Synonym for \"perturblim\".",
	pgradient_desc[] = "Pricing algorithm for primal simplex (default 0):\n\
		   -1 = reduced-cost pricing\n\
		    0 = hybrid reduced-cost and Devex pricing\n\
		    1 = Devex pricing\n\
		    2 = steepest-edge pricing\n\
		    3 = steepest-edge with slack initial norms\n\
		    4 = full pricing.",
#ifdef CPX_PARAM_PRECOMPRESS
	precompress_desc[] = "New in CPLEX 7.1:  whether to compress the original\n\
		problem after CPLEX's presolve:\n\
		   -1 = no\n\
		    0 = automatic choice (default)\n\
		    1 = yes\n\
		Withdrawn in CPLEX 10.0.",
#endif
	predual_desc[] = "Whether CPLEX's presolve phase should present the\n\
		CPLEX solution algorithm with the primal (-1) or\n\
		dual (1) problem or (default = 0) should decide\n\
		which automatically.  Specifying \"predual=1\" often\n\
		gives better performance than specifying just \"dual\",\n\
		but sometimes \"dual predual=1\" is still better.",
	prelinear_desc[] = "Whether CPLEX's presolve should do full reductions\n\
		or only linear ones.  Default = 1 = full.",
#ifdef CPX_PARAM_PREPASS
	prepass_desc[] = "Limit on number of CPLEX presolve passes.\n\
		Default = -1 = decide limit automatically.",
#endif
	prereduce_desc[] = "Kinds of reductions permitted during CPLEX presolve:\n\
		   0 = none\n\
		   1 = only primal\n\
		   2 = only dual\n\
		   3 = both primal and dual (default).",
	presolve_desc[] = "0 or 1 (default 1): Whether to run CPLEX's presolve\n\
		algorithm.",
	presolvedual_desc[] = "Synonym for \"predual\".",
	presolvenode_desc[] = "-1, 0, or 1 (default 0): Whether to run CPLEX's\n\
		presolve at each node of the MIP branch-and-bound\n\
		tree: -1 = no; 1 = yes; 0 = automatic choice.",
	prestats_desc[] = "0 or 1 (default 0):  Whether to include summary\n\
		statistics (if nonzero) for CPLEX's \"aggregate\" and\n\
		\"presolve\" algorithms in the solve_message.",
#ifdef CPX_TUNE_TILIM
	pretunefile_desc[] = "File to which nondefault keyword settings are written\n\
		before tuning; written whether or not tunefile or\n\
		tunefilecpx is specified.",
	pretunefileprm_desc[] = "File to which nondefault keyword settings are written\n\
		in CPLEX PRM format before tuning; written whether or\n\
		not tunefile or tunefileprm is specified.  Includes\n\
		some display settings suppressed by pretunefile.",
#endif
	pricing_desc[] = "Size of pricing candidate list (for partial pricing).\n\
		0 (default) means the algorithm decides.",
	primal_desc[] = "Single-word phrase:  solve the primal problem.",
	primalopt_desc[] = "Use the primal simplex algorithm.",
#ifdef CPXERR_QCP_SENSE
	qcdmax_desc[] = "Limit on k*n*n for computing duals for quadratically\n\
		constrained problems, where k = number of quadratic\n\
		constraints and n = number of variables. Default = 1e9.",
	qcdual_desc[] = "Whether to compute dual variable values for problems\n\
		with quadratic constraints.  Default = 1 (for \"yes\").\n\
		This may be expensive if there are many quadratic\n\
		constraints.  Specifying qcdual=0 suppresses the\n\
		computation.",
#endif
#ifdef CPX_PARAM_BARQCPEPCOMP
	qcpconvergetol_desc[] = "Convergence tolerance on relative complementarity for\n\
		problems with quadratic constraints.  Default = 1e-7.",
#endif
#ifdef CPXERR_QCP_SENSE
	qctol1_desc[] = "Tolerance on a quadratic inequality constraint's slack.\n\
		After CPLEX has returned a solution, dual values are\n\
		deduced for \"active\" quadratic constraints.\n\
		Default 1e-5; a negative value is quietly treated as 0.",
	qctol2_desc[] = "Tolerance on the maxnorm of the gradient of an\n\
		\"active\" quadratic constraint (see qctol1):  if the\n\
		maxnorm is no more than qctol2, the gradient is\n\
		considered to vanish and dual value 0 is deduced.\n\
		Default 1e-5; a negative value is quietly treated as 0.",
	qctol3_desc[] = "Tolerance on the reduction during QR factorization of\n\
		the maxnorm of an \"active\" constraint's gradient\n\
		(see qctol1) for the constraint to be considered\n\
		independent of the other active quadratic constraints.\n\
		Dual value 0 is deduced for dependent constraints.\n\
		Default 1e-5; a negative value is quietly treated as 0.",
#endif
#ifdef CPX_PARAM_QPMETHOD
	qpmethod_desc[] = "Choice of algorithm for a continuous quadratic\n\
		programming problem:\n\
		   0 = automatic choice (default)\n\
		   1 = primal simplex\n\
		   2 = dual simplex\n\
		   3 = network simplex\n\
		   4 = barrier algorithm\n\
		   6 = concurrent optimizer.",
#endif
#ifdef CPX_PARAM_QTOLININD
	qtolin_desc[] = "Whether to to linearize products of bounded variables\n\
			in quadratic objectives:\n\
		   -1 = automatic choice (default)\n\
		    0 = no\n\
		    1 = yes.",
#endif
	rays_desc[] = "Whether to return suffix .unbdd when the objective is\n\
		unbounded or suffix .dunbdd when the constraints are\n\
		infeasible:\n\
		   0 = neither\n\
		   1 = just .unbdd\n\
		   2 = just .dunbdd\n\
		   3 = both (default)\n\
		To get .dunbdd, you may need to specify presolve=0\n\
		in $cplex_options.",
	readbasis_desc[] = "BAS file containing starting basis.",
	readsol_desc[] = "File (previously written by an endsol directive) for\n\
		reading the starting point.  This is for debugging\n\
		and is normally not used.",
#if defined(BARRIER) && !defined(NO_DEPRECATED)
	readvector_desc[] = "VEC file containing starting point for barrier alg.\n\
		Deprecated; use \"readsol\" instead.",
#endif
#ifdef CPX_PARAM_RECORD
	record_desc[] = "Whether to record CPLEX library calls for debugging use\n\
		by IBM in a file with an automatically chosen name of\n\
		the form cplexXXXXXXX.db:\n\
			0 = no (default)\n\
			1 = yes.",
#endif
	refactor_desc[] = "LP iterations between refactorizing the basis.\n\
		0 (default) means the algorithm decides.",
	scale_desc[] = "How to scale the problem:\n\
		   -1 = no scaling\n\
		    0 (default) = equilibration\n\
		    1 = a more aggressive scheme that sometimes helps.",
#ifdef CPX_PARAM_RANDOMSEED
	seed_desc[] = "Seed for random number generator used internally\n\
		by CPLEX.  Use \"seed=?\" to see the default, which\n\
		depends on the CPLEX release.",
#endif
	sensitivity_desc[] = "Single-word phrase:  return sensitivity information for\n\
		the objective (in suffixes .up for the largest value\n\
		of a variable's cost coefficient or constraint's\n\
		right-hand side before the optimal basis changes,\n\
		.down for the smallest such value, and .current for\n\
		the current cost coefficient or right-hand side).",
	siftingopt_desc[] = "Synonym for \"siftopt\".",
	siftopt_desc[] = "Single-word phrase:  on LPs with CPLEX versions >= 8,\n\
		solve ever larger sequences of subproblems until the\n\
		whole LP is solved.",
#ifdef CPX_PARAM_SIFTSIM
	simplexsifting_desc[] = "Whether to allow the simplex algorithm to use sifting\n\
		when appropriate:\n\
		   0 = no\n\
		   1 = yes (default).",
#endif
#ifdef CPX_PARAM_SIMTHREADS
	simthreads_desc[] = "Maximum threads for the dual simplex algorithm\n\
		(currently always 1).  Withdrawn from CPLEX 9.0.",
#endif
	singular_desc[] = "Maximum number of times CPLEX should try to\n\
		repair the basis when it encounters singularities.\n\
		Default = 10.",
	singularlim_desc[] = "Synonym for \"singular\".",
#ifdef CPX_PARAM_SOLUTIONTYPE
	solutiontype_desc[] = "Whether to seek a basic solution when solving an LP:\n\
		   0 = automatic choice (default)\n\
		   1 = yes\n\
		   2 = no (just seek a primal-dual pair).",
#endif
	sparam_desc[] = "Used with syntax \"sparam=n=str\" (no spaces), where n\n\
		is a decimal integer, the number of a CPLEX string\n\
		parameter.  If str is ?, report the current value of\n\
		string parameter n.  Otherwise, if str is a quoted\n\
		string or a sequence of nonblank characters, assign\n\
		str to string parameter n.  This facility provides a\n\
		way to modify string parameters that have not (yet)\n\
		been assigned a keyword.",
#ifdef CPX_PARAM_THREADS
	threads_desc[] = "Default maximum number of threads for any of\n\
		the parallel CPLEX optimizers (limited also\n\
		by licensing).  Default = 1 prior to CPLEX 11,\n\
		or 0 (use maximum threads available) starting\n\
		with CPLEX 11.  May be overridden, prior to\n\
		CPLEX 11, by more specific limits, such as\n\
		barthreads or mipthreads.",
#endif
	time_desc[] = "Time limit in seconds; default = 1e75.",
	timelimit_desc[] = "Synonym for \"time\".",
	timing_desc[] = "Whether to write times in seconds or \"ticks\" to\n\
		stdout or stderr: sum of\n\
		    1 = write time in seconds to stdout\n\
		    2 = write time in seconds to stderr"
#if CPX_VERSION >= 12050000
			"\n\
		    4 = write time in \"ticks\" to stdout\n\
		    8 = write time in \"ticks\" to stderr\n\
		   16 = write number of logical cores to stdout\n\
		   32 = write number of logical cores to stderr"
#endif
		".\n\
		Default = 0.",
	tranopt_desc[] = "Synonym for \"dualopt\".",
#ifdef CPX_TUNE_TILIM
	tunedisplay_desc[] = "How much to print during tuning:\n\
		   0 = nothing\n\
		   1 = minimal printing (default)\n\
		   2 = show parameters being tried\n\
		   3 = exhaustive printing.",
	tunefile_desc[] = "Name of file for tuning results.  If specified, CPLEX\n\
		will experiment with parameter settings that would\n\
		make the solution faster.  This can significantly\n\
		increase execution time of the current invocation, but\n\
		the settings it finds might save time in future runs.",
	tunefileprm_desc[] = "Name of file for tuning results in CPLEX PRM format.\n\
		If specified, CPLEX will experiment with parameter\n\
		settings as described for \"tunefile\".",
	tunefix_desc[] = "List of keywords not to tune, enclosed in quotes\n\
		(\" or ') or separated by commas without white space\n\
		if more than one.",
	tunefixfile_desc[] = "Name of file containing keywords not to tune.\n\
		(There is no PRM format alternative.)  Merged with\n\
		tunefix specification (if any).",
	tunerepeat_desc[] = "How many times to perturb the problem during tuning.\n\
		Default = 1.",
	tunetime_desc[] = "Limit (in seconds) on tuning time; meaningful\n\
		if < time.  Default = 1e75.",
#ifdef CPX_PARAM_TUNINGDETTILIM
	tunetimedet_desc[] = "Limit (in \"ticks\") on tuning time; meaningful\n\
		if < time.  Default = 1e75.",
#endif
#endif
	upperobj_desc[] = "Stop maximizing when the objective value\n\
		goes above upperobj.  Default = 1e75.",
	upperobjlim_desc[] = "Synonym for \"upperobj\".",
	version_desc[] = "Single-word phrase:  show the current version.",
	wantsol_desc[] = "solution report without -AMPL: sum of\n\
		    1 = write .sol file\n\
		    2 = print primal variable values\n\
		    4 = print dual variable values\n\
		    8 = do not print solution message",
#ifdef CPXPARAM_Read_WarningLimit
	warninglimit_desc[] = "Limit on the number of warnings per issue\n\
		given when \"datacheck=2\" is specified.  Default = 10.",
#endif
	workfiledir_desc[] = "Directory where CPLEX creates a temporary\n\
		subdirectory for temporary files, e.g., for\n\
		node information and Cholesky factors.",
#ifdef CPX_PARAM_WORKMEM
	workfilelim_desc[] = "Maximum size in megabytes for in-core work \"files\".\n\
		Default 2048.",
#endif
	writebasis_desc[] = "Synonym for \"endbasis\".",
	writemipstart_desc[] = "[Debug option] The name of a file to which the MIP\n\
		starting guess (if any) is written in \".mst\" format.\n\
		If there is no MIP start, an empty file is written.",
	writeprob_desc[] = "Name of file to which the problem is written\n\
		in a format determined by the name's suffix:\n\
		     .sav = binary SAV file;\n\
		     .mps = MPS file, original names;\n\
		     .lp = LP file, original names;\n\
		     .rmp = MPS file, generic names;\n\
		     .rew = MPS file, generic names;\n\
		     .rlp = LP file, generic names.\n\
		SAV and LP formats are peculiar to CPLEX.",
	writesol_desc[] = "Synonym for \"endsol\"."
#if defined(BARRIER) && !defined(NO_DEPRECATED) && CPX_VERSION < 12070000
	,writevector_desc[] = "Synonym for \"endvector\"."
#endif
#ifdef CPX_PARAM_XXXIND
uxxxstart_desc[] = "Whether to read .xxx files to resume an\n\
		interrupted run of CPLEX (deprecated in 10.0,\n\
		removed in 12.6.1):\n\
		    0 = no (default)\n\
		    1 = yes."
#endif
	;

#define VP (void*)

 static keyword keywds[] = {	/* must be in alphabetical order */

#if 0	/* These can be accessed via iparam and dparam. {*/
	/* Undocumented keywords start with underscore... */

	{ "_aggsort",	sf_int,		VP 1061 /*CPX_PARAM_PREAGGSORT*/, undoc},
	{ "_aggtolerance", sf_dbl,	VP 1055 /*CPX_PARAM_EPSAGG*/, undoc},
	{ "_cancel",	sf_int,		VP 1071 /*CPX_PARAM_PRECANCEL*/, undoc},
	{ "_clique",	sf_int,		VP 1072 /*CPX_PARAM_PRECLIQUE*/, undoc},
	{ "_cliquetablesize", sf_dbl,	VP 2064 /*CPX_PARAM_CLIQUETABLESZ*/, undoc},
	{ "_domination", sf_int,	VP 2038 /*CPX_PARAM_COLDOMIND*/, undoc},
	{ "_effslack",	sf_int,		VP 1042 /*CPX_PARAM_EFFSLACKIND*/, undoc},
	{ "_factormem",	sf_dbl,		VP 3021 /*CPX_PARAM_BARFACTMEM*/, undoc},
	{ "_flip",	sf_int,		VP 1051 /*CPX_PARAM_FLIPIND*/, undoc},
	{ "_hfeasibility", sf_dbl,	VP 1050 /*CPX_PARAM_EPRHS_H*/, undoc},
	{ "_hoptimality", sf_dbl,	VP 1049 /*CPX_PARAM_EPOPT_H*/, undoc},
	{ "_insubtree",	sf_int,		VP 2063 /*CPX_PARAM_INSUBTREE*/, undoc},
	{ "_kernel",	sf_int,		VP 3020 /*CPX_PARAM_BARKERNEL*/, undoc},
	{ "_knapcoeff",	sf_int,		VP 1059 /*CPX_PARAM_KNAPCOERED*/, undoc},
	{ "_localcovers", sf_int,	VP 2061 /*CPX_PARAM_LOCALCOVERS*/, undoc},
	{ "_logparams",	sf_int,		VP 1075 /*CPX_PARAM_LOGPARAMS*/, undoc},
	{ "_memfact",	sf_dbl,		VP 1045 /*CPX_PARAM_PREMEMFACT*/, undoc},
	{ "_memsave",	sf_int,		VP 1060 /*CPX_PARAM_PREMEMSAVE*/, undoc},
	{ "_minstuck",	sf_int,		VP 3023 /*CPX_PARAM_BARMINSTUCK*/, undoc},
	{ "_oldpricing", sf_int,	VP 1054 /*CPX_PARAM_OLDPRICING*/, undoc},
	{ "_oldqpfactor", sf_int,	VP 3024 /*CPX_PARAM_OLDQPFACTOR*/, undoc},
	{ "_oldratio",	sf_int,		VP 1068 /*CPX_PARAM_OLDRATIO*/, undoc},
	{ "_orderthreads", sf_int,	VP 3022 /*CPX_PARAM_BARORDERTHREADS*/, undoc},
	{ "_primalstart", sf_dbl,	VP 3005 /*CPX_PARAM_BARPSTART*/, undoc},
	{ "_probe",	sf_int,		VP 1070 /*CPX_PARAM_PREPROBE*/, undoc},
	{ "_recurseheur", sf_int,	VP 2062 /*CPX_PARAM_RECURSEHEUR*/, undoc},
#ifdef CPX_PARAM_IISIND /* version < 9.2b */
	{ "_rowsdense",	sf_int,		VP 3015 /*CPX_PARAM_BARROWSDEN*/, undoc},
#endif
	{ "_splitrow",	sf_int,		VP 1079 /*CPX_PARAM_PRESPLITROW*/, undoc},
	{ "_svbound",	sf_int,		VP 1069 /*CPX_PARAM_SVBNDSTR*/, undoc},
#endif /*}*/

	/* Documented keywords... */

#ifdef CPLEX_MIP
	{ "absmipgap",	sf_dbl,		VP CPX_PARAM_EPAGAP, absmipap_desc },
#endif
	{ "advance",	sf_int1,	VP CPX_PARAM_ADVIND, advance_desc },
#ifdef CPLEX_MIP
	{ "aggcutlim",	sf_int,		VP CPX_PARAM_AGGCUTLIM, aggcutlim_desc },
#endif
	{ "aggfill",	sf_int2,	VP CPX_PARAM_AGGFILL, aggfill_desc },
	{ "agglim",	sf_int,		VP CPX_PARAM_AGGFILL, agglim_desc },
	{ "aggregate",	sf_int1,	VP CPX_PARAM_AGGIND, aggregate_desc },
	{ "aggtol",	sf_dbl,		VP CPX_PARAM_EPSAGG, aggtol_desc },
	{ "autoopt",	sf_known,	VP set_autoopt, autoopt_desc },
	{ "autopt",	sf_known,	VP set_autoopt, autopt_desc },
#ifdef CPLEX_MIP
#ifdef CPX_PARAM_AUXROOTTHREADS
	{ "auxrootthreads", sf_int,	VP CPX_PARAM_AUXROOTTHREADS, auxrootthreads_desc },
#endif
	{ "backtrack",	sf_dbl,		VP CPX_PARAM_BTTOL, backtrack_desc },
#endif
#ifdef BARRIER
	{ "baralg",	sf_int,		VP CPX_PARAM_BARALG, baralg_desc },
	{ "barcorr",	sf_int,		VP CPX_PARAM_BARMAXCOR, barcorr_desc },
	{ "bardisplay",	sf_int2,	VP CPX_PARAM_BARDISPLAY, bardisplay_desc },
	{ "bargrowth",	sf_dbl,		VP CPX_PARAM_BARGROWTH, bargrowth_desc },
	{ "bariterlim",	sf_int,		VP CPX_PARAM_BARITLIM, bariterlim_desc },
	{ "barobjrange", sf_dbl,	VP CPX_PARAM_BAROBJRNG, barobjrange_desc },
	{ "baropt",	sf_known,	VP set_barrier, baropt_desc },
#ifdef CPX_PARAM_BAROOC
	{ "baroutofcore", sf_int,	VP CPX_PARAM_BAROOC, baroutofcore_desc },
#endif
	{ "barstart",	sf_int,		VP CPX_PARAM_BARSTARTALG, barstart_desc },
	{ "barstartalg", sf_int,	VP CPX_PARAM_BARSTARTALG, barstartalg_desc },
#ifdef CPX_PARAM_BARTHREADS
	{ "barthreads",	sf_int,		VP CPX_PARAM_BARTHREADS, barthreads_desc },
#endif
#ifdef CPX_PARAM_BARVARUP
	{ "barvarup",	sf_dbl,		VP CPX_PARAM_BARVARUP, barvarup_desc },
#endif
#endif /* BARRIER */
#ifdef BASDEBUG
	{ "basdebug",	I_val,		VP &basdebug, basdebug_desc },
#endif
	{ "basis_cond",	sf_mint,	VP set_basis_cond, basis_cond_desc },
#ifdef CPX_PARAM_BASINTERVAL
	{ "basisinterval", sf_int,	VP CPX_PARAM_BASINTERVAL, basisinterval_desc },
#endif
#ifdef CPLEX_MIP
	{ "bbinterval",	sf_int,		VP CPX_PARAM_BBINTERVAL, bbinterval_desc },
#ifdef CPXPARAM_Benders_Strategy
	{ "benders_feascut_tol", sf_dbl,VP CPXPARAM_Benders_Tolerances_feasibilitycut, benders_feascut_desc },
	{ "benders_optcut_tol",	sf_dbl,	VP CPXPARAM_Benders_Tolerances_optimalitycut, benders_optcut_desc },
	{ "benders_strategy", sf_int,	VP CPXPARAM_Benders_Strategy, benders_strategy_desc },
	{ "bendersopt", sf_known,	VP set_benders, bendersopt_desc },
#endif
	{ "bestbound",	sf_known,	VP set_bestbound, bestbound_desc },
	{ "bestnode",	sf_known,	VP set_bestnode, bestnode_desc },
	{ "boundstr",	sf_int,		VP CPX_PARAM_BNDSTRENIND, boundstr_desc },
#ifdef CPX_PARAM_BQPCUTS
	{ "bqpcuts",	sf_int,		VP CPX_PARAM_BQPCUTS, bqpcuts_desc },
#endif
	{ "branch",	sf_int,		VP CPX_PARAM_BRDIR, branch_desc },
	{ "branchdir",	sf_int,		VP CPX_PARAM_BRDIR, branchdir_desc },
	{ "cliquecuts",	sf_int2,	VP CPX_PARAM_CLIQUES, cliquecuts_desc },
	{ "cliques",	sf_int,		VP CPX_PARAM_CLIQUES, cliques_desc },
#endif
	{ "clocktype",	sf_int,		VP CPX_PARAM_CLOCKTYPE, clocktype_desc },
#ifdef CPLEX_MIP
	{ "coeffreduce", sf_int1,	VP CPX_PARAM_COEREDIND, coeffreduce_desc },
#endif
#ifdef BARRIER
	{ "comptol",	sf_dbl,		VP CPX_PARAM_BAREPCOMP, comptol_desc },
#endif
	{ "concurrent",	sf_known,	VP set_concurrentopt, concurrent_desc },
	{ "concurrentopt", sf_known,	VP set_concurrentopt, concurrentopt_desc },
#ifdef CPXPARAM_Conflict_Algorithm
	{ "conflictalg", sf_int,	VP CPXPARAM_Conflict_Algorithm, conflictalg_desc },
#endif
#ifdef CPX_PARAM_CONFLICTDISPLAY
	{ "conflictdisplay", sf_int,	VP CPX_PARAM_CONFLICTDISPLAY, conflictdisplay_desc },
#endif
#ifdef CPLEX_MIP
	{ "covercuts",	sf_int2,	VP CPX_PARAM_COVERS, covercuts_desc },
	{ "covers",	sf_int,		VP CPX_PARAM_COVERS, covers_desc },
#endif
#ifdef CPX_PARAM_CPUMASK
	{ "cpumask",	sf_str,		VP CPX_PARAM_CPUMASK, cpumask_desc },
#endif
	{ "crash",	sf_int,		VP CPX_PARAM_CRAIND, crash_desc },
#ifdef BARRIER
	{ "crossover",	sf_mint,	VP set_crossover, crossover_desc },
#endif
#ifdef CPLEX_MIP
	{ "cutpass",	sf_int,		VP CPX_PARAM_CUTPASS, cutpass_desc },
	{ "cutsfactor",	sf_dbl,		VP CPX_PARAM_CUTSFACTOR, cutsfactor_desc },
#if CPX_VERSION >= 1100
	{ "cutstats",	sf_mint,	VP set_cutstats, cutstats_desc },
#endif
#endif
#if CPX_VERSION >= 12070000 || defined(WANT_DATACHECK)
	{ "datacheck",	sf_int,		VP CPX_PARAM_DATACHECK, datacheck_desc },
#endif
#ifdef BARRIER
	{ "dense",	sf_int2,	VP CPX_PARAM_BARCOLNZ, dense_desc },
	{ "densecol",	sf_int,		VP CPX_PARAM_BARCOLNZ, densecol_desc },
#endif
	{ "dependency",	sf_int1,	VP CPX_PARAM_DEPIND, dependency_desc },
#ifdef CPX_PARAM_DETTILIM
	{ "dettimelim",	sf_dbl,		VP CPX_PARAM_DETTILIM, dettimelim_desc },
#endif
	{ "dgradient",	sf_int,		VP CPX_PARAM_DPRIIND, dgradient_desc },
#ifdef CPLEX_MIP
	{ "disjcuts",	sf_int,		VP CPX_PARAM_DISJCUTS, disjcuts_desc },
#endif
	{ "display",	sf_int2,	VP CPX_PARAM_SIMDISPLAY, display_desc },
	{ "doperturb",	sf_int1,	VP CPX_PARAM_PERIND, doperturb_desc },
	{ "dparam",	sf_dpar,	0, dparam_desc },
#if CPX_VERSION >= 1100
	{ "droptol",	sf_mdbl,	VP set_droptol, droptol_desc },
#endif
	{ "dual",	sf_known,	VP set_dual, dual_desc },
	{ "dualopt",	sf_known,	VP set_dualopt, dualopt_desc },
	{ "dualratio",	sf_mdbl,	VP set_dual_ratio, dualratio_desc },
	{ "dualthresh",	sf_mint,	VP set_dualthresh, dualthresh_desc },
#ifdef CPX_PARAM_EACHCUTLIM
	{ "eachcutlim",	sf_int,		VP CPX_PARAM_EACHCUTLIM, eachcutlim_desc },
#endif
	{ "endbasis",	sf_char,	VP set_endbas, endbasis_desc },
	{ "endsol",	sf_char,	VP set_endsol, endsol_desc },
#ifdef CPLEX_MIP
#ifdef CPX_PARAM_IISIND /* version < 9.2b */
	{ "endtree",	sf_char,	VP set_endtree, endtree_desc },
#endif
#endif
#if defined(BARRIER) && !defined(NO_DEPRECATED) && CPX_VERSION < 12070000
	{ "endvector",	sf_char,	VP set_endvector, endvector_desc },
#endif
	{ "feasibility",	sf_dbl,		VP CPX_PARAM_EPRHS, feasibility_desc },
#ifdef CPX_PARAM_FEASOPTMODE /* >= 9.2b */
	{ "feasopt",	sf_mint,	VP set_feasopt, feasopt_desc },
	{ "feasoptobj",	sf_mint,	VP set_feasoptobj, feasoptobj_desc },
#endif
	{ "file",	sf_char,	VP set_wrtfname, file_desc },
#ifdef CPX_PARAM_FINALFACTOR
	{ "finalfactor", sf_int,	VP CPX_PARAM_FINALFACTOR, finalfactor_desc },
#endif
#ifdef CPLEX_MIP
	{ "flowcuts",	sf_int,		VP CPX_PARAM_FLOWCOVERS, flowcuts_desc },
	{ "flowpathcuts", sf_int,	VP CPX_PARAM_FLOWPATHS, flowpathcuts_desc },
#ifdef CPX_PARAM_FPHEUR
	{ "fpheur",	sf_int,		VP CPX_PARAM_FPHEUR, fpheur_desc },
#endif
#ifndef NO_CPLEX66 /* for versions prior to CPLEX 6.6 */
	{ "fraccand",	sf_int,		VP CPX_PARAM_FRACCAND, fraccand_desc },
	{ "fraccuts",	sf_int,		VP CPX_PARAM_FRACCUTS, fraccuts_desc },
	{ "fracpass",	sf_int,		VP CPX_PARAM_FRACPASS, fracpass_desc },
	{ "fractionalcuts",	sf_int,		VP CPX_PARAM_FRACCUTS, fractionalcuts_desc },
#endif
#endif
#ifdef BARRIER
	{ "growth",	sf_dbl,		VP CPX_PARAM_BARGROWTH, growth_desc },
#endif
#ifdef CPLEX_MIP
	{ "gubcuts",	sf_int,		VP CPX_PARAM_GUBCOVERS, gubcuts_desc },
#ifdef CPX_PARAM_HEURFREQ
	{ "heurfreq",	sf_int,		VP CPX_PARAM_HEURFREQ, heurfreq_desc },
#endif
#ifdef CPX_PARAM_HEURISTIC
	{ "heuristic",	sf_int,		VP CPX_PARAM_HEURISTIC, heuristic_desc },
#endif
#ifdef CPX_PARAM_HEURFREQ
	{ "heuristicfreq", sf_int,	VP CPX_PARAM_HEURFREQ, heuristicfreq_desc },
#endif
#endif
	{ "iisfind",	sf_mint,	VP set_iis, iisfind_desc },
#ifdef CPLEX_MIP
	{ "impliedcuts", sf_int,	VP CPX_PARAM_IMPLBD, impliedcuts_desc },
#endif
#ifdef CPXERR_PARAM_INCOMPATIBLE
	{ "incompat",	sf_mint,	VP set_incompat, incompat_desc },
#endif
#ifdef CPLEX_MIP
	{ "integrality",	sf_dbl,		VP CPX_PARAM_EPINT, integrality_desc },
	{ "intwarntol",	D_val,		VP &intwarn_tol, intwarntol_desc },
#endif
	{ "iparam",	sf_ipar,	0, iparam_desc },
	{ "iterations",	sf_int,		VP CPX_PARAM_ITLIM, iterations_desc },
	{ "iterlim",	sf_int,		VP CPX_PARAM_ITLIM, iterlim_desc },
#ifdef Uselazy
	{ "lazy",	sf_mint,	VP set_lazy, lazy_desc },
#endif
#ifdef CPLEX_MIP
#ifdef CPX_PARAM_LBHEUR
	{ "lbheur",	sf_int,		VP CPX_PARAM_LBHEUR, lbheur_desc },
#endif
#endif
	{ "limitperturb", sf_int2,	VP CPX_PARAM_PERLIM, limitperturb_desc },
#ifdef CPXPARAM_MIP_Cuts_LocalImplied
	{ "localimpliedcuts",	sf_int,	VP CPXPARAM_MIP_Cuts_LocalImplied, localimpliedcuts_desc },
#endif
	{ "logfile",	sf_char,	VP set_logname, logfile_desc },
#ifdef CPLEX_MIP
	{ "lowercutoff", sf_dbl,	VP CPX_PARAM_CUTLO, lowercutoff_desc },
#endif
	{ "lowerobj",	sf_dbl,		VP CPX_PARAM_OBJLLIM, lowerobj_desc },
	{ "lowerobjlim", sf_dbl,	VP CPX_PARAM_OBJLLIM, lowerobjlim_desc },
	{ "lpdisplay",	sf_int2,	VP CPX_PARAM_SIMDISPLAY, lpdisplay_desc },
	{ "lpiterlim",	sf_int2,	VP CPX_PARAM_ITLIM, lpiterlim_desc },
	{ "lptimelim",	sf_dbl2,	VP CPX_PARAM_TILIM, lptimelim_desc },
	{ "markowitz",	sf_dbl,		VP CPX_PARAM_EPMRK, markowitz_desc },
	{ "maximize",	sf_known,	VP set_max, maximize_desc },
#ifdef CPX_PARAM_MCFCUTS
	{ "mcfcuts",	sf_int,		VP CPX_PARAM_MCFCUTS, mcfcuts_desc },
#endif
#ifdef CPX_PARAM_MEMORYEMPHASIS
	{ "memoryemphasis", sf_int,	VP CPX_PARAM_MEMORYEMPHASIS, memoryemphasis_desc },
#endif
	{ "minimize",	sf_known,	VP set_min, minimize_desc },
#ifdef CPLEX_MIP /*{*/
	{ "mipalg",	sf_int,		VP CPX_PARAM_SUBALG, mipalg_desc },
	{ "mipalgorithm", sf_int,	VP CPX_PARAM_SUBALG, mipalgorithm_desc },
	{ "mipbasis",	sf_mint,	VP set_mipbasis, mipbasis_desc },
	{ "mipcrossover", sf_int,	VP CPX_PARAM_BARCROSSALG, mipcrossover_desc },
	{ "mipcuts",	sf_mint,	VP set_mipcuts, mipcuts_desc },
	{ "mipdisplay",	sf_int2,	VP CPX_PARAM_MIPDISPLAY, mipdisplay_desc },
	{ "mipemphasis", sf_int,	VP CPX_PARAM_MIPEMPHASIS, mipemphasis_desc },
	{ "mipgap",	sf_dbl,		VP CPX_PARAM_EPGAP, mipgap_desc },
	{ "mipinterval", sf_int,	VP CPX_PARAM_MIPINTERVAL, mipinterval_desc },
#ifdef CPX_PARAM_MIPKAPPASTATS
	{ "mipkappa",	sf_int,		VP CPX_PARAM_MIPKAPPASTATS, mipkappa_desc },
#endif
	{ "mipordertype", sf_int2,	VP CPX_PARAM_MIPORDTYPE, mipordertype_desc },
#ifdef CPX_PARAM_MIPSEARCH
	{ "mipsearch",	sf_int,		VP CPX_PARAM_MIPSEARCH, mipsearch_desc },
#endif
	{ "mipsolutions", sf_int,	VP CPX_PARAM_INTSOLLIM, mipsolutions_desc },
	{ "mipstart",	sf_mint,	VP set_mipstval, mipstart_desc },
	{ "mipstartalg", sf_int,	VP CPX_PARAM_STARTALG, mipstartalg_desc },
	{ "mipstartstatus", sf_mint,	VP set_mipststat, mipstartstatus_desc },
	{ "mipstartvalue", sf_mint,	VP set_mipstval, mipstartvalue_desc },
	{ "mipsubalg",	sf_int,		VP CPX_PARAM_SUBALG, mipsubalg_desc },
#ifdef CPX_PARAM_MIPTHREADS
	{ "mipthreads",	sf_int,		VP CPX_PARAM_MIPTHREADS, mipthreads_desc },
#endif
#ifdef CPX_PARAM_MIQCPSTRAT
	{ "miqcpstrat",	sf_int,		VP CPX_PARAM_MIQCPSTRAT, miqcpstrat_desc },
#endif
	{ "mircuts",	sf_int,		VP CPX_PARAM_MIRCUTS, mircuts_desc },
#endif /*}*/
#ifndef NO_MOkwf
	{ "modisplay",	sf_int,		VP CPXPARAM_MultiObjective_Display, modisplay_desc },
	{ "multiobj", sf_mint, VP set_multiobj, multiobj_desc },
#endif
	{ "nameround",	sf_mint,	VP set_namernd, nameround_desc },
	{ "netdisplay",	sf_int,		VP CPX_PARAM_NETDISPLAY, netdisplay_desc },
	{ "netfeasibility", sf_dbl,	VP CPX_PARAM_NETEPRHS, netfeasibility_desc },
	{ "netfind",	sf_int,		VP CPX_PARAM_NETFIND, netfind_desc },
	{ "netfinder",	sf_int,		VP CPX_PARAM_NETFIND, netfinder_desc },
	{ "netiterations", sf_int,	VP CPX_PARAM_NETITLIM, netiterations_desc },
	{ "netopt",	sf_mint,	VP set_netopt, netopt_desc },
	{ "netoptimality", sf_dbl,	VP CPX_PARAM_NETEPOPT, netoptimality_desc },
	{ "netpricing",	sf_int,		VP CPX_PARAM_NETPPRIIND, netpricing_desc },
#ifdef CPLEX_MIP
	{ "node",	sf_int2,	VP CPX_PARAM_NODELIM, node_desc },
	{ "nodefile",	sf_int,		VP CPX_PARAM_NODEFILEIND, nodefile_desc },
	{ "nodefiledir", sf_char,	VP set_workfiledir, nodefiledir_desc },
#ifdef CPX_PARAM_NODEFILELIM
	{ "nodefilelim", sf_dbl,	VP CPX_PARAM_NODEFILELIM, nodefilelim_desc },
	{ "nodefilesize", sf_dbl,	VP CPX_PARAM_NODEFILELIM, nodefilesize_desc },
#endif
	{ "nodelim",	sf_int2,	VP CPX_PARAM_NODELIM, nodelim_desc },
	{ "nodes",	sf_int,		VP CPX_PARAM_NODELIM, nodes_desc },
	{ "nodesel",	sf_int,		VP CPX_PARAM_NODESEL, nodesel_desc },
	{ "nodeselect",	sf_int,		VP CPX_PARAM_NODESEL, nodeselect_desc },
#endif
	{ "nosolve",	sf_known,	VP set_nosolve, nosolve_desc },
#ifdef CPX_PARAM_NUMERICALEMPHASIS
	{ "numericalemphasis", sf_int,	VP CPX_PARAM_NUMERICALEMPHASIS, numericalemphasis_desc },
#endif
#ifdef CPLEX_MIP
	{ "objdifference", sf_dbl,	VP CPX_PARAM_OBJDIF, objdifference_desc },
#endif
	{ "objno",	sf_mint,	VP set_objno, objno_desc },
	{ "objrep",	sf_mint,	VP set_objrep, objrep_desc },
	{ "optimality",	sf_dbl,		VP CPX_PARAM_EPOPT, optimality_desc },
	{ "optimize",	sf_known,	VP set_primalopt, optimize_desc },
#ifdef BARRIER
	{ "ordering",	sf_int,		VP CPX_PARAM_BARORDER, ordering_desc },
#endif
#ifdef CPLEX_MIP
	{ "ordertype",	sf_int,		VP CPX_PARAM_MIPORDTYPE, ordertype_desc },
#endif
	{ "outlev",	sf_int2,	VP CPX_PARAM_SIMDISPLAY, outlev_desc },
#ifdef CPX_PARAM_PARALLELMODE
	{ "parallelmode", sf_int,	VP CPX_PARAM_PARALLELMODE, parallelmode_desc },
#endif
	{ "paramfile", sf_par,		VP set_paramfile, paramfile_desc },
#if CPX_VERSION >= 1000
	{ "paramfileprm", sf_parm,	VP set_paramfile, paramfileprm_desc },
#endif
#ifdef CPX_PARAM_PDSWITCH
	{ "pdswitch",	sf_int,		VP CPX_PARAM_PDSWITCH, pdswitch_desc },
#endif
	{ "perturb",	sf_int1,	VP CPX_PARAM_PERIND, perturb_desc },
	{ "perturbation", sf_dbl,	VP CPX_PARAM_EPPER, perturbation_desc },
	{ "perturbconst", sf_dbl,	VP CPX_PARAM_EPPER, perturbconst_desc },
	{ "perturblim",	sf_int,		VP CPX_PARAM_PERLIM, perturblim_desc },
	{ "perturblimit", sf_int,	VP CPX_PARAM_PERLIM, perturblimit_desc },
	{ "pgradient",	sf_int,		VP CPX_PARAM_PPRIIND, pgradient_desc },
#ifdef CPLEX_MIP /*{*/
	{ "plconpri",	sf_mint,	VP set_conpri, plconpri_desc },
	{ "plobjpri",	sf_mint,	VP set_objpri, plobjpri_desc },
#ifdef CPX_PARAM_POLISHAFTEREPAGAP
	{ "polishafter_absmipgap", sf_dbl, VP CPX_PARAM_POLISHAFTEREPAGAP, polishafter_absmipgap_desc },
	{ "polishafter_intsol",	sf_int,	VP CPX_PARAM_POLISHAFTERINTSOL, polishafter_intsol_desc },
	{ "polishafter_mipgap",	sf_dbl,	VP CPX_PARAM_POLISHAFTEREPGAP, polishafter_mipgap_desc },
	{ "polishafter_nodes", sf_int,	VP CPX_PARAM_POLISHAFTERNODE, polishafter_nodes_desc },
	{ "polishafter_time", sf_dbl,	VP CPX_PARAM_POLISHAFTERTIME, polishafter_time_desc },
#endif
#ifdef  CPX_PARAM_POLISHAFTERDETTIME
	{ "polishafter_timedet", sf_dbl, VP CPX_PARAM_POLISHAFTERDETTIME, polishafter_timedet_desc },
#endif
#ifdef CPX_PARAM_POLISHTIME
	{ "polishtime",	sf_dbl,		VP CPX_PARAM_POLISHTIME, polishtime_desc },
#endif
#ifdef CPX_PARAM_POPULATELIM
	{ "poolagap",	sf_dbl,		VP CPX_PARAM_SOLNPOOLAGAP, poolagap_desc },
	{ "poolcapacity", sf_int2,	VP CPX_PARAM_SOLNPOOLCAPACITY, poolcapacity_desc },
	{ "pooldual",	sf_mint,	VP set_pooldual, pooldual_desc },
	{ "poolgap",	sf_dbl,		VP CPX_PARAM_SOLNPOOLGAP, poolgap_desc },
	{ "poolintensity", sf_int,	VP CPX_PARAM_SOLNPOOLINTENSITY, poolintensity_desc },
	{ "poolreplace", sf_int,	VP CPX_PARAM_SOLNPOOLREPLACE, poolreplace_desc },
	{ "poolstub",	sf_char,	VP set_poolstub, poolstub_desc },
	{ "populate",	sf_mint,	VP set_populate, populate_desc },
	{ "populatelim", sf_int,	VP CPX_PARAM_POPULATELIM, populatelim_desc },
#endif
#endif /*}CPLEX_MIP*/
#ifdef CPX_PARAM_PRECOMPRESS
	{ "precompress", sf_int,	VP CPX_PARAM_PRECOMPRESS, precompress_desc },
#endif
	{ "predual",	sf_int,		VP CPX_PARAM_PREDUAL, predual_desc },
	{ "prelinear",	sf_int,		VP CPX_PARAM_PRELINEAR, prelinear_desc },
#ifdef CPX_PARAM_PREPASS
	{ "prepass",	sf_int,		VP CPX_PARAM_PREPASS, prepass_desc },
#endif
	{ "prereduce",	sf_int,		VP CPX_PARAM_REDUCE, prereduce_desc },
#ifdef CPLEX_MIP
	{ "prerelax",	sf_int1,	VP CPX_PARAM_RELAXPREIND, prerelax_desc },
#endif
	{ "presolve",	sf_int1,	VP CPX_PARAM_PREIND, presolve_desc },
	{ "presolvedual", sf_int,	VP CPX_PARAM_PREDUAL, presolvedual_desc },
#ifdef CPLEX_MIP
	{ "presolvenode", sf_int,	VP CPX_PARAM_PRESLVND, presolvenode_desc },
#endif
	{ "prestats",	sf_mint,	VP set_prestats, prestats_desc },
#ifdef CPX_TUNE_TILIM
	{ "pretunefile", sf_char,	VP set_pretunefile, pretunefile_desc },
	{ "pretunefileprm", sf_char,	VP set_pretunefileprm, pretunefileprm_desc },
#endif
	{ "pricing",	sf_int,		VP CPX_PARAM_PRICELIM, pricing_desc },
	{ "primal",	sf_known,	VP set_primal, primal_desc },
	{ "primalopt",	sf_known,	VP set_primalopt, primalopt_desc },
#ifdef CPLEX_MIP
	{ "priorities",	sf_int1,	VP CPX_PARAM_MIPORDIND, priorities_desc },
	{ "probe",	sf_int,		VP CPX_PARAM_PROBE, probe_desc },
#ifdef CPX_PARAM_PROBETIME
	{ "probetime",	sf_dbl,		VP CPX_PARAM_PROBETIME, probetime_desc },
#endif
#ifdef CPX_PARAM_PROBEDETTIME
	{ "probetimedet", sf_dbl,	VP CPX_PARAM_PROBEDETTIME, probetimedet_desc },
#endif
#endif /*CPLEX_MIP*/
#ifdef CPXERR_QCP_SENSE
	{ "qcdmax",	D_val,		VP &qcdmax, qcdmax_desc },
	{ "qcdual",	sf_mint,	VP set_qcdual, qcdual_desc },
#endif
#ifdef CPX_PARAM_BARQCPEPCOMP
	{ "qcpconvergetol", sf_dbl,	VP CPX_PARAM_BARQCPEPCOMP, qcpconvergetol_desc },
#endif
#ifdef CPXERR_QCP_SENSE
	{ "qctol1",	D_val,		VP &qctol1, qctol1_desc },
	{ "qctol2",	D_val,		VP &qctol2, qctol2_desc },
	{ "qctol3",	D_val,		VP &qctol3, qctol3_desc },
#endif
#ifdef CPX_PARAM_QPMETHOD
	{ "qpmethod",	sf_int,		VP CPX_PARAM_QPMETHOD, qpmethod_desc },
#endif
#ifdef CPX_PARAM_QTOLININD
	{ "qtolin",	sf_int,		VP CPX_PARAM_QTOLININD, qtolin_desc },
#endif
#ifdef Want_Distmipopt
	{ "rampup_duration", sf_int,	VP CPX_PARAM_RAMPUPDURATION, rampup_duration_desc },
	{ "rampup_timelim", sf_dbl,	VP CPX_PARAM_RAMPUPDETTILIM, rampup_timelim_desc },
	{ "rampup_walltimelim",	sf_dbl,	VP CPX_PARAM_RAMPUPTILIM, rampup_walltimelim_desc },
#endif
	{ "rays",	sf_mint,	VP set_rays, rays_desc },
	{ "readbasis",	sf_char,	VP set_startbas, readbasis_desc },
	{ "readsol",	sf_char,	VP set_startsol, readsol_desc },
#if defined(BARRIER) && !defined(NO_DEPRECATED)
	{ "readvector",	sf_char,	VP set_startvector, readvector_desc },
#endif
#ifdef CPX_PARAM_RECORD
	{ "record",	sf_int,		VP CPX_PARAM_RECORD, record_desc },
#endif
	{ "refactor",	sf_int,		VP CPX_PARAM_REINV, refactor_desc },
#ifdef CPLEX_MIP
	{ "relax",	sf_known,	VP set_relax, relax_desc },
	{ "relaxpresolve",	sf_int,		VP CPX_PARAM_RELAXPREIND, relaxpresolve_desc },
	{ "relobjdif",	sf_dbl2,	VP CPX_PARAM_RELOBJDIF, relobjdif_desc },
	{ "relobjdiff",	sf_dbl,		VP CPX_PARAM_RELOBJDIF, relobjdiff_desc },
	{ "relpresolve", sf_int,	VP CPX_PARAM_RELAXPREIND, relpresolve_desc },
#ifdef CPX_PARAM_REPAIRTRIES
	{ "repairtries", sf_int,	VP CPX_PARAM_REPAIRTRIES, repairtries_desc },
#endif
#ifdef CPX_PARAM_REPEATPRESOLVE
	{ "repeatpresolve", sf_int,	VP CPX_PARAM_REPEATPRESOLVE, repeatpresolve_desc },
#endif
#ifdef REQCONVEX
	{ "reqconvex",	sf_int,		VP REQCONVEX, reqconvex_desc },
#endif
	{ "resolve",	sf_mint,	VP set_resolve, resolve_desc },
	{ "return_mipgap",	sf_mint,	VP set_retmipgap, return_mipgap_desc },
#ifdef CPX_PARAM_RINSHEUR
	{ "rinsheur",	sf_int,		VP CPX_PARAM_RINSHEUR, rinsheur_desc },
#endif
#ifdef CPXPARAM_MIP_Cuts_RLT
	{ "rltcuts",	sf_int,		VP CPXPARAM_MIP_Cuts_RLT, rltcuts_desc },
#endif
#ifdef CPX_PARAM_HEURISTIC
	{ "rootheuristic",	sf_int,		VP CPX_PARAM_HEURISTIC, rootheuristic_desc },
#endif
	{ "round",	sf_mint,	VP set_round, round_desc },
#endif /*CPLEX_MIP*/
	{ "scale",	sf_int,		VP CPX_PARAM_SCAIND, scale_desc },
#ifdef CPX_PARAM_RANDOMSEED
	{ "seed",	sf_int,		VP CPX_PARAM_RANDOMSEED, seed_desc },
#endif
	{ "sensitivity", sf_known,	VP set_sens, sensitivity_desc },
	{ "siftingopt",	sf_known,	VP set_siftopt, siftingopt_desc },
	{ "siftopt",	sf_known,	VP set_siftopt, siftopt_desc },
#ifdef CPX_PARAM_SIFTSIM
	{ "simplexsifting", sf_int,	VP CPX_PARAM_SIFTSIM, simplexsifting_desc},
#endif
#ifdef CPX_PARAM_SIMTHREADS
	{ "simthreads",	sf_int,		VP CPX_PARAM_SIMTHREADS, simthreads_desc },
#endif
	{ "singular",	sf_int,		VP CPX_PARAM_SINGLIM, singular_desc },
	{ "singularlim", sf_int,	VP CPX_PARAM_SINGLIM, singularlim_desc },
#ifdef CPLEX_MIP /*{*/
	{ "solutionlim", sf_int,	VP CPX_PARAM_INTSOLLIM, solutionlim_desc },
#endif /*}*/
#ifdef CPX_PARAM_SOLUTIONTYPE
	{ "solutiontype", sf_int,	VP CPX_PARAM_SOLUTIONTYPE, solutiontype_desc },
#endif
#ifdef CPLEX_MIP /*{*/
	{ "sos",	sf_mint,	VP set_sos, sos_desc },
	{ "sos2",	sf_mint,	VP set_sos2, sos2_desc },
	{ "sparam", 	sf_spar,	0, sparam_desc },
#ifdef CPX_PARAM_LANDPCUTS
	{ "splitcuts",	sf_int,		VP CPX_PARAM_LANDPCUTS, splitcuts_desc },
#endif
	{ "startalg",	sf_int,		VP CPX_PARAM_STARTALG, startalg_desc },
	{ "startalgorithm",	sf_int,		VP CPX_PARAM_STARTALG, startalgorithm_desc },
	{ "startbasis",	sf_char,	VP set_startbas, startbasis_desc },
	{ "startsol",	sf_char,	VP set_startsol, startsol_desc },
#ifdef CPX_PARAM_IISIND /* version < 9.2b */
	{ "starttree",	sf_char,	VP set_starttree, starttree_desc },
#endif
#if defined(BARRIER) && !defined(NO_DEPRECATED)
	{ "startvector",	sf_char,	VP set_startvector, "synonym for readvector" },
#endif
	{ "strongcand",	sf_int,		VP CPX_PARAM_STRONGCANDLIM, strongcand_desc },
	{ "strongit",	sf_int,		VP CPX_PARAM_STRONGITLIM, strongit_desc },
#ifdef CPX_PARAM_STRONGTHREADLIM
	{ "strongthreads", sf_int,	VP CPX_PARAM_STRONGTHREADLIM, strongthreads_desc },
#endif
	{ "subalg",	sf_int,		VP CPX_PARAM_SUBALG, subalg_desc },
	{ "subalgorithm", sf_int,	VP CPX_PARAM_SUBALG, subalgorithm_desc },
#ifdef CPX_PARAM_SUBMIPSUBALG
	{ "submipalg",	sf_int,		VP CPX_PARAM_SUBMIPSUBALG, submipalg_desc },
#endif
#ifdef CPXPARAM_MIP_SubMIP_NodeLimit
	{ "submipnodelim", sf_int,	VP CPXPARAM_MIP_SubMIP_NodeLimit, submipnodelim_desc },
#endif
#ifdef CPX_PARAM_SUBMIPSCAIND
	{ "submipscale", sf_int,	VP CPX_PARAM_SUBMIPSCAIND, submipscale_desc },
#endif
#ifdef CPX_PARAM_SUBMIPSTARTALG
	{ "submipstart", sf_int,	VP CPX_PARAM_SUBMIPSTARTALG, submipstart_desc },
#endif
#ifdef CPX_PARAM_SYMMETRY
	{ "symmetry",	sf_int,		VP CPX_PARAM_SYMMETRY, symmetry_desc },
#endif
#endif /*}*/
#ifdef CPX_PARAM_THREADS
	{ "threads",	sf_int,		VP CPX_PARAM_THREADS, threads_desc },
#endif
	{ "time",	sf_dbl,		VP CPX_PARAM_TILIM, time_desc },
	{ "timelimit",	sf_dbl,		VP CPX_PARAM_TILIM, timelimit_desc },
	{ "timing",	sf_mint,	VP set_timing, timing_desc },
	{ "tranopt",	sf_known,	VP set_dualopt, tranopt_desc },
#ifdef CPLEX_MIP
	{ "treelimit",	sf_dbl2,	VP CPX_PARAM_TRELIM, treelimit_desc },
	{ "treememlim",	sf_dbl2,	VP CPX_PARAM_TRELIM, treememlim_desc },
	{ "treememory",	sf_dbl,		VP CPX_PARAM_TRELIM, treememory_desc },
#endif
#ifdef CPX_TUNE_TILIM
	{ "tunedisplay", sf_int,	VP CPX_PARAM_TUNINGDISPLAY, tunedisplay_desc },
	{ "tunefile",	sf_char,	VP set_tunefile, tunefile_desc },
	{ "tunefileprm", sf_char,	VP set_tunefileprm, tunefileprm_desc },
	{ "tunefix",	sf_char,	VP set_tunefix, tunefix_desc },
	{ "tunefixfile", sf_char,	VP set_tunefixfile, tunefixfile_desc },
	{ "tunerepeat",	sf_int,		VP CPX_PARAM_TUNINGREPEAT, tunerepeat_desc },
	{ "tunetime",	sf_dbl,		VP CPX_PARAM_TUNINGTILIM, tunetime_desc },
#ifdef CPX_PARAM_TUNINGDETTILIM
	{ "tunetimedet", sf_dbl,	VP CPX_PARAM_TUNINGDETTILIM, tunetimedet_desc },
#endif
#endif
#ifdef CPLEX_MIP
	{ "uppercutoff", sf_dbl,	VP CPX_PARAM_CUTUP, uppercutoff_desc },
#endif
	{ "upperobj",	sf_dbl,		VP CPX_PARAM_OBJULIM, upperobj_desc },
	{ "upperobjlim", sf_dbl,	VP CPX_PARAM_OBJULIM, upperobjlim_desc },
#ifdef CPLEX_MIP
	{ "varsel",	sf_int,		VP CPX_PARAM_VARSEL, varsel_desc },
	{ "varselect",	sf_int,		VP CPX_PARAM_VARSEL, varselect_desc },
#endif
	{ "version",	Ver_val,	0, version_desc },
#ifdef Want_Distmipopt
	{ "vmconf",	sf_char,	VP set_vmconfig, vmconf_desc },
#endif
	{ "wantsol",	WS_val,		0, wantsol_desc },
#ifdef CPXPARAM_Read_WarningLimit
	{ "warninglimit", sf_int,	VP CPXPARAM_Read_WarningLimit, warninglimit_desc },
#endif
	{ "workfiledir", sf_char,	VP set_workfiledir, workfiledir_desc },
#ifdef CPX_PARAM_WORKMEM
	{ "workfilelim", sf_dbl,		VP CPX_PARAM_WORKMEM, workfilelim_desc },
#endif
	{ "writebasis",	sf_char,	VP set_endbas, writebasis_desc },
	{ "writemipstart", sf_char,	VP set_wrtmipstart, writemipstart_desc },
	{ "writeprob",	sf_char,	VP set_wrtfname, writeprob_desc },
	{ "writesol",	sf_char,	VP set_endsol, writesol_desc },
#if defined(BARRIER) && !defined(NO_DEPRECATED) && CPX_VERSION < 12070000
	{ "writevector", sf_char,	VP set_endvector, writevector_desc },
#endif
#ifdef CPX_PARAM_XXXIND
	{ "xxxstart",	sf_int,		VP CPX_PARAM_XXXIND, xxxstart_desc },
#endif
#ifdef CPX_PARAM_ZEROHALFCUTS
	{ "zerohalfcuts", sf_int,	VP CPX_PARAM_ZEROHALFCUTS, zerohalfcuts_desc },
#endif
	};

 static Option_Info Oinfo = { "cplex", 0, "cplex_options", keywds, nkeywds,
				ASL_OI_keep_underscores, cplex_version, 0,
				MOkwf,0,0,0, 20190501, 0,0,0,0,0,0,0,
				ASL_OI_tabexpand | ASL_OI_addnewline };

 static void
badlic(int rc, int status)
{
	char buf[4096];

	if (!CPXgeterrorstring(Env, status, buf))
		Snprintf(buf, sizeof(buf),
		"CPLEX licensing problem: error code %d from CPXopenCPLEX.",
			status);
	badretfmt(rc, "%s", buf);
	}

 static void
nonlin(int n, int rc, char *what)
{
	if (n) {
		badretfmt(rc, "%s contains %s.\n", asl->i.filename_, what);
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
#if !defined(NO_DEPRECATED) && CPX_VERSION < 12070000
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
#if CPX_VERSION < 12070000 /*{*/
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
#endif /*}*/
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
qmatadj(int k, int nr, int os, CPXNNZ *colq, int *colqcnt, double **qmatp, char *ctype)
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
			if (t*osd < 0. && (!ctype || ctype[k] == 'C')) {
				if (nbad < 3) {
					badd[nbad] = t;
					badk[nbad] = k;
					}
				nbad++;
				}
			}
		qmat[k] = t;
		}
#ifdef REQCONVEX
	if (nbad && parval(REQCONVEX) >= 2)
		nbad = 0;
#endif
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

#ifdef CPXERR_QCP_SENSE /*{*/

 static void
Surprise(ASL *asl, int n, char *who)
{
	char msgbuf[64];

	sprintf(msgbuf, "Surprise return %d from %s.", n, who);
	nonlin(1, 561, msgbuf);
	}

 static void**
linadj(ASL *asl, CStype *ka, int *ia, double *a, int **kaqp, int **kalqp, int **iaqp, double **aqp,
	CPXNNZ **nelqcp, int ***rowqcp, int ***colqcp, double ***qcmatp, dims *d, void **v)
{
	/* Move linear parts of nqc quadratic constraints from ka, ia, a */
	/* to kaq, kalq, iaq, aq and corresponding Cgrad pointers and return a */
	/* pointer to the relevant allocation. */

	CPXNNZ *nelqc;
	QPinfo *qpi;
	cde *consave;
	cgrad **cgp, *cg, *cga, *cg0;
	char *what;
	double *aq, **qcmat, *x, *y;
	int i, i1, i2, j, k, m, m0, m1, n, ncol, nqc;
	int *colno, **colqc, *iaq, *kalq, *kaq, *kaq0, *rowno, **rowqc, *z, *z1;
	size_t *colbeg, is, js, ks, ms, nelq, nl, nqcnl;
	ssize_t nq;

	m = n_con;
	m1 = asl->i.n_con1;
	n = n_var;
	nqc = nlc;

	d->consave = consave = (cde*)M1alloc(m1*sizeof(cde));
	memcpy(consave, ((ASL_fg*)asl)->I.con_de_, m1*sizeof(cde));
	kaq0 = (int*)Malloc(nqc*sizeof(int));
	if ((cgp = Cgrad)) {
		nl = 0;
		for(i = 0; i < nqc; ++i) {
			for(j = 0, cg = cgp[i]; cg; cg = cg->next)
				++j;
			kaq0[i] = j;
			nl += j;
			}
		cga = (cgrad*)M1alloc(nl*sizeof(cgrad));
		}
	else {
		js = nzc;
		for(is = nl = 0; is < js; ++is)
			if (ia[is] < nqc)
				nl++;
		cg = cga = (cgrad*)M1alloc(m1*sizeof(cgrad*) + nl*sizeof(cgrad));
		cgp = (cgrad**)(cga + nl);

		memset(kaq0, 0, nqc*sizeof(int));
		memset(cgp, 0, m*sizeof(cgrad*));
		Cgrad = cgp;
		for(js = i = 0; i < n; i++) {
			for(ks = ka[i+1]; js < ks; js++)
				if ((i1 = ia[js]) < nqc) {
					kaq0[i1]++;
					cg->varno = i;
					cg->coef = a[js];
					cg->next = cgp[i1];
					cgp[i1] = cg++;
					}
			}
		}
	x = LUrhs;
	y = Urhsx;
	for(nl = nqcnl = i = 0; i < nqc; i++) {
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
		nq = mqpcheckv(-(i+1), 0, v);
		if (nq <= 0) {
			what = con_name(i);
			if (nq == -2)
				badretfmt(558, "Constraint %s is a quadratic constraint "
					"involving division by 0.", what);
			else if (nq < 0)
				badretfmt(550, "Constraint %s is a nonquadratic "
						"nonlinear constraint.", what);
			else
				badretfmt(559, "CPLEX driver bug: no quadratic terms "
					"in \"nonlinear\" constraint %s.", what);
			exit(4);
			}
		nqcnl += nq;
		}

	*aqp = aq = (double*)Realloc(kaq0, nqc*(sizeof(double*) + 2*sizeof(int*))
				+ (2*(nqc + nqcnl) + nl)*sizeof(int)
				+ nqc*sizeof(CPXNNZ)
				+ (nl+nqcnl)*sizeof(double));
	x = aq + nl;
	*qcmatp = qcmat = (double**)(x + nqcnl);
	*rowqcp = rowqc = (int**)(qcmat + nqc);
	*colqcp = colqc = rowqc + nqc;
	*nelqcp = nelqc = (CPXNNZ*)(colqc + nqc);
	*kaqp = kaq = (int*)(nelqc + nqc);
	*kalqp = kalq = kaq + nqc;
	z = kalq + nqc;
	*iaqp = iaq = z + 2*nqcnl;

	kaq0 = (int*)aq;
	for(i = nqc; i > 0; ) {
		--i;
		kaq[i] = kaq0[i];
		}

	memset(cgp, 0, nqc*sizeof(cgrad*));

	for(ms = i = j = 0; j < n; j++) {
		ka[j] = ms;
		for(k = ka[j+1]; i < k; i++)
			if ((i1 = ia[i]) < nqc) {
				iaq[i2 = kaq[i1]++] = j;
				aq[i2] = a[i];
				}
			else {
				ia[ms] = i1 - nqc;
				a[ms++] = a[i];
				}
		}
	nzc = ka[n] = ms;
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
		nelqc[i] = nelq = mqpcheckv(-(i+1), &qpi, v);
		qcmat[i] = x;
		y = qpi->delsq;
		for(js = 0; js < nelq; ++js)
			*x++ = 0.5*y[js];
		rowqc[i] = z1 = z;
		colqc[i] = z += nelq;
		rowno = qpi->rowno;
		colbeg = qpi->colbeg;
		colno = qpi->colno;
		ncol = qpi->nc;
		is = colbeg[0];
		for(j = 0; j < ncol;) {
			k = colno[j];
			for(js = colbeg[++j]; is < js; ++is) {
				*z1++ = rowno[is];
				*z++ = k;
				}
			}
		free(qpi);
		kaq[i] = m0 = m;
		for(cg = cgp[i]; cg; cg = cg->next)
			if (cg->coef) {
				iaq[m] = cg->varno;
				aq[m++] = cg->coef;
				}
		kalq[i] = m - m0;
		}

	return M1record((void*)aq);
	}

 static real
vmax(real *x, size_t n)	/* returns max{i in 0..n-1} |x[i]|, assuming n >= 1 */
{
	real mn, mx, t, *xe;

	mn = mx = *x;
	xe = x + n;
	while(++x < xe) {
		if (mn > (t = *x))
			mn = t;
		else if (mx < t)
			mx = t;
		}
	if (mx < 0.)
		mx = -mn;
	else if ((t = -mn) > mx)
		mx = t;
	return mx;
	}

 static real
dot(real *x, real *y, size_t n)
{
	real rv, *xe;

	xe = x + n;
	rv = 0.;
	while(x < xe)
		rv += *x++ * *y++;
	return rv;
	}

 static real
ssq(real *x, size_t n)
{
	real rv, t;
	size_t i;

	rv = 0.;
	for(i = 0; i < n; ++i) {
		t = x[i];
		rv += t*t;
		}
	return rv;
	}

 static void
LSsol(size_t n, size_t m, real *A, real *x, real *b, int *bt)
{
	/* Solve A*x = b in the least-squares sense, with A consisting of m */
	/* columns of length n, stored column-wise, subject to x[i] >= 0 */
	/* if bt[i] > 0 and to x[i] <= 0 if bt[i] < 0. */
	/* We clobber A and b. */

	/* We use column pivoting and Householder transformations. */
	/* Given vector x = (x1, x2)^T with x1 a scalar and x2 a vector, */
	/* H = I - 2*y*y^T/y^T*y with y = (x1 + sigma*norm(x), x2)^T, sigma = +-1, */
	/* norm(x) = sqrt(x^T*x) is a reflection (Householder transformation) */
	/* with H*x = (-sigma*norm(x), 0)^T. */
	/* Note that y^T*y = 2*(x^T*x + sigma*x1*norm(x)). */
	/* For stability, sigma = sign(x1).  */

	int move;
	real *ai, *aie, *aj, *bnew, di, *rd, *scale, t, t2, tb, x1, *xnew;
	size_t i, j, jb, k, k1, k2, m1, ni, ni1, nq, *p, pj, *q, *ql;

	if ((k = m) > n)
		k = n;
	if (k <= 0) {
		for(i = 0; i < m; ++i)
			x[i] = 0.;
		return;
		}
	rd = (real*)Malloc((3*m + n)*sizeof(real) + (m+2*k)*sizeof(size_t));
	scale = rd + m;
	bnew = scale + m;
	xnew = bnew + n;
	p = (size_t*)(xnew + m);
	q = p + m;
	ql = q + k;
	for(i = 0; i < m; ++i) {
		p[i] = i;
		x[i] = 0.;
		}
	ai = A + n*m;
	do {
		--i;
		t = vmax(ai -= n, n);
		if (t <= qctol2)
			p[i] = p[--m];
		else if (t != 1.) {
			t = 1. / t;
			for(j = 0; j < n; ++j)
				ai[j] *= t;
			}
		scale[i] = t;
		} while(i > 0);
	m1 = m;
	k1 = move = 0;
	while(i < n) {
 findbest:
		/* find "best" column to process next */
		for(jb = i;;) {
			pj = p[jb];
			tb = dot(A + n*pj, b, n);
			if (bt[pj] < 0 || (bt[pj] == 0 && tb < 0.))
				tb = -tb;
			if (tb > 0.) {
				ai = A + n*pj;
				for(j = i; j < n; ++j)
					if (ai[j] != 0.)
						goto have_jb0;
				}
			if (++jb >= m1) {
				if (!move || m1 >= m)
					goto break2;
				m1 = m;
				move = 0;
				}
			}
		/* now jb = candidate "best" next column; look for better ones */
 have_jb0:
		for(j = jb; ++j < m1; ) {
			pj = p[j];
			t = dot(A + n*pj, b, n);
			if (bt[pj] < 0 || (bt[pj] == 0 && t < 0.))
				t = -t;
			if (tb < t) {
				tb = t;
				jb = j;
				}
			}
		/* now jb = "best" next column */
		pj = p[jb];
		if (jb > i) {
			p[jb] = p[i];
			p[i] = pj;
			}
		ni = n - i;
		ni1 = ni - 1;
		aie = A + n*pj;
		ai = aie + i;
		aie += n;
		t = sqrt(t2 = ssq(ai, ni));
		if ((x1 = *ai) < 0.)
			t = -t;
		/* now t == sigma*norm(x) */
		di = 1. / (t2 + x1*t);
		rd[i] = -t;
		*ai += t;
		/* apply the just-constructed Householder transformation */
		for(j = i + 1; j < m; ) {
			pj = p[j];
			aj = A + i + n*pj;
			if ((t = dot(ai, aj, ni))) {
				t *= di;
				for(k = 0; k < ni; ++k)
					aj[k] -= t*ai[k];
				if (ni1 > 0) {
					t2 = vmax(aj+1, ni1);
					if (t2 <= qctol3) {
						p[j] = p[--m];
						p[m] = pj;
						continue;
						}
					}
				}
			++j;
			}
		/* apply the Householder to the rhs */
		aj = b + i;
		if ((t = dot(ai, aj, ni))) {
			t *= di;
			for(k = 0; k < ni; ++k)
				aj[k] -= t*ai[k];
			}
		/* update column i -- later we may move x[i] to 0 and still later */
		/* we may need to apply future Householder transforms to it */
		*ai = rd[i];
		while(++ai < aie)
			*ai = 0.;

		/* compute xnew */

		for(j = 0; j <= i; ++j)
			bnew[j] = b[j];
		for(;;) {
			pj = p[--j];
			t = bnew[j] /= rd[j];
			xnew[pj] = x[pj] + t*scale[pj];
			if (j == 0)
				break;
			if (t) {
				aj = A + n*pj;
				k = j;
				do {
					--k;
					bnew[k] -= t*aj[k];
					}
					while(k > 0);
				}
			}

		/* see if xnew is feasible */

		for(jb = 0; jb < i; ++jb) {
			pj = p[jb];
			if (xnew[pj] * bt[pj] < 0.) {
				/* no: find the first x component to be zeroed */
				/* as we move from x to xnew */
				tb = x[pj] / (x[pj] - xnew[pj]);
				for(j = jb; ++j < i; ) {
					pj = p[j];
					if (xnew[pj] * bt[pj] < 0.) {
						t = x[pj] / (x[pj] - xnew[pj]);
						if (tb > t) {
							tb = t;
							jb = j;
							}
						}
					}
				if (tb > 0.)
					++move;

				/* move tb of the way from x to xnew, */
				/* noting nq zeroed x components in q */

				nq = 0;
				for(j = k = 0; j <= i; ++j) {
					pj = p[j];
					x[pj] += tb * (xnew[pj] - x[pj]);
					if (j == jb || (bt[pj] && bt[pj]*x[pj] <= 0.)) {
						if (!nq)
							k1 = j;
						q[nq++] = pj;
						}
					else {
						p[k] = pj;
						ql[k++] = j;
						}
					t = tb * bnew[j];
					aj = A + pj*n;
					for(k2 = 0; k2 <= j; ++k2)
						b[k2] -= t*aj[k2];
					}
				while(j < m1)
					p[k++] = p[j++];
				for(j = 0; j < nq; ++j)
					p[k++] = q[j];
				m1 -= nq;
				i -= nq;

				/* bring R matrix back to triangular form */

				for(j = k1; j < i; ++j) {
					pj = p[j];
					ai = A + pj*n + j;
					ni = ql[j] - j + 1;
					t = sqrt(t2 = ssq(ai, ni));
					if ((x1 = *aj) < 0.)
						t = -t;
					di = 1. / (t2 + x1*t);
					rd[j] = -t;
					for(k = j + 1; k < m; ++k) {
						aj = A + n*p[k] + j;
						if ((t = dot(ai, aj, ni))) {
							t *= di;
							for(k1 = 0; k1 < ni; ++k1)
								aj[k1] -= t*ai[k1];
							}
						}
					aj = b + j;
					if ((t = dot(ai, aj, ni))) {
						t *= di;
						for(k1 = 0; k1 < ni; ++k1)
							aj[k1] -= t*ai[k1];
						}
					*ai = rd[j];
					for(k1 = 1; k1 < ni; ++k1)
						ai[k1] = 0.;
					}
				goto findbest;
				}
			}
		/* we can take the full step */
		++move;
		for(j = 0; j <= i; ++j) {
			pj = p[j];
			x[pj] = xnew[pj];
			aj = A + pj*n;
			if ((t = bnew[j])) /* the increment */
				for(k2 = 0; k2 <= j; ++k2)
					b[k2] -= t*aj[k2];
			}
		if (++i >= n)
			break;
		if (i >= m1) {
			if (m1 >= m)
				break;
			move = 0;
			m1 = m;
			}
		}
 break2:
	free(rd);
	}

 static void
qcduals(ASL *asl, PBuf *B, cpxlp *cpx, dims *d, real *x, real *y)
{
	CStype *cs;
	cde *Csave, *Osave;
	int bs, *bt, gmap, i1, k, *p, *q, *qinv, *rn;
	ograd **Ogsave;
	real *A, *A1, *Av, *L, *U, *b, *c, *r, t, tm, tn, tw;
	size_t Ls, i, j, je, m, m1, m2, n, n0, n1;
	ssize_t on;

	tm = n_conjac[1] = m = nlc;
	tn = n = n_var;
	t = (tn*(tm+2) + tm)*sizeof(real) + tn*(3*sizeof(int)) + tm*sizeof(int);
	Ls = (n*(m+2) + m)*sizeof(real) + n*(3*sizeof(int)) + m*sizeof(int);
	if (Ls != t)
		return;
	if ((tw = tn*tm*tm) > qcdmax) {
		if (solve_result_num == 0)
			solve_result_num = 5;
		Bpf(B, "\nComputing of dual variables for quadratic constraints suppressed since "
			"qcdmax = %.g < k*n*n = %.g\n(k = number of quadratic constraints = %d,"
			" n = number of variables = %d).", qcdmax, tw, m, n);
		return;
		}
	A = (real*)Malloc(Ls);
	b = A + n*m;
	r = b + n;
	c = r + n;
	p = (int*)(c+m);
	q = p + n;
	qinv = q + n;
	bt = qinv + n;

	L = LUv;
	U = Uvx;
	for(i = n0 = n1 = 0; i < n; ++i) {
		if ((t = x[i]) > L[i] && t < U[i]) {
			qinv[i] = n1;
			q[n1++] = i;
			}
		else
			qinv[i] = -1;
		}
	if ((gmap = n1 < n)) {
		for(; n0 < n1 && q[n0] == n0; ++n0);
		}
	Csave = ((ASL_fg*)asl)->I.con_de_;
	Osave = ((ASL_fg*)asl)->I.obj_de_;
	((ASL_fg*)asl)->I.con_de_ = d->consave;
	((ASL_fg*)asl)->I.obj_de_ = d->objsave;
	qp_opify();
	bs = 0;
	if ((on = obj_no) >= 0) {
		bs = objtype[on] ? -1 : 1;
		Ogsave = Ograd;
		Ograd = d->ogsave;
		objgrd(on, x, b, 0);
		Ograd = Ogsave;
		if (gmap)
			for(i = n0; i < n1; ++i)
				b[i] = b[q[i]];
		}
	else
		memset(b, 0, n*sizeof(real));
	if ((m2 = n_con - m)) {
		cs = A_colstarts;
		rn = A_rownos;
		Av = A_vals;
		for(i = 0; i < n; ++i) {
			if ((i1 = qinv[i]) >= 0) {
				j = cs[i];
				t = b[i1];
				for(je = cs[i+1]; j < je; ++j) {
					k = rn[j];
					t -= y[k]*Av[j];
					}
				b[i1] = t;
				}
			}
		for(i = m2; i-- > 0;)
			y[i+m] = y[i];
		}
	conval(x, c, 0);
	A1 = A;
	L = LUrhs;
	U = Urhsx;
	if (qctol1 < 0.)
		qctol1 = 0.;
	for(i = m1 = 0; i < m; ++i) {
		t = c[i];
		if (t <= L[i] + qctol1) {
			k = bs;
			if (t >= U[i] - qctol1)
				k = 0;
			}
		else if (t >= U[i] - qctol1)
			k = -bs;
		else {
			y[i] = 0.;
			continue;
			}
		congrd(i, x, A1, 0);
		if (gmap) {
			for(j = n0; j < n1; ++j)
				A1[j] = A1[q[j]];
			}
		A1 += n1;
		p[m1] = i;
		bt[m1++] = k;
		}
	((ASL_fg*)asl)->I.con_de_ = Csave;
	((ASL_fg*)asl)->I.obj_de_ = Osave;
	if (m1) {
		if (qctol2 < 0.)
			qctol2 = 0.;
		if (qctol3 < 0.)
			qctol3 = 0.;
		LSsol(n1, m1, A, r, b, bt);
		for(i = 0; i < m1; ++i)
			y[p[i]] = r[i];
		}
	free(A);
	}
#endif /*CPXERR_QCP_SENSE }*/

#ifdef CPXPARAM_Benders_Strategy /*{*/
 static void
benders_sufcheck(ASL *asl, cpxlp *cpx)
{
	CPXLONG *Z;
	SufDesc *b;
	int a, *ind, i, j, k, n, nm, ns, *z, zi;

	k = 0;
	CPXgetintparam(Env, CPXPARAM_Benders_Strategy, &k);
	if (k < 0 || k >= 3)
		return;
	b = suf_get("benders", ASL_Sufkind_var);
	if (!(z = b->u.i)) {
		if (k) {
			fprintf(Stderr,
				"Ignoring bendersopt because suffix benders is not present\n"
				"and benders_strategy = %d, which requires .benders.\n", k);
			Optimize = CPXmipopt;
			}
		return;
		}
	n = n_var;
	for(i = nm = ns = 0; i < n; ++i) {
		if (z[i])
			++ns;
		else
			++nm;
		}
	if (!nm || !ns) {
		fprintf(Stderr, "Ignoring .benders because %s variables have .benders = 0.\n",
			nm ? "all" : "no");
		if (k)
			Optimize = CPXmipopt;
		return;
		}
	if ((j = CPXnewlongannotation(Env, cpx, CPX_BENDERS_ANNOTATION, 0)))
		fprintf(Stderr, "\n*** CPXnewlongannotation() returned %d\n", j);
	a = 0;
	j = CPXgetlongannotationindex(Env, cpx, CPX_BENDERS_ANNOTATION, &a);
	if (j) {
		fprintf(Stderr,
		"\n*** CPXgetlongannotationindex(CPX_BENDERS_ANNOTATION) returned %d\n", j);
		return;
		}
	Z = (CPXLONG*)Malloc(ns*(sizeof(int) + sizeof(CPXLONG)));
	ind = (int*)(Z + ns);
	for(i = j = 0; i < n; ++i) {
		if ((zi = z[i]) > 0) {
			ind[j] = i;
			Z[j++] = zi;
			}
		}
	j = CPXsetlongannotations(Env, cpx, a, CPX_ANNOTATIONDATA_LONG, ns, ind, Z);
	if (j)
		fprintf(Stderr, "\n*** CPXsetlongannotations returned %d\n", j);
	free(Z);
	}
#endif /*}*/

#ifdef Uselazy
 typedef struct
LazyInfo {
	CStype nz[2], *rowbeg[2];
	char *sense[2], **rowname[2];
	int *colno[2], rows[2];
	double *matval[2], *rhs[2];
	} LazyInfo;

 void **
lazyadj(ASL *asl, LazyInfo *LI, int n, int nint, int nqc, int *mp, double *b, char *senx,
	CStype *ka, int *kal, int *ia, double *a, double *rngvec, char **rname)
{
	CStype i, i2, j, j0, j1, k, m, m0, nq, nrt, nz;
	CStype *rowbeg;
	SufDesc *lp;
	char *nn, *nn1, *sense, **rowname;
	const char *s1, *s2;
	double *matval, *rhs;
	int *colno, i1, nr[3], nt[2], *rl, *z, *z0;
	size_t Lnn, Lrl;
	void *rv;

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
	matval = (double*)Malloc((nz+nrt)*sizeof(double)
				+ nz*sizeof(int)
				+ (nrt+2)*sizeof(CStype)
				+ Lrl*sizeof(int) + nrt + Lnn);
	rv = (void*)matval;
	rhs = matval + nz;
	if (rname) {
		rowname = (char**)(rhs + nrt);
		rowbeg = (CStype*)(rowname + nrt);
		}
	else {
		rowname = 0;
		rowbeg = (CStype*)(rhs + nrt);
		}
	colno = (int*)(rowbeg + nrt + 2);
	rl = colno + nz;
	sense = (char*)(rl + Lrl);
	nn = sense + nrt;
	memset(rowbeg, 0, nrt*sizeof(CStype));
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
		if ((i1 = z[i]) > 0) {
			rhs[--i1] = b[i];
			if ((sense[i1] = senx[i]) == 'R') {
				sense[i1] = 'G';
				sense[++i1] = 'L';
				rhs[i1] = b[i] + rngvec[i];
				if (rname) {
					nn1 = strcpy1(rowname[i1] = nn, rname[i]);
					nn = strcpy(nn1, ".u");
					}
				}
			}
		else {
			b[i1 = -i1] = b[i];
			rngvec[i1] = rngvec[i];
			senx[i1] = senx[i];
			}
		}

	if (rname) {
		LI->rowname[0] = rowname;
		LI->rowname[1] = rowname + nr[0];
		for(i = 0; i < m; i++)
			if ((i1 = z[i]) > 0)
				rowname[i1-1] = rname[i];
			else
				rname[-i1] = rname[i];
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
	CStype *ka, *ka1, *ka2, nz, nz1, nzc0, nzr;
	int *ia, *ia1, *ia2, *ja, *kal, *kal1, *kal2, *vmi;
	double *a, *a1, *b, *b1, *c, *c1, *l, *l1,
		*le, *lrhs, *lx, *rngvec, *rv1, *u, *u1, *urhs, *ux;
	double nos, os;
	ograd *og, *og1, **ogp;
	int i, j, j0, k, k0, m, m0, mq, n, n0, nbv1, nextra, nint, nint1,
		no, nqc, nr, nsos, nzextra;
	char **cname, **rname, *s, *senx;
	void **zap[5];
	size_t L, nelqf;
#ifdef BARRIER
	CPXNNZ *colq;
	int *colno, *colqcnt;
	int havestart, *rowq;
	size_t *colqf;
	ssize_t nq;
	double *cdual, *cprim, *qmat, *rdual, *rprim;
#undef Want_mqpcheck
#ifdef BARRIER_FOR_QP
#ifdef BARRIER
#define Want_mqpcheck
#endif
#else /*!BARRIER_FOR_QP, i.e., >= 8.0*/
#define Want_mqpcheck
#endif
#ifdef Want_mqpcheck
	QPinfo *qpi;
	void *v;
#else
#undef CPXERR_QCP_SENSE
#endif
#ifdef CPXERR_QCP_SENSE	/* CPLEX version >= 9.0 */
	CPXNNZ *nelqc;
	double *aq, **qcmat;
	int **colqc, *iaq, *kaq, *kalq, **rowqc;
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
	static int repmap[4] = { 0, ASL_obj_replace_ineq,  ASL_obj_replace_eq,
				    ASL_obj_replace_ineq | ASL_obj_replace_eq };

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

	if (!(mint_val[set_objno].U = no = n_obj))
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
		multiobj = 0;
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
#ifndef NO_MOkwf /*{*/
	if (no > 1 && multiobj && obj_no >= 0 && !Optimize) {
		Optimize = Optimizemo;
		method = 1;
		use_netopt = 0;
		algname = baralgname = "";
		}
	else
#endif /*}*/
	if (!Optimize || nqc) {	/* check nqc in case baropt was specified */
				/* with CPLEX 12.2, at least, we cannot use */
				/* CPXhybbaropt on problems with quadratic constraints. */
		Optimize = CPXlpopt;
		if (!lpoptalg)
			lpoptalg = CPX_ALG_AUTOMATIC;
#ifdef CPXERR_QCP_SENSE
		if (nqc | nlo) {
			Optimize = CPXqpopt;
			lpoptalg = CPX_ALG_AUTOMATIC;
			}
#endif
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

	b = LUrhs = (double *)M1alloc(2*sizeof(double)*(m+n));
	Urhsx = b + m;
	LUv = Urhsx + m;
	Uvx = LUv + n;
	L = nzr*(sizeof(double) + sizeof(int)) + (n+1)*sizeof(CStype);
	a = A_vals = (double *)Malloc(L);
	zap[0] = M1record(a);
	zap[1] = zap[2] = zap[3] = zap[4] = 0;
	ka = A_colstarts = (CStype*)(a + nzr);
	ia = A_rownos = (int*)(ka + n + 1);
	nelqf = 0;
#ifdef Want_mqpcheck /*{{*/
	v = 0;
	nint = nbv + niv + nlvbi + nlvci + nlvoi;
	if (Optimize == Optimize2 || Optimize == Optimizebar || nlo || nqc
	    || (nint && mipstval))
		want_xpi0 = 3;
	qp_read(*nl, ALLOW_CLP | repmap[objrep] ForceZ);
	nqc = nlc;	/* may have been adjusted */
	*nl = 0;	/* Indicate closed file -- in case we're */
			/* subsequently interrupted. */
	if (obj_no >= 0 && obj_no < no) {
		if (nqc) {
			og = og1 = Ograd[obj_no];
			for(i = 0; og; og = og->next)
				++i;
			d->ogsave = (ograd**)M1alloc(
				no*(sizeof(ograd*) + sizeof(cde)) + i*sizeof(ograd));
			memcpy(d->ogsave, Ograd, no*sizeof(ograd*));
			d->objsave = (cde*)(d->ogsave + no);
			memcpy(d->objsave, ((ASL_fg*)asl)->I.obj_de_, no*sizeof(cde));
			ogp = &d->ogsave[obj_no];
			og = (ograd*)(d->objsave + no);
			while(og1) {
				og->varno = og1->varno;
				og->coef = og1->coef;
				*ogp = og;
				ogp = &og->next;
				og1 = og1->next;
				++og;
				}
			*ogp = 0;
			}
		if (nlo && (nelqf = nq = mqpcheckv(obj_no, &qpi, &v))) {
			if (nq < 0) {
				nonlin(nq == -2, 555,
				 "a quadratic objective involving division by 0");
				nonlin(1, 551, "a nonlinear objective");
				}
#ifdef BARRIER_FOR_QP
			nonlin(nint, 553,
				"integer variables and a quadratic objective");
#endif
			crossover = 0;
#ifndef BARRIER_FOR_QP
#ifdef REQCONVEX
			if (nint) {
				i = 0;
				CPXgetintparam(Env, REQCONVEX, &i);
				if (i == 2) {
					printf("Changing reqconvex from 2 to 3 for a "
						"QP with integer variables.\n");
					CPXsetintparam(Env, REQCONVEX, 3);
					}
				}
#endif
			if (Optimize == Optimize2)
#endif
				set_baropt();
			method = 1;
			}
		}
#ifdef BARRIER_FOR_QP
	else
		nonlin(nlvoi, 552, "nonlinear integer variables");
#endif
	if (mipbasis == -1)
		mipbasis = endbas ? 1 : nlc ? want_qcdual : nlo ? 0 : 1;
#ifdef CPXERR_QCP_SENSE /* if CPLEX version >= 9.0 */

	iaq = kaq = kalq = 0; /* silence buggy warning */
	nelqc = 0;		/* ditto */
	colqc = rowqc = 0;	/* ditto */
	aq = 0; /* ditto */
	qcmat = 0; /* ditto */
	if (nqc > 0)
		zap[0] = linadj(asl, ka, ia, a, &kaq, &kalq, &iaq, &aq,
				&nelqc, &rowqc, &colqc, &qcmat, d, &v);

#endif /*CPXERR_QCP_SENSE*/
#else /*}{Want_mqpcheck*/
	nonlin(nlvoi, 552, "nonlinear integer variables");
	f_read(*nl, 0 ForceZ);
	*nl = 0;
#endif/*}}Want_mqpcheck*/

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
	vmi = 0;
	if (obj_no >= 0 && obj_no < no) {
		og = Ograd[obj_no];
		if (asl->i.vmap)
			vmi = get_vminv_ASL(asl);
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
		a = (double *)Malloc((nz1+2*n+m)*sizeof(double)
					+ nz1*sizeof(int) + (n+1)*sizeof(CStype));
		*zap[0] = a;
		l = a + nz1;
		u = l + n;
		b = u + n;
		ka = (CStype*)(b + m);
		ia = (int *)(ka + n + 1);
		nr = n;
		}
	else {
		m += objadj;
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

	rngvec = (double *)Malloc(m*sizeof(double) + nr*sizeof(int));
	zap[1] = M1record((void*)rngvec);
	kal = (int *)(rngvec + m);
#ifdef CPLEX_MIP
#ifdef CPXPARAM_Benders_Strategy
	if (Optimize == CPXbendersopt) {
		baralgname = "";
		if (nint <= 0 || nint >= n0) {
			fprintf(Stderr, "Ignoring bendersopt since all variables are %s.",
				nint <= 0 ? "continuous" : "integer");
			if (nint1)
				Optimize = CPXmipopt;
			else {
				Optimize = CPXlpopt;
				if (!lpoptalg)
					lpoptalg = CPX_ALG_AUTOMATIC;
				}
			}
		}
	else
#endif
#ifndef NO_MOkwf
	if (Optimize == Optimizemo)
		MOcheck(asl);
	else
#endif
	if (nint1) {
		Optimize = CPXmipopt;
		baralgname = "";
		}
#endif
	if (method > 0) {
		b = (double*)Malloc(m*sizeof(double));
		zap[2] = M1record((void*)b);
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
		if (vmi)
			for(; og; og = og->next)
				b[vmi[og->varno]] = os * og->coef;
		else
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
		if (vmi)
			for(; og; og = og->next)
				c[vmi[og->varno]] = og->coef;
		else
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
#ifdef CPX_PARAM_XXXIND
	 || parval(CPX_PARAM_XXXIND)
#endif
		) {
		make_names(&cname, &rname);
		zap[3] = M1record((void*)cname);
		}
	mq = m - nqc;
#ifdef Uselazy
	i = mq;
	zap[4] = lazyadj(asl, &LI, n, nint, nqc, &mq, b+nqc, senx+nqc, ka, kal,
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
#if CPX_VERSION >=  12060200
	if (obj_adj && (i = CPXchgobjoffset(Env, cpx, obj_adj)))
		badret("CPXchgobjoffset", i, 531);
#endif

#ifdef CPXERR_QCP_SENSE /* if CPLEX version >= 9.0 */
	lrhs = LUrhs;
	urhs = Urhsx;
	for(i = 0; i < nqc; i++) {
		j = CPXXaddqconstr(Env, cpx, kalq[i], nelqc[i], b[i], senx[i],
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
#ifdef CPXPARAM_Benders_Strategy
		if (Optimize == CPXbendersopt)
			benders_sufcheck(asl, cpx);
#endif
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
		colq = (CPXNNZ*)Malloc(L = (sizeof(int)+sizeof(CPXNNZ))*nr);
		colqcnt = (int*)(colq + nr);
		rowq = qpi->rowno;
		colqf = qpi->colbeg;
		qmat = qpi->delsq;
		colno = qpi->colno;
		k = qpi->nc;
		for(i = j0 = k0 = 0; i < k; i++) {
			j = colno[i];
			while(j0 < j) {
				colqcnt[j0] = 0;
				colq[j0++] = k0;
				}
			colq[j] = colqf[i];
			colqcnt[j] = (int)((k0 = colqf[i+1]) - colqf[i]);
			j0 = j + 1;
			}
		while(j0 < nr) {
			colqcnt[j0] = 0;
			colq[j0++] = k0;
			}
		k = n_var;
		for(i = 0; ; i++) {
			if (i == k) {
				i = qmatadj(k, nr, objsen, colq, colqcnt, &qmat, ctype);
				if (!i) {
					for(i = k; i < nr; i++)
						qmat[i++] = 0.;
					i = CPXcopyqpsep(Env,cpx,qmat);
					if (i)
						badret("CPXcopyqpsep",i,531);
					}
				free(qpi);
				if (i)
					exit(1);
				baralgname = "separable QP ";
				break;
				}
			if (!(j = colqcnt[i]))
				continue;
			if (j > 1 || rowq[colq[i]] != i) {
				i = CPXXcopyquad(Env,cpx,colq,colqcnt,rowq,qmat);
				free(qpi);
				if (i)
					badret("CPXcopyquad", i, 531);
				baralgname = "QP ";
				break;
				}
			}
		}
	if (v)
		mqpcheckv_free(&v);
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
	for(i = 5; i > 0; )
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
#if CPX_VERSION >= 1100
	if (droptol > 0.)
		CPXcleanup(Env, cpx, droptol);
#endif
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
	if (stat > 15 || stat <= 0) {
#ifndef CPXERR_ALGNOTLICENSED
#define CPXERR_ALGNOTLICENSED 32024
#endif
		switch(stat) {
		  case CPXERR_ALGNOTLICENSED:
			return (obj_no >= 0 && obj_no < nlo) || nlc ? 34 : 32;
#ifdef CPXMIP_POPULATESOL_LIM
		  case CPXMIP_POPULATESOL_LIM:
		  case CPXMIP_OPTIMAL_POPULATED:
		  case CPXMIP_OPTIMAL_POPULATED_TOL:
			stat = 36;
			break;
#endif
#ifdef CPX_STAT_FIRSTORDER
		  case CPX_STAT_FIRSTORDER:
			stat = 35;
			break;
#endif
		  case CPXMIP_DETTIME_LIM_FEAS:
			stat = 63;
			break;
		  case CPXMIP_DETTIME_LIM_INFEAS:
			stat = 64;
			break;
#ifdef CPX_PARAM_FEASOPTMODE /* >= 9.2b {*/
		  case CPX_STAT_FEASIBLE_RELAXED_SUM:
			stat = 65; break;
		  case CPX_STAT_OPTIMAL_RELAXED_SUM:
			stat = 66; break;
		  case CPX_STAT_FEASIBLE_RELAXED_INF:
			stat = 67; break;
		  case CPX_STAT_OPTIMAL_RELAXED_INF:
			stat = 68; break;
		  case CPX_STAT_FEASIBLE_RELAXED_QUAD:
			stat = 69; break;
		  case CPX_STAT_OPTIMAL_RELAXED_QUAD:
			stat = 70; break;
		  case CPXMIP_FEASIBLE_RELAXED_SUM:
			stat = 71; break;
		  case CPXMIP_OPTIMAL_RELAXED_SUM:
			stat = 72; break;
		  case CPXMIP_FEASIBLE_RELAXED_INF:
			stat = 73; break;
		  case CPXMIP_OPTIMAL_RELAXED_INF:
			stat = 74; break;
		  case CPXMIP_FEASIBLE_RELAXED_QUAD:
			stat = 75; break;
		  case CPXMIP_OPTIMAL_RELAXED_QUAD:
			stat = 76; break;
#ifdef CPX_STAT_MULTIOBJ_OPTIMAL /*{*/
		  case CPX_STAT_MULTIOBJ_OPTIMAL:
			stat = 77; break;
		  case CPX_STAT_MULTIOBJ_INFEASIBLE:
			stat = 78; break;
		  case CPX_STAT_MULTIOBJ_INForUNBD:
			stat = 79; break;
		  case CPX_STAT_MULTIOBJ_UNBOUNDED:
			stat = 80; break;
		  case CPX_STAT_MULTIOBJ_NON_OPTIMAL:
			stat = 81; break;
		  case CPX_STAT_MULTIOBJ_STOPPED:
			stat = 82; break;
#endif /*}*/
#endif /*}*/
		  default:
			stat = 0;
		  }
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

 static void
retfrom(PBuf *B, int k, const char *what)
{ Bpf(B, "\nReturn %d from %s.", k, what); }

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

 static void
send_confiis_ext(ASL *asl, cpxlp *cpx, PBuf *B)
{
	char *grptype;
	double *L, *U, *gp;
	int *ci, *grpbeg, *grpind, *gs, i, i1, j, j1, k, lc0, n, nc, nv, *vi;
	static int cmap[8]  = {0, 5, 6, 7, 4, 1, 3, 8 };

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
	gp = (double*)Malloc(n*(sizeof(double)+3*sizeof(int)+1)+2*k);
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
	lc0 = nlc;
	for(i = 0, k += n_con - lc0; j < k; j++) {
		grpind[j] = i++;
		grptype[j] = CPX_CON_LINEAR;
		}
	for(i = 0, k += nlc; j < k; j++) {
		grpind[j] = i++;
		grptype[j] = CPX_CON_QUADRATIC;
		}
	if ((i = CPXrefineconflictext(Env, cpx, n, n, gp, grpbeg, grpind, grptype))) {
		retfrom(B, i, "CPXrefineconflictext");
		goto ret;
		}
	if ((i = CPXgetconflictext(Env, cpx, gs, 0, n-1))) {
		retfrom(B, i, "CPXgetconflictext");
		goto ret;
		}
	for(i = 0; i < n; i++)
		if ((j = ++gs[i]) < 0 || j > 6)
			gs[i] = 7;
	vi = (int*)M1zapalloc((n_var+n_con)*sizeof(int));
	ci = vi + n_var;
	k = n_var;
	for(i = j = nv = 0; i < k; i++) {
		i1 = 0;
		if (L[i] > negInfinity) {
			if ((vi[i] = cmap[gs[j++]]))
				i1 = 1;
			}
		if (U[i] < Infinity && (j1 = cmap[gs[j++]])) {
			i1 = 1;
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
	Bpf(B, "\nReturning an IIS involving %d variables and %d constraints.", nv, nc);
 ret:
	free(gp);
	}

 static void
send_confiis(ASL *asl, cpxlp *cpx, PBuf *B)
{
	int cs, k, mc, mc1, nc, nc1;
	int *colind, *colst, *rowind, *rowst;
	static char *abort_reason[7] = {
		"contradiction", "time limit", "iteration limit", "node limit",
		"objective limit", "memory limit", "user request" };

	if ((k = CPXrefineconflict(Env, cpx, &mc, &nc))) {
		retfrom(B, k, "CPXrefineconflict");
		return;
		}

	k = mc + nc;
	if (mc < 0 || nc < 0 || k == 0)
		Bpf(B, "\nSurprise values mc = %d, nc = %d from CPXrefineconflict.",
			mc, nc);

	colind = (int*)Malloc(2*sizeof(int)*k);
	colst = colind + nc;
	rowind = colst + nc;
	rowst = rowind + mc;

	if ((k = CPXgetconflict(Env, cpx, &cs, rowind, rowst, &mc1, colind, colst, &nc1))) {
		retfrom(B, k, "CPXgetconflict");
		goto ret;
		}
	if (cs == CPX_STAT_CONFLICT_FEASIBLE) {
		Bpf(B, "\nNo IIS after all: problem is feasible!");
		goto ret;
		}
	if (cs != CPX_STAT_CONFLICT_MINIMAL) {
		if ((cs - CPX_STAT_CONFLICT_ABORT_CONTRADICTION) < 0
		 || cs > CPX_STAT_CONFLICT_ABORT_CONTRADICTION + 6)
			Bpf(B, "\nSurprise conflict status = %d from CPXgetconflict\n", cs);
		else
			Bpf(B, "\nSearch for conflicts aborted because of %s\n",
				abort_reason[cs - CPX_STAT_CONFLICT_ABORT_CONTRADICTION]);
		goto ret;
		}
	if (mc1 > mc || nc1 > nc) {
		Bpf(B, "\nSurprise confnumrows = %d (should be <= %d),"
				"\nconfnumcols = %d (should be <= %d) from CPXgetconflict.",
				mc1, mc, nc1, nc);
		goto ret;
		}
	nc = confiis_put(asl, ASL_Sufkind_var, nc1, colind, colst);
	mc = confiis_put(asl, ASL_Sufkind_con, mc1, rowind, rowst);
	Bpf(B, "\nReturning an IIS of %d variables and %d constraints.", nc, mc);

 ret:
	free(colind);
	}
#endif /*CPX_PARAM_FEASOPTMODE*/

 static void
send_iis(ASL *asl, cpxlp *cpx, PBuf *B)
{
	char *ckind;
	int i, j, m, m1, n, n1;
	int *ci, *cs, *ri, *rs;

	if (!want_iis)
		return;
#ifdef CPX_PARAM_FEASOPTMODE /* >= 9.2b */
	if (nlc) {
		send_confiis_ext(asl, cpx, B);
		return;
		}
	if (want_iis > 2
	|| (nbv + niv + nlvoi + nlvci + nlvbi > 0 && !relax)
	|| (Optimize != CPXprimopt && Optimize != CPXdualopt)) {
		send_confiis(asl, cpx, B);
		return;
		}
#else
	if (nlc) {
		printf("Ignoring iisfind request because of nonlinearities.\n");
		return;
		}
#endif
#ifdef CPX_PARAM_IISIND
	CPXsetintparam(Env, CPX_PARAM_IISIND, want_iis-1);
	if (i = CPXfindiis(Env, cpx, &m, &n)) {
		retfrom(B, i, "CPXfindiis");
		return;
		}
#else
	if ((i = CPXrefineconflict(Env, cpx, &m, &n))) {
		retfrom(B, i, "CPXrefineconflict");
		return;
		}
#endif
	ci = (int*)M1alloc((2*(m+n))*sizeof(int));
	cs = ci + n;
	ri = cs + n;
	rs = ri + m;
	ckind = "";
#ifdef CPX_PARAM_IISIND
	if (i = CPXgetiis(Env, cpx, &j, ri, rs, &m1, ci, cs, &n1)) {
		retfrom(B, i, "CPXgetiis");
		return;
		}
	iis_put(asl, ASL_Sufkind_var, n, ci, cs);
	iis_put(asl, ASL_Sufkind_con, m, ri, rs);
	if (j != 1)
		ckind = "partial ";
#else
	if ((i = CPXgetconflict(Env, cpx, &j, ri, rs, &m1, ci, cs, &n1))) {
		retfrom(B, i, "CPXgetconflict");
		return;
		}
	switch(j) {
	  case CPX_STAT_CONFLICT_MINIMAL:
		confiis_put(asl, ASL_Sufkind_var, n, ci, cs);
		confiis_put(asl, ASL_Sufkind_con, m, ri, rs);
		break;
	  case CPX_STAT_CONFLICT_FEASIBLE:
		Bpf(B, "\nCPXgetconflict claims the problem is feasible.");
		return;
	  default:
		Bpf(B, "\nSurprise conflict status %d from CPXgetconflict.", j);
		return;
	  }
#endif
	Bpf(B, "\nReturning %siis of %d variables and %d constraints.", ckind, n, m);
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
poolwrite(ASL *asl, cpxlp *cpx, dims *d, int nelqf, PBuf *B)
{
	PBuf B1;
	char buf[32], buf1[1200], *fname, *s;
	const char *whence;
	int j, k, m, n, nm1, npm, nprt, nsols, oprt, rft, srn, stat1, status;
	int *cstat, *ocstat, *orstat, *rstat;
	real feastol, inttol, obj, *x, *y, *y1;
	size_t L;
	static char
		errfmt1[] = "\nCPLEX solution status %d with fixed integers "
			"on solution pool\nmember %d (file %s):\n\t%s",
		errfmt2[] = "\nCPLEX status %d from CPXsolution "
			"on solution pool\nmember %d (file %s).";

	B1.s = buf1;
	B1.se = buf1 + sizeof(buf1);
	srn = solve_result_num;
	solve_result_num = 0;
	j = rft = 0;
	y = 0;
	ocstat = orstat = 0;
	if (populate == 1) {
		if ((k = CPXpopulate(Env, cpx))) {
			Bpf(B, "\nReturn %d from CPXpopulate.", k);
			goto done;
			}
		}
	nsols = CPXgetsolnpoolnumsolns(Env,cpx);
	if (nsols <= 0) {
		Bpf(B, "\n%d solutions in solution pool.", nsols);
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
		Bpf(&B1, "Solution pool member %d (of %d); objective %s",
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
				Bpf(&B1, errfmt1, stat1, j, fname,
					solinfo[statadjust(cpx,stat1)].msg);
				y1 = 0;
				}
			else if (CPXsolution(Env, cpx, &stat1, &obj, x, y, 0, 0)) {
				y1 = 0;
				status = CPXgetstat(Env, cpx);
				stat1 = statadjust(cpx, status);
				Bpf(&B1, errfmt2, status, j, fname);
				}
			else if (send_statuses(asl, cpx, d))
				Bpf(&B1, "\nNo basis.");
			if ((k = CPXchgprobtype(Env, cpx, oprt))) {
				whence = "CPXchgprobtype<restore>";
				break;
				}
			}
#ifdef CPXERR_QCP_SENSE
		if (nlc && y1 && x) {
			solve_result_num = 0;
			qcduals(asl, B, cpx, d, x, y1);
			}
#endif
		if (write_solf_ASL(asl, buf1, x, y1, 0, fname)) {
			whence = "";
			break;
			}
		B1.s = buf1;
		}
	Bpf(B, "\nWrote %d solution%s in solution pool", j, j == 1 ? "" : "s");
	s = L > 32 ? "\n" : " ";
	if (j > 0)
	  switch(j) {
		case 1:
			Bpf(B, "\nto file %s1.sol", poolstub);
			break;
		case 2:
			Bpf(B, "\nto files %s1.sol%sand %s2.sol",
				poolstub, s, poolstub);
			break;
		default:
			Bpf(B, "\nto files %s1.sol%s... %s%d.sol",
				poolstub, s, poolstub, j);
		}
	Bpf(B, ".");
	if (whence) {
		if (k)
			Bpf(B, "\nSolution pool writing stopped by return %d from %s.",
				k, whence);
		else
			Bpf(B, "\nCould not open \"%s\"", fname);
		}
	if (ocstat) {
		memcpy(cstat, ocstat, n*sizeof(int));
		memcpy(rstat, orstat, m*sizeof(int));
		}
	free(x);
 done:
	solve_result_num = srn;
	if (rft) /*restore "feasibility" tolerance*/
		CPXsetdblparam(Env, CPX_PARAM_EPRHS, feastol);
	d->npool = j;
	suf_iput("npool", ASL_Sufkind_obj, &d->npool);
	suf_iput("npool", ASL_Sufkind_prob, &d->npool);
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
tunewrite(char *fname, const char *fkeywd, PBuf *B)
{
	Element *e, *ee;
	FILE *f;
	ParInfo *PI;
	int n;
	keyword *kw;

	if (!(f = fopen(fname, "w"))) {
		Bpf(B, "\nCould not open %s file \"%s\"", fkeywd, fname);;
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
	Bpf(B, "\n%d settings written to %s file \"%s\"", n, fkeywd, fname);
	AVL_Tree_free(&PI->tree);
	free(PI);
	}
#endif /*}*/

#if CPX_VERSION >= 1000 /*{*/
 static void
tunewriteprm(char *fname, const char *fkeywd, PBuf *B)
{
	char *fmt;

	if (CPXwriteparam(Env, fname)) {
		fmt = "\nError writing %s file \"%s\"";
		solve_result_num = 565;
		}
	else
		fmt = "\n%s file \"%s\" written";
	Bpf(B, fmt, fkeywd, fname);
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
		badretfmt(564, "Could not open tunefixfile \"%s\"", tunefixfile);
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
tunerun(cpxlp *cpx, PBuf *B)
{
	Element *E;
	Tunefix TF;
	const char *what;
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
	if (k)
		Bpf(B, "\nSurprise return %d from CPXtuneparam.", k);
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
			Bpf(B, "\nTuning got surprise tunestat = %d.", j);
			goto no_what;
		 }
		Bpf(B, "\nTuning %s.", what);
		}
 no_what:
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

#ifdef CPXPROB_MIQCP
 static void
fixints(ASL *asl, cpxlp *cpx, real *x)
{
	char *t;
	int i, k, nint, *z, *z1;
	real *b, *b1;

	nint = nbv + niv + nlvoi + nlvci + nlvbi;
	b = b1 = (real*)Malloc(nint*(sizeof(real) + sizeof(int) + 1));
	z = z1 = (int*)(b + nint);
	t = (char*)(z + nint);
	memset(t, 'B', nint);
	i = n_var - nbv - niv;
	if ((k = nbv)) {
		memcpy(b1, x + i, k*sizeof(real));
		b1 += k;
		while(k--)
			*z1++ = i++;
		}
	if ((k = niv)) {
		memcpy(b1, x + i, k*sizeof(real));
		b1 += k;
		while(k--)
			*z1++ = i++;
		}
	if ((k = nlvbi)) {
		i = nlvb - k;
		memcpy(b1, x + i, k*sizeof(real));
		b1 += k;
		while(k--)
			*z1++ = i++;
		}
	if ((k = nlvci)) {
		i = nlvc - k;
		memcpy(b1, x + i, k*sizeof(real));
		b1 += k;
		while(k--)
			*z1++ = i++;
		}
	if ((k = nlvoi)) {
		i = nlvo - k;
		memcpy(b1, x + i, k*sizeof(real));
		b1 += k;
		while(k--)
			*z1++ = i++;
		}
	i = CPXchgbds(Env, cpx, nint, z, t, b);
	if (i)
		fprintf(Stderr, "CPxchbds() returned %d\n", i);
	free(b);
	}
#endif

 static void
amplout(ASL *asl, cpxlp **pcpx, dims *d, int status, int nelqf, int nint1, int *nosp)
{
	Optalg Contopt;
	PBuf B;
	Sol_info *SI;
	char buf[32], hbuf[4096];
	const char *wb, *what;
	cpxlp *cpx, *cpx1;
	int bitc, i, ii, itc, itci, j, m, mipstat, n, needsol, nint;
	int nodecnt, nodelim, nos, npt, nround, opt, stat, stat0, stat1, stat10;
	real *bb, *bn, *l, *le, *u, w[4], *x, *x2, *y, *z, *z1;
	real absmipgap, bbound, bcond, bobj, bobj0, feastol, inttol;
	real relmipgap, obj, t;
#ifndef NO_MOkwf
	int k, maxobjlen, objlen;
#endif
	static Sol_info solinfo1[] = {
	 { "QP Hessian is not positive semi-definite", 542, 0 },
	 { "primal optimal (no dual solution available)", 004, 1 },
	 { "ran out of memory", 522, 0 }
#ifdef CPXERR_Q_NEG_ZERO_COMP	/* CPLEX 7.1 */
	 ,{ "QP Hessian has diag. elements of the wrong sign", 541, 0 }
#endif
		};
#ifdef CPXERR_NOT_MIP
	static Sol_info solinfo2[1] = {
	 { "problem is not a MIP or has a nonconvex quadratic constraint", 543, 0 }};
#endif
	static Sol_info subprob_failed = { "Failed to solve a MIP subproblem", 513, 0 };
	static char *netmsg[] = {
		"no network",
		/*"an optimal network solution"*/ "",
		"an unbounded network",
		"an infeasible network",
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

	B.s = hbuf;
	B.se = hbuf + sizeof(hbuf);
	cpx = *pcpx;
	absmipgap = obj = relmipgap = Infinity;
	stat0 = stat10 = -1;
	stat1 = 1;
	ii = mipstat = nodecnt = npt = nround = opt = 0;
	m = n_con;
	n = n_var;
	x = d->x;
	y = d->y;
#ifdef USE_CHANNELS
	if (logfname)
		CPXsetlogfilename(Env, NULL, NULL);
#else
	if (Logf)
		CPXsetlogfile(Env, NULL); /* suppress idiotic Error 1217 */
#endif
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
	nos = nosp ? *nosp : 0;
	switch (status) {
	  case 0: break;
	  case CPXERR_SUBPROB_SOLVE:
		SI = &subprob_failed;
		x = y = 0;
		Bpf(&B, "%s: %s", Oinfo.bsname, SI->msg);
		if ((stat0 = CPXgetsubstat(Env,cpx)) && (stat = statadjust(cpx,stat0)))
			Bpf(&B, ": %s", solinfo[stat].msg);
		else
			Bpf(&B, "; CPXgetsubstat returned %d.", stat0);
		goto have_SI1;
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
#ifdef CPXERR_NOT_MIP
	  case CPXERR_NOT_MIP:
		SI = solinfo2;
		goto have_SI;
#endif
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
			  case CPX_ALG_MIP:
				/* possible with reqconvex == 3 */
				nint1 = 1;
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
			if (mipbasis) {
				opt = CPXgetprobtype(Env, cpx);
				Contopt = CPXprimopt;
				what = "CPXprimopt";
#ifdef CPXPROB_MIQCP
				npt = CPXPROB_FIXEDMILP;
				switch(opt) {
				  case CPXPROB_MIQP:
					npt = CPXPROB_FIXEDMIQP;
					goto use_qpopt;
				  case CPXPROB_MIQCP:
					if (!x)
						goto no_qcbasis;
					npt = CPXPROB_QCP;
 use_qpopt:
					Contopt = CPXqpopt;
					what = "CPXqpopt";
					}
#else
				npt = nelqf ? 8 : 3;
#endif
				if (npt != opt) {
					i = CPXchgprobtype(Env, cpx, npt);
					if (i) {
						fprintf(Stderr,
						 "CPXchgprobtype(...,%d) returned %d\n",
							npt, i);
						goto no_qcbasis;
						}
					}
#ifdef CPXPROB_MIQCP
				if (npt == CPXPROB_QCP)
					fixints(asl, cpx, x);
#endif
				if ((i = !CPXgetdblparam(Env,CPX_PARAM_EPRHS,&feastol)
				 && !CPXgetdblparam(Env,CPX_PARAM_EPINT,&inttol)
				 && inttol > feastol))
					CPXsetdblparam(Env, CPX_PARAM_EPRHS, inttol);
				status = Contopt(Env, cpx);
				stat1 = CPXgetstat(Env, cpx);
				if (i)	/*restore "feasibility" tolerance*/
					CPXsetdblparam(Env, CPX_PARAM_EPRHS, feastol);
				if (status) {
					fprintf(Stderr, "\n%s returned %d\n", what, status);
					y = 0;
					}
				else if (stat1 != 1) {
					fprintf(Stderr,
					"CPLEX solution status %d with fixed integers:\n\t%s\n",
						stat1, solinfo[statadjust(cpx,stat1)].msg);
					y = 0;
					}
				else if (x)
					needsol = 1;
				stat1 = 1;
				if (endbas && !nelqf)
					write_basis(cpx);
				}
			else if (!nlc)
				y = 0;
 no_qcbasis:
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
	Bpf(&B, "%s: %s", Oinfo.bsname, SI->msg);
	if (!stat)
		Bpf(&B, ": %s", failstat(stat0));
#ifdef CPLEX_MIP
	else if (stat == 16 || stat == 17) {
		CPXgetintparam(Env, CPX_PARAM_NODELIM, &nodelim);
		Bpf(&B, ".\nCurrent node limit = %d", nodelim);
		}
#endif
	else if (!SI->wantobj)
		Bpf(&B, ".");

	if (SI->wantobj) {
#ifndef NO_MOkwf
		if (Optimize == Optimizemo) {
			maxobjlen = 0;
			for(i = 0; i < nmo; ++i) {
				objlen = strlen(obj_name(indobj[i]));
				if (maxobjlen < objlen)
					maxobjlen = objlen;
				}
			++maxobjlen;
			Bpf(&B, "; %d objectives:\n", nmo);
			for(i = 0; i < nmo; ++i) {
				k = objind[i];
				if ((j = CPXmultiobjgetobjval(Env, cpx, k, &obj)))
					badret("CPXmultiobjgetobjval", j, 1);
				g_fmtop(buf, obj);
				Bpf(&B, "\t%-*s= %s\n", maxobjlen, obj_name(indobj[k]), buf);
				}
			Bpf(&B, "Blended objective value%s:\n", nmopri == 1 ? "" : "s");
			for(i = nmopri; i > 0; --i) {
				if ((j = CPXmultiobjgetobjvalbypriority(Env, cpx, i, &obj)))
					badret("CPXmultiobjgetobjvalbypriority", j, 1);
				g_fmtop(buf, obj);
				Bpf(&B, "\tPriority %d (.pri = %d): %s\n", i, pri0[i-1], buf);
				}
			}
		else
#endif
			{
			if (obj == Infinity)
				CPXgetobjval(Env, cpx, &obj);
			g_fmtop(buf, obj);
			Bpf(&B, "; objective %s", buf);
			}
		}
 have_SI1:
	solve_result_num = SI->code;
	if (!SI->wantobj && SI->code >= 300 && SI->code < 400)
		x = y = 0;
	if (nosp) {
		if (nos < 0 || nos > 12)
			nos = 11;
		if (nos != 1)
			Bpf(&B, "\nnetopt found %s.", netmsg[nos]);
		if (net_nodes != lnc || net_arcs != nwv)
			Bpf(&B, "\nNetwork extractor found %d nodes and %d arcs.",
				net_nodes, net_arcs);
		Bpf(&B, "\n%d network simplex iterations.", netiters);
		}
	if (nint1 && !mipstat)
		Bpf(&B, "\n%d MIP simplex iterations\n%d branch-and-bound nodes",
			ii, nodecnt);
	if (itc > 0 || itci > 0
	 || Optimize == CPXprimopt || Optimize == CPXdualopt) {
		Bpf(&B, "\n%d %s%ssimplex iterations (%d in phase I)%s",
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
		Bpf(&B, "\n%d sifting subproblems solved (%d in phase I)", itc, itci);
		}
#endif
	if (bitc > 0) {
		Bpf(&B, "\n%d %sbarrier iterations", bitc, baralgname);
		if (cbi.nx[0] | cbi.nx[1])
			Bpf(&B, "\n%d push, %d exchange %s crossover iterations",
				cbi.nx[0], cbi.nx[1], cbi.xkind);
		}
	if (cbi.np[0] + cbi.np[1] + cbi.np[3])
		Bpf(&B, "\nCPLEX's %spresolve eliminated totals of %d "
			"constraints and %d variables\n(perhaps repeatedly) "
			"and made %d coefficient changes.",
			nint1 ? "MIP " : "", cbi.np[0], cbi.np[1], cbi.np[3]);
	if (cbi.np[2])
		Bpf(&B, "\nCPLEX's aggregator made %d substitutions.", cbi.np[2]);
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
		if (bobj >= 1e75)
			bobj = Infinity;
		else if (bobj <= -1e75)
			bobj = negInfinity;
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
				Bpf(&B, "\n%d of %d clique inequalities used", cl, cl0);
			if (cov)
				Bpf(&B, "\n%d cover cuts added", cov);
			}
		else
			Bpf(&B, "\nNo clique or cover cuts used.");
		}
#endif
	if (absmipgap > 0. && absmipgap < Infinity && !(retmipgap & 4))
			Bpf(&B, "\nabsmipgap = %g, relmipgap = %g",
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
		Bpf(&B, "\nStatus recovering solution: %s", SI->msg);
		if (!stat1)
			Bpf(&B, ": stat = %d", stat10);
		if (solve_result_num < SI->code)
			solve_result_num = SI->code;
		}
	if (aggtries > 1)
		Bpf(&B, "\nTried aggregator %d times", aggtries);
	if (solve_result_num >= 200 && solve_result_num < 300)
		send_iis(asl, cpx, &B);
	if (send_statuses(asl, cpx, d))
		Bpf(&B, "\nNo basis.");
	else if (asl->i.flags & 1 && method > 0) {
		j = 1;
		switch(stat) {
			/* Note: stat values have been adjusted to accord with */
			/* earlier CPLEX versions -- don't use CPX_STAT_INFEASIBLE */
			/* or CPX_STAT_UNBOUNDED here... */
		  case 2: /*infeasible*/
			if (wantray(1,&cpx,&itc,&itci,nelqf)) {
				j = send_dray(asl,cpx,nelqf);
				Bpf(&B, dray_msg[j]);
				}
			goto j_check;
		  case 3: /*unbounded*/
			if (wantray(0,&cpx,&itc,&itci,0))
			  Bpf(&B, (j = send_ray(asl,cpx))
				? "\nfailed to compute variable.unbdd"
				: "\nvariable.unbdd returned");
 j_check:
			*pcpx = cpx;
			if (itc) {
				Bpf(&B, "\n%d extra %ssimplex iterations for ray",
					itc, algname);
				if (itci)
				  Bpf(&B, " (%d in phase I)",itci);
				}
			if (!j)
				solve_result_num += 10;
		  }
		}
	if (basis_cond) {
		if (CPXgetdblquality(Env, cpx, &bcond, CPX_KAPPA)) {
			Bpf(&B, "\nBasis condition is unavailable.");
			bcond = 0;
			}
		else
			Bpf(&B, "\nBasis condition = %g", bcond);
		if (n_obj > 0)
			suf_rput("basis_cond", ASL_Sufkind_obj, &bcond);
		suf_rput("basis_cond", ASL_Sufkind_prob, &bcond);
		}
#ifdef CPX_PARAM_MIPKAPPASTATS
	if (haveqmet) {
		Bpf(&B, "\n%s = %.3g", QualInfo[0].desc, qmet[0]);
		for(j = 1; j < haveqmet; ++j)
			if (qmet[j] > 0.)
				Bpf(&B, "\n%s = %.3g", QualInfo[j].desc, qmet[j]);
		}
#endif
	if (nround) {
		wb = "";
		if (nround < 0) {
			nround = -nround;
			wb = "would be ";
			}
		Bpf(&B, "\n%d integer variables %srounded (maxerr = %g).",
			nround, wb, w[0]);
		if (w[0] > intwarn_tol && 0.5*w[0] > w[2])
			Bpf(&B, "\nAssigning integrality = %.1g might help."
			  "\nCurrently integrality = %g.", 0.5*w[0], w[1]);
		}
#if CPX_VERSION >= 1100
	if (npt != opt)
		CPXchgprobtype(Env, cpx, opt);
	if (cutstats && nint1)
	    for(CI = Cut_Info; CI->cutname; ++CI)
		if (!CPXgetnumcuts(Env, cpx, CI->cuttype, &j) && j > 0)
			Bpf(&B, "\n%d %s cut%s", j, CI->cutname, j == 1 ? "" : "s");
#endif
#ifdef CPX_PARAM_POPULATELIM
	if (Optimize == CPXmipopt && poolstub && x)
		poolwrite(asl, cpx, d, nelqf, &B);
	if (pretunefile)
		tunewrite(pretunefile, "pretunefile", &B);
#endif
#if CPX_VERSION >= 1000 /*{*/
	if (pretunefileprm)
		tunewriteprm(pretunefileprm, "pretunefileprm", &B);
#ifdef CPX_TUNE_TILIM /* CPLEX 11 */
	if (tunefile || tunefileprm) {
		tunerun(cpx, &B);
		if (tunefile)
			tunewrite(tunefile, "tunefile", &B);
		if (tunefileprm)
			tunewriteprm(tunefileprm, "tunefileprm", &B);
		}
#endif
#endif /*}*/
#ifdef CPXERR_QCP_SENSE
	if (nlc) {
		if (y && x && want_qcdual && SI->wantobj)
			qcduals(asl, &B, cpx, d, x, y);
		else
			y = 0;
		}
#endif
	write_sol(hbuf, x, y, &Oinfo);
#ifdef USE_CHANNELS
	if (logfname) {
		CPXsetlogfilename(Env, logfname, "a");
		CPXgetchannels(Env, NULL, NULL, NULL, &cpxlog);
		Bpf(&B, "\n");
		CPXmsgstr(cpxlog, hbuf);
		CPXsetlogfilename(Env, NULL, NULL);
		logfname = 0;
		}
#else
	if (Logf) {
		FILE *f;

		CPXsetlogfile(Env, NULL);
		CPXfclose(Logf);
		Logf = 0;
		if ((f = fopen(*file_name[set_logname], "a"))) {
			fprintf(f, "%s\n", hbuf);
			fclose(f);
			}
		}
#endif
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
#if CPX_VERSION >= 12070100
	if (!strncmp(msg, "CPXPARAM", 8))
		return;
#endif
#ifdef USE_CHANNELS
	if (logfname)
		CPXmsgstr(cpxlog, msg);
#else
	if (Logf)
		CPXfputs(msg, Logf);
#endif
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
#ifdef CPXPARAM_Benders_Strategy
	{ "benders", 0, ASL_Sufkind_var },
#endif
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
	{ "npool", 0, ASL_Sufkind_obj | ASL_Sufkind_outonly},
	{ "npool", 0, ASL_Sufkind_prob | ASL_Sufkind_outonly},
#endif
#ifndef NO_MOkwf
	{ "objabstol", 0, ASL_Sufkind_obj | ASL_Sufkind_real },
	{ "objpriority", 0, ASL_Sufkind_obj },
	{ "objreltol", 0, ASL_Sufkind_obj | ASL_Sufkind_real },
	{ "objweight", 0, ASL_Sufkind_obj | ASL_Sufkind_real },
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
#ifdef Want_Distmipopt
	const char *cname;
	char *s, *t;
	int (CPXPUBLIC *vmconf)(CPXENVptr, const char*);
	size_t L;

	if ((s = vmconfig) && *s) {
		if (*s == '@') {
			++s;
			cname = "CPXreadcopyvmconfig";
			vmconf = CPXreadcopyvmconfig;
			}
		else {
			cname = "CPXcopyvmconfig";
			vmconf = CPXcopyvmconfig;
			}
		if ((k = vmconf(e,s))) {
			solve_result_num = 580;
			asl->i.uinfo = t = (char*)M1alloc(L = strlen(s) + 80);
			Snprintf(t, L, "Return %d from %s(e, \"%s\").", k, cname, s);
			longjmp(Jb,1);
			}
		rc = CPXdistmipopt(e, c);
		CPXdelvmconfig(e);
		}
	else
#endif
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
		  case CPXMIP_INFEASIBLE:
			rc = feasopt(e, c);
		  }
#endif /*CPX_PARAM_FEASOPTMODE*/
	return rc;
	}

 int
main(int argc, char **argv)
{
	CPXCHANNELptr cpxresults;
	char	*s, *solmsg, *stub;
	dims	d;
	const char *fmt;
	int	nelqf, nint1, nms, nos, *nosp, rc, status, z;
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
			L = strlen(Oinfo.bsname);
			fmt = "%s";
			if ((s = (char*)asl->i.uinfo)) {
				L += strlen(s);
				fmt = "%s: %s";
				}
			solmsg = (char*)M1alloc(L+3);
			sprintf(solmsg, fmt, Oinfo.bsname, s);
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
	if (wrtfname && (status = CPXwriteprob(Env, cpx, wrtfname, NULL)))
		badwrite(581, "writeprob=\"%s\" failed with status %d.\n",
			wrtfname, status);
	if (startbas)
		basread(asl, cpx);
#ifdef BARRIER
#if !defined(NO_DEPRECATED) && CPX_VERSION < 12070000
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
	if (wrtmipstart) {
		if (!(nms = CPXgetnummipstarts(Env, cpx))) {
			if (!(nl = fopen(wrtmipstart, "w")))
				badwrite(582, "writemipstart=\"%s\" failed: could not create an\n"
				 "empty file to indidate that no MIP start is available.\n",
					wrtmipstart, status);
			fclose(nl);
			nl = 0;
			}
		else if ((status = CPXwritemipstarts(Env, cpx, wrtmipstart, 0, nms-1)))
			badwrite(582, "writemipstart=\"%s\" failed:\n"
				"\tCPXwritemipstarts() returned %d.\n", wrtmipstart, status);
		}
	if (nosolve) {
		badretfmt(600, "Not solved because of \"nosolve\".\n");
		exit(1);
		}
	Times[1] = xectim_();
#if CPX_VERSION >= 12050000
	if (!DTimes_failed)
		DTimes_failed = CPXgetdettime(Env, &DTimes[1]);
#endif
	nosp = 0;

	CPXdisconnectchannel(Env, cpxresults);
	CPXaddfuncdest(Env, cpxresults, 0, mymsgfunc);
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
