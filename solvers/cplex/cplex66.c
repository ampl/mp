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

#ifdef __cplusplus
#define CPX_PROTOTYPE_ANSI
#else
#ifdef KR_headers
#define CPX_PROTOTYPE_CKR
#define size_t int
#else
#define CPX_PROTOTYPE_ANSI
#endif
#endif

#include "cplex.h"
#define STDIO_H_included
#include "nlp.h"
#include "getstub.h"
#include "signal.h"

#ifdef Sig_ret_type
#ifndef SigRet
#define SigRet return 0
#endif
#else
#define Sig_ret_type void
#define SigRet /*nothing*/
#endif

#ifndef Sig_func_type
typedef void sig_func_type ANSI((int));
#endif

static ASL *asl;

typedef struct cpxlp cpxlp;

 typedef struct
mint_values {
	int L;
	int U;
	int val;
	} mint_values;

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
	set_fcpri	= 12,
	set_sos		= 13,
	set_mipcuts	= 14,
	set_round	= 15
	};

 static mint_values
mint_val[16] = {
	/* set_crossover */	0, 2, 1,
	/* set_dualthresh */	-0x7fffffff, 0x7fffffff, 32000,
	/* set_netopt */	0, 2, 1,
	/* set_objno */		0, 0/*n_obj*/,	1,
	/* set_conpri */	0, 0x7fffffff, 1,
	/* set_objpri */	0, 0x7fffffff, 2,
	/* set_prestats */	0, 1, 0,
	/* set_sos2 */		0, 1, 1,
	/* set_timing */	0, 3, 0,
	/* set_iis */		0, 2, 0,
	/* set_mipststat */	0, 1, 1,
	/* set_mipstval */	0, 1, 1,
	/* set_fcpri */		0, 0x7fffffff, 3,
	/* set_sos */		0, 1, 1,
	/* set_mipcuts */	-1, 2, 0,
	/* set_round */		0, 15, 1
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
#define fcpri		mint_val[12].val
#define sos		mint_val[13].val
#define Round		mint_val[15].val

 static int hybmethod = CPX_ALG_PRIMAL;
 static int netiters = -1;
 static CPXFILEptr Logf;
 extern char cplex_version[];
 static char *algname, *baralgname, *endbas, *endtree, *endvec, *logfname;
 static char *startbas, *starttree, *startvec, *wrtfname, *version;
 static int bestnode, breaking, costsens, dispzap, mbas, method;
 static int nbas, netopting, objadj, objsen, relax, zap_lpcbf;
 static int aggtries, net_status, net_nodes, net_arcs;
 static double obj_adj;
 typedef struct
Cbfinfo {
	int agg, disp, m, mipdisp, n, n0[4], nx[2], nint, npr, pres;
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
	} dims;

#ifdef DEBUG
 static void
Basedebug(char *name, int *cs, int *rs) /*!!!!*/
{
	FILE *f;
	int i;

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
#ifdef KR_headers
bascopy(cpx, d, n, m, cmap, ncs, nrs, rmap)
	CPXLPptr cpx; dims *d; int n, m, *cmap, *ncs, *nrs, *rmap;
#else
bascopy(CPXLPptr cpx, dims *d, int n, int m, int *cmap, int *ncs, int *nrs, int *rmap)
#endif
{
	int *cs, i, j, *rs;

	if (cs = d->cs)
		rs = d->rs;
	else {
		cs = (int*)Malloc((mbas+nbas)*sizeof(int));
		rs = cs + nbas;
		for (i = 0, j = nbas + mbas; i < j; i++)
			cs[i] = CPX_AT_LOWER;;
		}
	for(i = 0; i < n; i++)
		if ((j = cmap[i]) >= 0)
			cs[j] = ncs[i];
	for(i = 0; i < m; i++)
		if ((j = rmap[i]) >= 0)
			rs[j] = nrs[i];
	Basedebug("net.bas",cs,rs);
	CPXcopybase(Env, cpx, cs, rs);
	if (!d->cs)
		free(cs);
	}

 static int
#ifdef KR_headers
netopt(cpx, d, stat, nodes, arcs, its)
	CPXLPptr cpx; dims *d;
	int *stat, *nodes, *arcs, *its;
#else
netopt(CPXLPptr cpx, dims *d, int *st, int *nodes, int *arcs, int *its)
#endif
{
	CPXNETptr net;
	int i, j, k, m, n;
	int *cmap, *ncs, *nrs, *rmap, *s;

	k = 0;
	if (!(net = CPXNETcreateprob(Env,&k,"embedded")) || k)
		return 1;
	m = mbas + 1;
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
	if (s = d->cs) {
		n = *arcs;
		for(i = 0; i < n; i++)
			ncs[i] = (j = cmap[i]) >= 0 ? s[j] : CPX_BASIC;
		n = *nodes;
		s = d->rs;
		for(i = 0; i < n; i++)
			nrs[i] = (j = rmap[i]) >= 0 ? s[j] : CPX_BASIC;
		if (i = CPXNETcopybase(Env, net, ncs, nrs))
			fprintf(Stderr, "Return %d from CPXNETcopybase.\n", i);
		}
	if (i = CPXNETprimopt(Env, net)) {
		fprintf(Stderr, "Return %d from CPXNETprimopt\n", i);
		*st = 11;
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
#ifdef DEBUG
		CPXNETbasewrite(Env,net,"net.bw");
#endif
		if (!CPXNETgetbase(Env, net, ncs, nrs)) {
			bascopy(cpx, d, *arcs, *nodes, cmap, ncs, nrs, rmap);
#ifdef DEBUG
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

#ifdef KR_headers
 static int (CPXPUBLIC *Optimize)();
#else
 static int (CPXPUBLIC *Optimize)(CPXENVptr,  cpxlp*);
#endif

 static jmp_buf Jb;

 Sig_ret_type
#ifdef KR_headers
intcatch(n) int n;
#else
intcatch(int n)
#endif
{
	printf("\n<BREAK> (cplex)\n", n);
	fflush(stdout);
	if (++breaking > 3)
		longjmp(Jb, 2);
	signal(SIGINT, intcatch);
	SigRet;
	}

 static int CPXPUBLIC
#ifdef KR_headers
breakcallback(e, c, w, h) CPXENVptr e; void *c; int w; void *h;
#else
breakcallback(CPXENVptr e, void *c, int w, void *h)
#endif
{
	Not_Used(e);
	Not_Used(c);
	Not_Used(w);
	Not_Used(h);
	return breaking > 1;
	}

#ifdef BARRIER

 static int CPXPUBLIC
#ifdef KR_headers
Optimize2(e, c) CPXENVptr e; cpxlp *c;
#else
Optimize2(CPXENVptr e, cpxlp *c)
#endif
{
	return CPXhybbaropt(e, c, hybmethod);
	}

 static int CPXPUBLIC
Optimizebar(CPXENVptr e, cpxlp *c)
{
	int rv, s;

	if (!(rv = CPXbaropt(e, c)) && endvec
	 && (s = CPXvecwrite(e, c, endvec)))
		printf("\n*** return %d from CPXvecwrite.\n", s);
	return rv;
	}

 static void
set_baropt(VOID)
{
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
		algname = baralgname = "";
		}
	}

 static int
#ifdef KR_headers
startcomp(asl,n0,m,nextra,ka,kal,ia,a,b,c,cdualp,cprimp,rdualp,rprimp)
	ASL *asl;
	int n0, m, nextra, *ka, *kal, *ia;
	double *a, *b, *c, **cdualp, **cprimp, **rdualp, **rprimp;
#else
startcomp(ASL *asl, int n0, int m, int nextra, int *ka, int *kal, int *ia,
	double *a, double *b, double *c,
	double **cdualp, double **cprimp, double **rdualp, double **rprimp)
#endif
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
			if (t = x[i]) {
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

#ifdef MSDOS
#undef Stand_alone
#endif
#ifdef Student_Edition
#undef Stand_alone
#endif

 extern int cpxmain ANSI((int, char**));
typedef char *sfunc ANSI((char *, char *, int));

#ifdef CPLEX_MIP

 static void
#ifdef KR_headers
treeio(cpx, fname, what, tio) cpxlp *cpx; char *fname, *what;
	int (CPXPUBLIC *tio)();
#else
treeio(cpxlp *cpx, char *fname, char *what,
	int (CPXPUBLIC *tio)(CPXENVptr,cpxlp*,char*))
#endif
{
	if ((*tio)(Env,cpx,fname)) {
		printf("Could not %s tree file \"%s\".\n", what, fname);
		need_nl = 0;
		}
	}
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
	set_version	= 8,
	set_sens	= 9,
	set_bestnode	= 10
	};

 static char *
#ifdef KR_headers
sf_known(oi, kw, v) Option_Info *oi; keyword *kw; char *v;
#else
sf_known(Option_Info *oi, keyword *kw, char *v)
#endif
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

	 case set_version:
#ifdef Student_Edition
#define Fmt "%sStudent %s\n"
#else
#define Fmt "%s%s\n"
#endif
		printf(Fmt, oi->nnl ? "\n" : "", version);
		oi->option_echo &= ~ASL_OI_echothis;
#undef Fmt
		if (oi->option_echo & ASL_OI_clopt)
			exit(0);
		break;

	 case set_sens:
		costsens = 1;
		break;

	 case set_bestnode:
		bestnode = 1;
	 }
	return v;
	}

 static void
#ifdef KR_headers
badival(oi, kw, t, L, U) Option_Info *oi; keyword *kw; int t, L, U;
#else
badival(Option_Info *oi, keyword *kw, int t, int L, int U)
#endif
{
	printf("rejecting %s %d; must be between %d and %d\n",
		kw->name, t, L, U);
	badopt_ASL(oi);
	}

 static char *
#ifdef KR_headers
sf_mint(oi, kw, v) Option_Info *oi; keyword *kw; char *v;
#else
sf_mint(Option_Info *oi, keyword *kw, char *v)
#endif
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
		static int op[5] = {
			CPX_PARAM_CLIQUES,
			CPX_PARAM_COVERS,
			CPX_PARAM_FLOWCOVERS,
			CPX_PARAM_GUBCOVERS,
			CPX_PARAM_IMPLBD
			};
		int f;
		for(f = 0; f < 5; f++)
			CPXsetintparam(Env, op[f], t);
		}
#endif
	return rv;
	}

 static char *
#ifdef KR_headers
sf_int(oi, kw, v) Option_Info *oi; keyword *kw; char *v;
#else
sf_int(Option_Info *oi, keyword *kw, char *v)
#endif
{
	int f, t, z[3];
	char *rv, *what = kw->name;

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
	else if (CPXsetintparam(Env, f, t)) {
		z[2] = 0;
		z[1] = 1;
		CPXinfointparam(Env, f, z, z+1, z+2);
		badival(oi, kw, t, z[1], z[2]);
		}
	return rv;
	}

/* In case some reason surfaces to distinguish the old setzzzpar()	*/
/* calls from the old setzzzind calls (which had sfunc value sf_int and	*/
/* sf_int1, respectively, in the keywds table, we use the following	*/
/* #define and retain appearances of sf_int1 in the keywds table.	*/

#define sf_int1 sf_int

 static char *
#ifdef KR_headers
sf_dbl(oi, kw, v) Option_Info *oi; keyword *kw; char *v;
#else
sf_dbl(Option_Info *oi, keyword *kw, char *v)
#endif
{
	double t, z[3];
	char *rv, *what = kw->name;
	int f = Intcast kw->info;

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
	else if (CPXsetdblparam(Env, f, t)) {
		z[2] = 0;
		z[1] = 1;
		CPXinfodblparam(Env, f, z, z+1, z+2);
		printf("rejecting %s %g; must be between %g and %g\n",
			what, t, z[1], z[2]);
		badopt_ASL(oi);
		}
	return rv;
	}

 static char **file_name[8] = { &endbas, &endtree, &startbas, &starttree,
				&endvec, &startvec, &logfname, &wrtfname };

enum {	/* sf_char f values */
	set_endbas	= 0,
	set_endtree	= 1,
	set_startbas	= 2,
	set_starttree	= 3,
	set_endvector	= 4,
	set_startvector	= 5,
	set_logname	= 6,
	set_wrtfname	= 7
	};

 static char *
#ifdef KR_headers
sf_char(oi, kw, v) Option_Info *oi; keyword *kw; char *v;
#else
sf_char(Option_Info *oi, keyword *kw, char *v)
#endif
{
	char *rv, *t;
	int f, n;

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
	for(rv = v; *++rv > ' ';);
	t = M1alloc(n = rv - v + 1);
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

#define VP (Char*)

 struct option_word {
	char *name;
	sfunc *sf;
	int ival;
	};
 typedef struct option_word option_word;

 static keyword keywds[] = {	/* must be in alphabetical order */
#ifdef CPLEX_MIP
	{ "absmipgap",	sf_dbl,		VP CPX_PARAM_EPAGAP },
#endif
	{ "advance",	sf_int1,	VP CPX_PARAM_ADVIND },
	{ "aggfill",	sf_int,		VP CPX_PARAM_AGGFILL },
	{ "agglim",	sf_int,		VP CPX_PARAM_AGGFILL },
	{ "aggregate",	sf_int1,	VP CPX_PARAM_AGGIND },
#ifdef CPLEX_MIP
	{ "backtrack",	sf_dbl,		VP CPX_PARAM_BTTOL },
#endif
#ifdef BARRIER
	{ "baralg",	sf_int,		VP CPX_PARAM_BARALG },
	{ "barcorr",	sf_int,		VP CPX_PARAM_BARMAXCOR },
	{ "bardisplay",	sf_int,		VP CPX_PARAM_BARDISPLAY },
	{ "bargrowth",	sf_dbl,		VP CPX_PARAM_BARGROWTH },
	{ "bariterlim",	sf_int,		VP CPX_PARAM_BARITLIM },
	{ "barobjrange", sf_dbl,	VP CPX_PARAM_BAROBJRNG },
	{ "baropt",	sf_known,	VP set_barrier },
	{ "barstart",	sf_int,		VP CPX_PARAM_BARSTARTALG },
	{ "barthreads",	sf_int,		VP CPX_PARAM_BARTHREADS },
	{ "barvarup",	sf_dbl,		VP CPX_PARAM_BARVARUP },
#endif
	{ "basisinterval", sf_int,	VP CPX_PARAM_BASINTERVAL },
#ifdef CPLEX_MIP
	{ "bbinterval",	sf_int,		VP CPX_PARAM_BBINTERVAL },
	{ "bestnode",	sf_known,	VP set_bestnode },
	{ "boundstr",	sf_int,		VP CPX_PARAM_BNDSTRENIND },
	{ "branch",	sf_int,		VP CPX_PARAM_BRDIR },
	{ "branchdir",	sf_int,		VP CPX_PARAM_BRDIR },
	{ "cliquecuts",	sf_int,		VP CPX_PARAM_CLIQUES },
	{ "cliques",	sf_int,		VP CPX_PARAM_CLIQUES },
#endif
	{ "clocktype",	sf_int,		VP CPX_PARAM_CLOCKTYPE },
#ifdef CPLEX_MIP
	{ "coeffreduce", sf_int1,	VP CPX_PARAM_COEREDIND },
#endif
#ifdef BARRIER
	{ "comptol",	sf_dbl,		VP CPX_PARAM_BAREPCOMP },
#endif
#ifdef CPLEX_MIP
	{ "covercuts",	sf_int,		VP CPX_PARAM_COVERS },
	{ "covers",	sf_int,		VP CPX_PARAM_COVERS },
#endif
	{ "crash",	sf_int,		VP CPX_PARAM_CRAIND },
#ifdef BARRIER
	{ "crossover",	sf_mint,	VP set_crossover },
#endif
#ifdef CPLEX_MIP
	{ "cutsfactor",	sf_dbl,		VP CPX_PARAM_CUTSFACTOR },
#endif
#ifdef BARRIER
	{ "dense",	sf_int,		VP CPX_PARAM_BARCOLNZ },
	{ "densecol",	sf_int,		VP CPX_PARAM_BARCOLNZ },
#endif
	{ "dependency",	sf_int1,	VP CPX_PARAM_DEPIND },
	{ "dgradient",	sf_int,		VP CPX_PARAM_DPRIIND },
	{ "display",	sf_int,		VP CPX_PARAM_SIMDISPLAY },
	{ "doperturb",	sf_int1,	VP CPX_PARAM_PERIND },
	{ "dual",	sf_known,	VP set_dual },
	{ "dualopt",	sf_known,	VP set_dualopt },
	{ "dualthresh",	sf_mint,	VP set_dualthresh },
	{ "endbasis",	sf_char,	VP set_endbas },
#ifdef CPLEX_MIP
	{ "endtree",	sf_char,	VP set_endtree },
#endif
#ifdef BARRIER
	{ "endvector", 	sf_char,	VP set_endvector },
#endif
#ifdef CPLEX_MIP
	{ "fcvarpri",	sf_mint,	VP set_fcpri },
#endif
	{ "feasibility", sf_dbl,	VP CPX_PARAM_EPRHS },
	{ "file",	sf_char,	VP set_wrtfname },
#ifdef CPLEX_MIP
	{ "flowcuts",	sf_int,		VP CPX_PARAM_FLOWCOVERS },
#ifndef NO_CPLEX66 /* for versions prior to CPLEX 6.6 */
	{ "fraccuts",   sf_int,	        VP CPX_PARAM_FRACCUTS },
	{ "fractionalcuts", sf_int,	VP CPX_PARAM_FRACCUTS },
#endif
#endif
#ifdef BARRIER
	{ "growth",	sf_dbl,		VP CPX_PARAM_BARGROWTH }, /*== bargrowth*/
#endif
#ifdef CPLEX_MIP
	{ "gubcuts",	sf_int,		VP CPX_PARAM_GUBCOVERS },
	{ "heuristic",	sf_int,		VP CPX_PARAM_HEURISTIC },
	{ "heuristicfreq", sf_int,	VP CPX_PARAM_HEURFREQ },
#endif
	{ "iisfind",	sf_mint,	VP set_iis },
#ifdef CPLEX_MIP
	{ "impliedcuts", sf_int,	VP CPX_PARAM_IMPLBD },
	{ "integrality", sf_dbl,	VP CPX_PARAM_EPINT },
#endif
	{ "iterations",	sf_int,		VP CPX_PARAM_ITLIM },
	{ "iterlim",	sf_int,		VP CPX_PARAM_ITLIM },
	{ "logfile",	sf_char,	VP set_logname },
#ifdef CPLEX_MIP
	{ "lowercutoff", sf_dbl,	VP CPX_PARAM_CUTLO },
#endif
	{ "lowerobj",	sf_dbl,		VP CPX_PARAM_OBJLLIM },
	{ "lowerobjlim", sf_dbl,	VP CPX_PARAM_OBJLLIM },
	{ "lpdisplay",  sf_int,		VP CPX_PARAM_SIMDISPLAY },
	{ "lpiterlim",  sf_int,		VP CPX_PARAM_ITLIM },
	{ "lptimelim",	sf_dbl,		VP CPX_PARAM_TILIM },
	{ "markowitz",	sf_dbl,		VP CPX_PARAM_EPMRK },
	{ "maximize",	sf_known,	VP set_max },
	{ "minimize",	sf_known,	VP set_min },
#ifdef CPLEX_MIP
	{ "mipalg",	sf_int,		VP CPX_PARAM_SUBALG },
	{ "mipalgorithm", sf_int,	VP CPX_PARAM_SUBALG },
	{ "mipcrossover",sf_int,	VP CPX_PARAM_MIPHYBALG },
	{ "mipcuts",	sf_mint,	VP set_mipcuts },
	{ "mipdisplay",	sf_int,		VP CPX_PARAM_MIPDISPLAY },
	{ "mipgap",	sf_dbl,		VP CPX_PARAM_EPGAP },
	{ "mipinterval", sf_int,	VP CPX_PARAM_MIPINTERVAL },
	{ "mipsolutions", sf_int,	VP CPX_PARAM_INTSOLLIM },
	{ "mipstartalg", sf_int,	VP CPX_PARAM_STARTALG },
	{ "mipstartstatus", sf_mint,	VP set_mipststat },
	{ "mipstartvalue", sf_mint,	VP set_mipstval },
	{ "mipsubalg",	sf_int,		VP CPX_PARAM_SUBALG },
	{ "mipthreads",	sf_int,		VP CPX_PARAM_MIPTHREADS },
#endif
	{ "netfind",	sf_int,		VP CPX_PARAM_NETFIND },
	{ "netfinder",	sf_int,		VP CPX_PARAM_NETFIND },
	{ "netopt",	sf_mint,	VP set_netopt },
#ifdef CPLEX_MIP
	{ "node",	sf_int,		VP CPX_PARAM_NODELIM },
	{ "nodefile",	sf_int,		VP CPX_PARAM_NODEFILEIND },
	{ "nodefilelim", sf_dbl,	VP CPX_PARAM_NODEFILELIM },
	{ "nodelim",	sf_int,		VP CPX_PARAM_NODELIM },
	{ "nodes",	sf_int,		VP CPX_PARAM_NODELIM },
	{ "nodesel",	sf_int,		VP CPX_PARAM_NODESEL },
	{ "nodeselect",	sf_int,		VP CPX_PARAM_NODESEL },
	{ "objdifference", sf_dbl,	VP CPX_PARAM_OBJDIF },
#endif
	{ "objno",	sf_mint,	VP set_objno },
	{ "optimality",	sf_dbl,		VP CPX_PARAM_EPOPT },
	{ "optimize",	sf_known,	VP set_primalopt },
#ifdef BARRIER
	{ "ordering",	sf_int,		VP CPX_PARAM_BARORDER },
#endif
#ifdef CPLEX_MIP
	{ "ordertype",	sf_int,		VP CPX_PARAM_MIPORDTYPE },
#endif
	{ "outlev",	sf_int,		VP CPX_PARAM_SIMDISPLAY },
	{ "perturb",	sf_int1,	VP CPX_PARAM_PERIND },
	{ "perturbation", sf_dbl,	VP CPX_PARAM_EPPER },
	{ "perturbconst", sf_dbl,	VP CPX_PARAM_EPPER },
	{ "perturblimit", sf_int,	VP CPX_PARAM_PERLIM },
	{ "pgradient",	sf_int,		VP CPX_PARAM_PPRIIND },
#ifdef CPLEX_MIP
	{ "plconpri",	sf_mint,	VP set_conpri },
	{ "plobjpri",	sf_mint,	VP set_objpri },
#endif
	{ "predual",	sf_int,		VP CPX_PARAM_PREDUAL },
#ifdef CPLEX_MIP
	{ "prerelax",	sf_int1,	VP CPX_PARAM_RELAXPREIND },
#endif
	{ "presolve",	sf_int1,	VP CPX_PARAM_PREIND },
	{ "prestats",	sf_mint,	VP set_prestats },
	{ "pricing",	sf_int,		VP CPX_PARAM_PRICELIM },
	{ "primal",	sf_known,	VP set_primal },
	{ "primalopt",	sf_known,	VP set_primalopt },
#ifdef CPLEX_MIP
	{ "priorities",	sf_int1,	VP CPX_PARAM_MIPORDIND },
	{ "probe",	sf_int,		VP CPX_PARAM_PROBE },
#endif
	{ "readbasis",	sf_char,	VP set_startbas },
#ifdef BARRIER
	{ "readvector",	sf_char,	VP set_startvector },
#endif
	{ "refactor",	sf_int,		VP CPX_PARAM_REINV },
#ifdef CPLEX_MIP
	{ "relax",	sf_known,	VP set_relax },
	{ "relobjdif",	sf_dbl,		VP CPX_PARAM_RELOBJDIF },
	{ "relobjdiff",	sf_dbl,		VP CPX_PARAM_RELOBJDIF },
	{ "rootheuristic", sf_int,	VP CPX_PARAM_HEURISTIC },
	{ "round",	sf_mint,	VP set_round },
#endif
	{ "scale",	sf_int,		VP CPX_PARAM_SCAIND },
	{ "sensitivity", sf_known,	VP set_sens },
	{ "simthreads",	sf_int,		VP CPX_PARAM_SIMTHREADS },
	{ "singular",	sf_int,		VP CPX_PARAM_SINGLIM },
	{ "singularlim", sf_int,	VP CPX_PARAM_SINGLIM },
#ifdef CPLEX_MIP
	{ "solutionlim", sf_int,	VP CPX_PARAM_INTSOLLIM },
	{ "sos",	sf_mint,	VP set_sos },
	{ "sos2",	sf_mint,	VP set_sos2 },
	{ "startalg",	sf_int,		VP CPX_PARAM_STARTALG },
	{ "startalgorithm", sf_int,	VP CPX_PARAM_STARTALG },
	{ "startbasis",	sf_char,	VP set_startbas },
	{ "starttree",	sf_char,	VP set_starttree },
#ifdef BARRIER
	{ "startvector", sf_char,	VP set_startvector },
#endif
#ifdef CPLEX_MIP
	{ "strongcand",	sf_int,		VP CPX_PARAM_STRONGCANDLIM},
	{ "strongit",	sf_int,		VP CPX_PARAM_STRONGITLIM},
	{ "strongthreads", sf_int,	VP CPX_PARAM_STRONGTHREADLIM},
#endif
	{ "subalg",	sf_int,		VP CPX_PARAM_SUBALG },
	{ "subalgorithm", sf_int,	VP CPX_PARAM_SUBALG },
#endif
	{ "time",	sf_dbl,		VP CPX_PARAM_TILIM },
	{ "timelimit",	sf_dbl,		VP CPX_PARAM_TILIM },
	{ "timing",	sf_mint, 	VP set_timing },
	{ "tranopt",	sf_known,	VP set_dualopt },
#ifdef CPLEX_MIP
	{ "treelimit",	sf_dbl,		VP CPX_PARAM_TRELIM },
	{ "treememlim",	sf_dbl,		VP CPX_PARAM_TRELIM },
	{ "treememory",	sf_dbl,		VP CPX_PARAM_TRELIM },
	{ "uppercutoff", sf_dbl,	VP CPX_PARAM_CUTUP },
#endif
	{ "upperobj",	sf_dbl,		VP CPX_PARAM_OBJULIM },
	{ "upperobjlim",sf_dbl,		VP CPX_PARAM_OBJULIM },
#ifdef CPLEX_MIP
	{ "varsel",	sf_int,		VP CPX_PARAM_VARSEL },
	{ "varselect",	sf_int,		VP CPX_PARAM_VARSEL },
#endif
	{ "version",	sf_known,	VP set_version },
	{ "wantsol",	WS_val,		0 },
	{ "writebasis",	sf_char,	VP set_endbas },
	{ "writeprob",	sf_char,	VP set_wrtfname },
#ifdef BARRIER
	{ "writevector", sf_char,	VP set_endvector },
#endif
	{ "xxxstart",	sf_int,		VP CPX_PARAM_XXXIND }
	};

 static Option_Info Oinfo = { "cplex", 0/*cplex_bsname*/, "cplex_options",
				keywds, nkeywds /*, 0, cplex_version+1*/};

 static void
badretfmt /*{*/
#ifdef KR_headers
	(va_alist)
 va_dcl
{
	char *fmt;
	int rc;
	va_list ap;
	va_start(ap);
	rc = va_arg(ap, int);
	fmt = va_arg(ap, char*);
	/*}*/
#else
	(int rc, char *fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
	/*}*/
#endif
{
	char buf[128], *s;
	int k = vsprintf(buf, fmt, ap) + 1;
	if (rc) {
		solve_result_num = rc;
		memcpy(s = (char*)M1alloc(k), buf, k);
		asl->i.uinfo = s;
		}
	if (!amplflag || !rc)
		fprintf(Stderr, "%s\n", buf);
	}}

 static void
#ifdef KR_headers
badret(what, i, rc) char *what; int i, rc;
#else
badret(char *what, int i, int rc)
#endif
{
	badretfmt(rc, "%s failed; error code %d.", what, i);
	}

 static void
#ifdef KR_headers
nonlin(n, rc, what) int n, rc; char *what;
#else
nonlin(int n, int rc, char *what)
#endif
{
	if (n) {
		badretfmt(rc, "%s contains %s.\n", filename, what);
		exit(4);
		}
	}

 static unsigned
#ifdef KR_headers
namemem(n) int n;
#else
namemem(int n)
#endif
{
	int k, L;
	unsigned rv;

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
#ifdef KR_headers
namegen(stub, name, s, n) char *stub, **name, *s; int n;
#else
namegen(char *stub, char **name, char *s, int n)
#endif
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

 static void
#ifdef KR_headers
make_names(cnamep, rnamep) char ***cnamep, ***rnamep;
#else
make_names(char ***cnamep, char ***rnamep)
#endif
{
	char **cname, **rname, *s;
	unsigned cstorsz, rstorsz;

	cstorsz = namemem(nbas);
	rstorsz = namemem(mbas);
	*cnamep = cname = (char **)Malloc((mbas+nbas)*sizeof(char *)
				+ cstorsz + rstorsz);
	*rnamep = rname = cname + nbas;
	s = namegen("x", cname, (char *)(rname + mbas), nbas);
	namegen("c", rname, s, mbas);
	}

#ifdef CPLEX_MIP
 static int
#ifdef KR_headers
nzeros(p, n) int *p; int n;
#else
nzeros(int *p, int n)
#endif
{
	int *pe = p + n;
	n = 0;
	while(p < pe)
		if (*p++)
			n++;
	return n;
	}

 static void
#ifdef KR_headers
mip_priorities(asl, cpx) ASL *asl; cpxlp *cpx;
#else
mip_priorities(ASL *asl, cpxlp *cpx)
#endif
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
	else if (nnames = mip_pri(&start, &num, &pri, 2147483647)) {
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
#ifdef KR_headers
set_mipinit(asl, cpx, nint) ASL *asl; cpxlp *cpx; int nint;
#else
set_mipinit(ASL *asl, cpxlp *cpx, int nint)
#endif
{
	int i, m, *r;
	real *L, *U, *v;

	r = (int*)Malloc(nint*sizeof(int));
	m = n_var - (niv + nbv);
	v = X0 + m;
	L = LUv + m;
	U = Uvx + m;
	for(i = 0; i < nint; i++) {
		r[i] = m + i;
		if (v[i] < L[i])
			v[i] = L[i];
		else if (v[i] > U[i])
			v[i] = U[i];
		}
	if (i = CPXcopymipstart(Env, cpx, nint, r, v))
		badret("CPXcopymipstart", i, 0);
	else if (i =  CPXsetintparam(Env, CPX_PARAM_MIPSTART, 1))
		badret("CPXsetintparam(CPX_PARAM_MIPSTART)", i, 0);
	free(r);
	}
#endif

 static int
#ifdef KR_headers
parval(par) int par;
#else
parval(int par)
#endif
{
	int rv;
	CPXgetintparam(Env, par, &rv);
	return rv;
	}

 static int CPXPUBLIC
#ifdef KR_headers
lpcbf(env, lp, wf, cbh) CPXENVptr env; char *lp, *cbh; int wf;
#else
lpcbf(CPXENVptr env, void *lp, int wf, void *cbh)
#endif
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
	static char pfmt[] =
		"%sLP Presolve eliminated %d rows and %d columns.\n";
	static int wf0;
	char *pfx, *prob;
	Not_Used(cbh);

	switch(wf) {
	  case CPX_CALLBACK_PRESOLVE:
		if (!prestats || !dispzap)
			goto done;
		for(i = 0; i < 4; i++)
			CPXgetcallbackinfo(env, lp, wf, wi[i], &n[i]);
		wf0 = wf;
		goto done;
	  case CPX_CALLBACK_PRIMAL_CROSSOVER:
		cbi.xkind = "primal";
		for(i = 0; i < 2; i++)
			CPXgetcallbackinfo(env, lp, wf, wi[i+4], &cbi.nx[i]);
		goto done;
	  case CPX_CALLBACK_DUAL_CROSSOVER:
		cbi.xkind = "dual";
		for(i = 0; i < 2; i++)
			CPXgetcallbackinfo(env, lp, wf, wi[i+6], &cbi.nx[i]);
		goto done;
	  }
	if (wf == wf0 || !prestats || !dispzap)
		goto done;
	wf0 = wf;
	if (cbi.npr++ && !memcmp(n, cbi.n0, 4*sizeof(int)))
		goto done;
	memcpy(cbi.n0, n, 4*sizeof(int));
	if (cbi.nint) {
		pfx = "MIP ";
		prob = "MIP";
		}
	else {
		pfx = "";
		prob = "LP";
		}
	if (!(n[0] + n[1] + n[2])) {
		printf("No %s presolve or aggregator reductions.\n", prob);
		goto done;
		}
	if (cbi.pres)
		printf(pfmt, pfx, n[0], n[1]);
	if (n[2])
		printf("Aggregator did %d substitutions.\n", n[2]);
	if (n[3])
		printf("%d coefficients modified.\n", n[3]);
 done:
	return breaking > 1;
	}

 static void
#ifdef KR_headers
stat_map(stat, n, map, mx, what)
	int *stat; int n; int *map; int mx; char *what;
#else
stat_map(int *stat, int n, int *map, int mx, char *what)
#endif
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
#ifdef KR_headers
get_statuses(asl, cpx, d) ASL *asl; cpxlp *cpx; dims *d;
#else
get_statuses(ASL *asl, cpxlp *cpx, dims *d)
#endif
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
	if (!havestats || !mipststat && niv + nbv)
		return;
	if (method > 0) {
		stat_map(d->cs = d->cstat, n, map, 6, "incoming cstat");
		stat_map(d->rs = d->rstat, m, map, 6, "incoming rstat");
		if (objadj) {
			d->cstat[n] = CPX_BASIC;
			d->rstat[m] = CPX_FREE_SUPER;
			}
		if (i = CPXcopybase(Env, cpx, d->cstat, d->rstat)) {
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
#ifdef KR_headers
qmatadj(k, nr, os, colq, colqcnt, qmatp)
	int k, nr, os, *colq, *colqcnt; double **qmatp;
#else
qmatadj(int k, int nr, int os, int *colq, int *colqcnt, double **qmatp)
#endif
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

 static void
#ifdef KR_headers
sos_kludge(nsos, sosbeg, sosref) int nsos; int *sosbeg; double *sosref;
#else
sos_kludge(int nsos, int *sosbeg, double *sosref)
#endif
{
	/* Adjust sosref if necessary to accommodate CPLEX's */
	/* undocumented requirement that sosref values differ */
	/* by at least 1e-10. */
	int i, j, k;
	double t, t1;
	for(i = j = 0; i++ < nsos; ) {
		k = sosbeg[i];
		t = sosref[j];
		while(++j < k) {
			t1 = sosref[j];
			t += 1e-10;
			if (t1 <= t)
				sosref[j] = t1 = t + 1e-10;
			t = t1;
			}
		}
	}

 static void
#ifdef KR_headers
amplin(asl, pcpx, nl, d, nintp, av)
	ASL *asl; cpxlp **pcpx, FILE **nl; dims *d; int *nintp; char **av;
#else
amplin(ASL *asl, cpxlp **pcpx, FILE **nl, dims *d, int *nintp, char **av)
#endif
{
	int *ia, *ia1, *ia2, *ja, *ka, *ka1, *ka2, *kal, *kal1, *kal2;
	double *a, *a1, *b, *b1, *c, *c1, *l, *l1,
		*le, *lrhs, *lx, *rngvec, *rv1, *u, *u1, *urhs, *ux;
	double nos, os;
	ograd *og, *og1;
	int i, j, k, m, m0, n, n0, nbv1, nextra, nint, nint1,
		nr, nsos, nz, nz1, nzc0, nzextra, nzr;
	int copri[2];
	char **cname, **rname, *s, *senx;
	Char **zap[4];
#ifdef BARRIER
	fint nelqf, *colqf, *rowqf;
	int havestart, nelq, *colq, *colqcnt, *rowq;
	double *cdual, *cprim, *qmat, *rdual, *rprim;
#endif
	char probname[48];
	cpxlp *cpx;
#ifdef CPLEX_MIP
	char *ctype;
	char *sostype;
	int nsosnz, *sosbeg, *sosind, *sospri;
	double *sosref;
#endif
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
	Optimize = CPXprimopt;
	algname = baralgname = "";
	aggtries = 0;

	if (!(mint_val[set_objno].U = n_obj))
		objno = 0;
	nonlin(nlc, 550, "nonlinear constraints");
	nonlin(nlvoi, 552, "nonlinear integer variables");
	nonlin(plterms, 554, "piecewise-linear terms");

	if (getopts(av, &Oinfo)) {
		if (amplflag)
			badretfmt(560, "Error in $cplex_options.");
		exit(1);
		}
	obj_no = objno - 1;

	m = n_con + 1;	/* allow elbow room for obj_adj */
	n = n_var + 1;
	nzr = nzc + 1;

	d->cstat = (int*)M1zapalloc((m+n+2)*sizeof(int));
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
	zap[1] = zap[2] = zap[3] = 0;
	ia = A_rownos = (int *)(a + nzr);
	ka = A_colstarts = ia + nzr;
#ifdef BARRIER
	nint = nbv + niv;
	if (Optimize == Optimize2 || Optimize == Optimizebar || nlo
	    || nint && mipstval)
		want_xpi0 = 3;
	want_deriv = 0;
	qp_read(*nl,0);
	*nl = 0;	/* Indicate closed file -- in case we're */
			/* subsequently interrupted. */
	nelqf = 0;
	if (obj_no >= 0 && obj_no < n_obj
	 && (nelqf = mqpcheck(obj_no, &rowqf, &colqf, &qmat))) {
		if (nelqf < 0) {
			nonlin(nelqf == -2, 555,
			 "a quadratic objective involving division by 0");
			nonlin(1, 551, "a nonlinear objective");
			}
		nonlin(nint, 553,
			"integer variables and a quadratic objective");
		crossover = 0;
		set_baropt();
		method = 1;
		}
#else
	f_read(*nl, 0);
	*nl = 0;
#endif

	nsos = 0;
#ifdef CPLEX_MIP
	if (!relax) {
		copri[0] = objpri;
		copri[1] = conpri;
		i = ASL_suf_sos_explict_free;
		if (!sos)
			i |= ASL_suf_sos_ignore_sosno;
		if (!sos2)
			i |= ASL_suf_sos_ignore_amplsos;
		nsos = suf_sos(i, &nsosnz, &sostype, &sospri, copri,
				&sosbeg, &sosind, &sosref);
		}
#endif
	m = m0 = n_con;
	n = n0 = nr = n_var;
	nzc0 = nzc;
	nint = nbv + niv;
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
	*nintp = nint1 = nint + nsos;
	if (use_netopt == 1 && !lnc)
		use_netopt = 0;
	if (nint1 || m <= 0) {
		method = 1;
		use_netopt = 0;
		}
	else if (!method)
		method = m - n > dual_thresh && !use_netopt ? -1 : 1;
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
	if (!(cpx = CPXcreateprob(Env, &i, probname))) {
		badret("CPXcreateprob", i, 531);
		return;
		}
	*pcpx = cpx;
	n = n0;
	m = m0;
	nextra = 0;
	if (method < 0) { /* solve dual */
		l = LUv;
		u = Uvx;
		for(le = l + n0; l < le; l++, u++) {
			if (*l > negInfinity && *l)
				nextra++;
			if (*u < Infinity && *u)
				nextra++;
			}
		nzextra = 0;
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
		nr += objadj;
		nz1 = nz;
		l = LUv;
		u = Uvx;
		if (nint1) {
#ifdef CPLEX_MIP
			ctype = M1alloc(nr);
			if (i = n - nbv1 - niv - nextra)
				memset(ctype, 'C', i);
			if (nbv1)
				memset(ctype+i, 'B', nbv1);
			if (niv)
				memset(ctype+i+nbv1, 'I', niv);
			if (nextra)
				memset(ctype+i+nbv1+niv, 'C', nextra);
#else
			badretfmt(557,
		"CPLEX's MIP option is needed to solve problems with\n%s.\n",
			"nonconvex (nonconcave) piecewise-linear terms");
			exit(1);
#endif
			}
		}

	d->x = (double *)M1zapalloc((m+2L*nr)*sizeof(double) + m);
	d->y = d->x + nr;
	d->c = c = d->y + m;
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
			if (*--ux < Infinity)
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
			if (*--lx > negInfinity)
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
			*--ia1 = m0;
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
	if (startbas || endbas || startvec || endvec
	 || Optimize == CPXprimopt || Optimize == CPXdualopt
	 || parval(CPX_PARAM_XXXIND)) {
		make_names(&cname, &rname);
		zap[3] = M1record((Char*)cname);
		}
	i = CPXcopylpwnames(Env, cpx, n, m, objsen, c, b, senx, ka,
				kal, ia, a, l, u, rngvec, cname, rname);
	if (i) {
		badret("CPXcopylpwnames", i, 531);
 bailout:
		CPXfreeprob(Env, pcpx);
		return;
		}
#ifdef CPLEX_MIP
	if (nint1) {
		if (i = CPXcopyctype(Env, cpx, ctype)) {
			badret("CPXcopyctype", i, 531);
			goto bailout;
			}
		if (nsos) {
			sos_kludge(nsos, sosbeg, sosref);
			i = CPXcopysos(Env,cpx,nsos,nsosnz,sostype,sospri,
					sosbeg,sosind,sosref);
			free(sosref);
			if (i) {
				badret("CPXcopysos", i, 531);
				goto bailout;
				}
			}
		if (nint) {
			mip_priorities(asl, cpx);
			if (mipstval && X0 && (i = nbv1 + niv))
				set_mipinit(asl, cpx, i);
			}
		}
	else
#endif
		{
#ifdef BARRIER
		havestart = 0;
		if (Optimize == Optimize2 || Optimize == Optimizebar)
			havestart = startcomp(asl, n, m, nextra, ka, kal, ia,
				a, b, c, &cdual, &cprim, &rdual, &rprim);
#endif
#ifdef BARRIER
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
						goto bailout;
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
					if (i) {
						badret("CPXcopyquad",
							i, 531);
						goto bailout;
						}
					baralgname = "QP ";
					break;
					}
				}
			}
		if (havestart) {
			if (i = CPXcopystart(Env,cpx,0,0,cprim,rprim,cdual,rdual))
				badret("CPXcopystart", i, 0);
			}
#endif
		}
	for(i = 4; i > 0; )
		if (zap[--i]) {
			free(*zap[i]);
			*zap[i] = 0;
			}
	get_statuses(asl, cpx, d);
	cbi.disp = parval(CPX_PARAM_SIMDISPLAY);
#ifdef CPLEX_MIP
	cbi.mipdisp = parval(CPX_PARAM_MIPDISPLAY);
#endif
	if (prestats && (dispzap = !cbi.disp || nint && !cbi.mipdisp)
#ifndef NO_BARRIER
	 || Optimize == Optimize2
#endif
				 ) {
		zap_lpcbf = 1;
		cbi.agg = parval(CPX_PARAM_AGGIND);
		cbi.m = m;
		cbi.n = n;
		cbi.nint = nint;
		cbi.pres = parval(CPX_PARAM_PREIND);
		CPXsetlpcallbackfunc(Env, lpcbf, 0);
		CPXsetmipcallbackfunc(Env, lpcbf, 0);
		}
	else {
		CPXsetlpcallbackfunc(Env, breakcallback, 0);
		CPXsetmipcallbackfunc(Env, breakcallback, 0);
		}
	breaking = 1;
	}

 static int
#ifdef KR_headers
round(x, n, assign, dp) double *x, *dp; fint n; int assign;
#else
round(double *x, fint n, int assign, double *w)
#endif
{
	double d, dx, *xe, y;
	int m = 0;

	dx = 0;
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
	if (m) {
		CPXinfodblparam(Env, CPX_PARAM_EPINT, w+1, w+2, w+3);
		if (dx <= w[2] && !(Round & 8))
			m = 0;
		else
			CPXgetdblparam(Env, CPX_PARAM_EPINT, w+1);
		}
	return m;
	}

 static char *
#ifdef KR_headers
failstat(j) int j;
#else
failstat(int j)
#endif
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
	  case CPXERR_BOUNDS_INFEAS:
		return "contradictory bounds";
	  case CPXERR_PRESLV_INForUNBD:
		return "presolve: infeasible or unbounded problem";
	  case CPXERR_INDEX_RANGE:
		return "illegal list length";
	  case CPXERR_NET_IMBALANCE:
		return "an infeasible or unbounded network";
	  case CPXERR_BAD_METHOD:
		return "unknown method after network optimization";
	  }
	sprintf(buf, "CPLEX error # %d.", j);
	return buf;
	}

 static int
#ifdef KR_headers
statadjust(stat) int stat;
#else
statadjust(int stat)
#endif
{
	if (stat > 100 && stat < 120)
		stat -= 69;
	else if (stat > 18 && stat < 44)
		stat -= 13;
	else if (stat == 1101)
		stat = 31;
	else if (stat < 0 || stat > 18)
		stat = 0;
	return stat;
	}

 static void
#ifdef KR_headers
write_basis(cpx) cpxlp *cpx;
#else
write_basis(cpxlp *cpx)
#endif
{
	if (CPXmbasewrite(Env, cpx, endbas)) {
		printf("Could not write ending basis to \"%s\"\n", endbas);
		need_nl = 0;
		}
	}

 static void
#ifdef KR_headers
cud_cadjust(asl) ASL *asl;
#else
cud_cadjust(ASL *asl)
#endif
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
#ifdef KR_headers
Kludge_getbase(d, rstat) dims *d; int *rstat;
#else
Kludge_getbase(dims *d, int *rstat)
#endif
{ /* Kludge around CPLEX bug in rstat values for <= constraints... */
	int i, m = mbas;
	char *rt = d->rtype;
	for(i = 0; i < m; i++)
		if (rt[i] == 'L' && rstat[i] == 0)
			rstat[i] = 2;
	}

 static int
#ifdef KR_headers
send_statuses(asl, cpx, d) ASL *asl; cpxlp *cpx; dims *d;
#else
send_statuses(ASL *asl, cpxlp *cpx, dims *d)
#endif
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
		if (i = CPXgetbase(Env, cpx, cstat, rstat)) {
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
			if (i = CPXobjsa(Env, cpx, 0, nbas-1, L, U))
				badret("CPXobjsa", i, 0);
			else {
				suf_rput("down", ASL_Sufkind_var, L);
				suf_rput("up", ASL_Sufkind_var, U);
				suf_rput("current", ASL_Sufkind_var, d->c);
				}
			L = U + nbas;
			U = L + mbas;
			if (i = CPXrhssa(Env, cpx, 0, mbas-1, L, U))
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
	if (i = CPXgetbase(Env, cpx, cs, rs))
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
#ifdef KR_headers
send_ray(asl, cpx) ASL *asl; cpxlp *cpx;
#else
send_ray(ASL *asl, cpxlp *cpx)
#endif
{
	int i, j, n, rc;
	int *cs, *rs;
	real *r, scale, *x;

	if (CPXgetijdiv(Env, cpx, &i, &j))
		return 1;
	if ((n = n_var) < nbas)
		n = nbas;
	r = (real*)M1zapalloc(n*sizeof(real));
	rc = 1;
	x = (real*)Malloc(mbas*sizeof(real) + (mbas+nbas)*sizeof(int));
	cs = (int*)(x + mbas);
	rs = cs + nbas;
	if (CPXgetbase(Env, cpx, cs, rs))
		goto done;
	scale = -1.;
	if (i >= 0) {
		if (CPXbinvcol(Env, cpx, i, x))
			goto done;
		if (rs[i] == CPX_AT_UPPER)
			scale = 1.;
		}
	else if (j >= 0) {
		if (CPXbinvacol(Env, cpx, j, x))
			goto done;
		if (cs[j] == CPX_AT_UPPER)
			scale = 1.;
		r[j] = -scale;
		}
	else
		goto done;
	for(i = 0; i < nbas; i++)
		if (cs[i] == CPX_BASIC) {
			if (CPXgetijrow(Env, cpx, -1, i, &j))
				goto done;
			r[i] = scale*x[j];
			}
	suf_rput("unbdd", ASL_Sufkind_var, r);
	rc = 0;
 done:
	free(x);
	return rc;
	}

 static void
#ifdef KR_headers
iis_put(asl, kind, n, ind, stat) ASL *asl; int kind; int n; int *ind; int *stat;
#else
iis_put(ASL *asl, int kind, int n, int *ind, int *stat)
#endif
{
	int i, *x;

	x = (int*)M1zapalloc((&asl->i.n_var_)[kind]*sizeof(int));
	suf_iput("iis", kind, x);
	for(i = 0; i < n; i++)
		x[ind[i]] = stat[i] + 1;
	}

 static int
#ifdef KR_headers
send_iis(asl, cpx, hbuf) ASL *asl; cpxlp *cpx; char *hbuf;
#else
send_iis(ASL *asl, cpxlp *cpx, char *hbuf)
#endif
{
	int i, j, m, m1, n, n1;
	int *ci, *cs, *ri, *rs;

	if (!want_iis)
		return 0;
	CPXsetintparam(Env, CPX_PARAM_IISIND, want_iis-1);
	if (i = CPXfindiis(Env, cpx, &m, &n))
		return Sprintf(hbuf, "\nReturn %d from CPXfindiis.", i);
	ci = (int*)M1alloc((2*(m+n))*sizeof(int));
	cs = ci + n;
	ri = cs + n;
	rs = ri + m;
	if (i = CPXgetiis(Env, cpx, &j, ri, rs, &m1, ci, cs, &n1))
		return Sprintf(hbuf, "\nReturn %d from CPXgetiis.", i);
	iis_put(asl, ASL_Sufkind_var, n, ci, cs);
	iis_put(asl, ASL_Sufkind_con, m, ri, rs);
	return Sprintf(hbuf,
		"\nReturning %siis of %d variables and %d constraints.",
		j == 1 ? "" : "partial ", n, m);
	}

 static void
#ifdef KR_headers
amplout(asl, cpx, d, status, nint1, nosp)
	ASL *asl; cpxlp *cpx; dims *d; int status, nint1, *nosp;
#else
amplout(ASL *asl, cpxlp *cpx, dims *d, int status, int nint1, int *nosp)
#endif
{
	char buf[32], hbuf[640], *wb;
	double bobj, bobj0, obj, *l, *le, *u, w[4], *x, *x2, *y, *z, *z1;
	int bitc, ccpr, cl0, cl, cov, i, ii, itc, itci, m, mipstat, n;
	int nint, nodecnt, nodelim, nos, nround, stat, stat0, stat1, stat10;
	typedef struct { char *msg; int code; int wantobj; } Sol_info;
	Sol_info *SI;
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
	 { "aborted in phase II", 501, 0 },
	 { "aborted in phase I", 502, 0 },
	 { "aborted in barrier, dual infeasible", 503, 0 },	/* 14* */
	 { "aborted in barrier, primal infeasible", 504, 0 },
	 { "aborted in barrier, primal and dual infeasible", 505, 0 },
	 { "aborted in barrier, primal and dual feasible", 506, 0 },
	 { "aborted in crossover", 507, 0 },			/* 18 */

	 { "converged, dual feasible, primal infeasible", 204, 1 }, /*32 */
	 { "converged, primal feasible, dual infeasible", 301, 1 },
	 { "converged, primal and dual infeasible", 205, 1 },
	 { "primal objective limit reached", 405, 1 },		/* 35 */
	 { "dual objective limit reached", 406, 1 },
	 { "primal has unbounded optimal face", 001, 1 },
	 { "best solution found, primal-dual feasible", 100, 1 },/* 38 */
	 { "best solution found, primal infeasible", 206, 1 },
	 { "best solution found, dual infeasible", 302, 1 },
	 { "best solution found, primal-dual infeasible", 207, 1 }, /* 41 */
	 { "solution found, numerical difficulties", 508, 1},
	 { "solution found, inconsistent equations", 509, 1},  /* 43 */

	 { "infeasible or unbounded in presolve", 208, 0 },	/*1101*/
#ifdef CPLEX_MIP
	 { "optimal integer solution", 002, 1 },		/* 101 */
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
	 { "node file limit with integer solution", 424, 1 },
	 { "node file limit with no integer solution", 413, 1 }, /* 119 */
	 { "optimal (non-)integer solution", 102, 1 },
	 { "optimal (non-)integer solution within mipgap or absmipgap", 103, 1 }
,
#endif
		};
	static Sol_info solinfo1[] = {
	 { "QP Hessian has diag. elements of the wrong sign", 541, 0 },
	 { "QP Hessian is not positive semi-definite", 542, 0 }
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
		"a bug" };

	stat1 = 1;
	ii = mipstat = nodecnt = nround = 0;
#ifdef CPLEX_MIP
	ccpr = 0;
	if (Optimize == CPXmipopt) {
		if (prestats && nint1
		 && (parval(CPX_PARAM_CLIQUES) || parval(CPX_PARAM_COVERS))) {
			cl0 = CPXgetgenclqcnt(Env, cpx);
			cl  =  CPXgetclqcnt(Env, cpx);
			cov =  CPXgetcovcnt(Env, cpx);
			ccpr = 1;
			}
		if (status > 0) {
			mipstat = statadjust(stat0 = status);
			status = 0;
			}
		}
#endif
	m = n_con;
	n = n_var;
	x = d->x;
	y = d->y;
#ifdef BARRIER
	bitc = CPXgetbaritcnt(Env, cpx);
#else
	bitc = 0;
#endif
	itc = CPXgetitcnt(Env, cpx);
	itci = CPXgetphase1cnt(Env, cpx);
	nos = nosp ? *nosp : 0;
	switch (status) {
	  case 0: break;
	  case CPXERR_Q_NEG_ZERO_COMP:
		SI = solinfo1;
		goto have_SI;
	  case CPXERR_Q_NOT_POS_DEF:
		SI = solinfo1 + 1;
		goto have_SI;
	  case CPXERR_PRESLV_INForUNBD:
		SI = solinfo + 31;
		goto have_SI;
	  default:
		stat0 = status;
		stat = 0;
		x = y = 0;
		goto have_stat;
		}
#ifdef CPLEX_MIP
	bobj = bobj0 = method * objsen * Infinity;
	if (nint1) {
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
		else
			stat = statadjust(stat0 = CPXgetstat(Env, cpx));
		ii = CPXgetmipitcnt(Env, cpx);
		nodecnt = CPXgetnodecnt(Env, cpx);
 have_ii:
		if (CPXgetbestobjval(Env, cpx, &bobj))
			bobj = bobj0;
		if (solinfo[stat].wantobj || endbas) {
			if (zap_lpcbf) {
				zap_lpcbf = 0;
				CPXsetlpcallbackfunc(Env, breakcallback, 0);
				CPXsetmipcallbackfunc(Env, breakcallback, 0);
				}
			CPXgetmipobjval(Env, cpx, &obj);
			CPXchgprobtype(Env, cpx, 3);
			status = CPXprimopt(Env, cpx);
			if (status)
				fprintf(Stderr,
				"CPLEX status %d with fixed integers:\n\t%s\n",
					status, solinfo[statadjust(status)].msg);
			if (endbas)
				write_basis(cpx);
			if (!solinfo[stat].wantobj) {
				x = y = 0;
				goto have_stat;
				}
			CPXsolution(Env, cpx, &stat1, &obj, x, y, 0, 0);
			stat1 = statadjust(stat10 = stat1);
			/* round to integer */
			if (!relax && Round >= 0) {
				nint = niv + nbv;
				x2 = x + n_var - nint;
				if (nround = round(x2, nint, Round & 1, w)) {
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
		CPXsolution(Env, cpx, &stat, &obj, x, y, 0, 0);
		stat = nos == 2 ? 2 : statadjust(stat0 = stat);
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
	i = Sprintf(hbuf, "%s: %s", Oinfo.bsname, SI->msg);
	if (!stat)
		i += Sprintf(hbuf+i, ": %s", failstat(stat0));
#ifdef CPLEX_MIP
	else if (stat == 16 || stat == 17) {
		CPXgetintparam(Env, CPX_PARAM_NODELIM, &nodelim);
		i += Sprintf(hbuf+i, ".\nCurrent node limit = %d", nodelim);
		}
#endif
	else if (!SI->wantobj)
		i += Sprintf(hbuf+i, ".");
	if (SI->wantobj) {
		g_fmtop(buf, obj);
		i += Sprintf(hbuf+i, "; objective %s", buf);
		}
	solve_result_num = SI->code;
	if (nosp) {
		if (nos < 0 || nos > 11)
			nos = 11;
		i += nos == 1
			? Sprintf(hbuf+i, "\n%d network simplex iterations.",
				netiters)
			: Sprintf(hbuf+i, "\nnetopt found %s.", netmsg[nos]);
		if (net_nodes != lnc || net_arcs != nwv)
			i += Sprintf(hbuf+i,
			"\nNetwork extractor found %d nodes and %d arcs.",
				net_nodes, net_arcs);
		}
	if (nint1 && !mipstat)
		i += Sprintf(hbuf+i,
		 "\n%d MIP simplex iterations\n%d branch-and-bound nodes",
			ii, nodecnt);
	if (itc > 0 || itci > 0
	 || Optimize == CPXprimopt || Optimize == CPXdualopt)
		i += Sprintf(hbuf+i,
			"\n%d %ssimplex iterations (%d in phase I)",
				itc, algname, itci);
	if (bitc > 0) {
		i += Sprintf(hbuf+i, "\n%d %sbarrier iterations",
			bitc, baralgname);
		if (cbi.nx[0] | cbi.nx[1])
			i += Sprintf(hbuf+i,
				"\n%d push, %d exchange %s crossover iterations",
				cbi.nx[0], cbi.nx[1], cbi.xkind);
		}
#ifdef CPLEX_MIP
	if (bestnode) {
		if (n_obj > 0) {
			u = (real*)M1zapalloc(n_obj*sizeof(real));
			if (obj_no >= 0)
				u[obj_no] = bobj;
			suf_rput("bestnode", ASL_Sufkind_obj, u);
			}
		suf_rput("bestnode", ASL_Sufkind_prob, &bobj);
		}
	if (ccpr) {
		if (cl0 | cl | cov) {
			if (cl0)
				i += Sprintf(hbuf+i,
			  "\n%d of %d clique inequalities used", cl, cl0);
			if (cov)
				i += Sprintf(hbuf+i,
			  "\n%d cover cuts added", cov);
			}
		else
			i += Sprintf(hbuf+i, "\nNo clique or cover cuts used.");
		}
#endif
	if (stat1 != 1) {
		i += Sprintf(hbuf+i, "\nStatus recovering solution: %s",
			 solinfo[stat1].msg);
		if (!stat1)
			i += Sprintf(hbuf+i, ": stat = %d", stat10);
		}
	if (aggtries > 1)
		i += Sprintf(hbuf+i, "\nTried aggregator %d times", aggtries);
	if (solve_result_num >= 200 && solve_result_num < 300)
		i += send_iis(asl, cpx, hbuf+i);
	if (send_statuses(asl, cpx, d))
		i += Sprintf(hbuf+i, "; no basis.");
	else if (stat == CPX_UNBOUNDED && asl->i.flags & 1 && method > 0)
		i += Sprintf(hbuf+i, send_ray(asl,cpx)
			? "\nfailed to compute variable.unbdd"
			: "\nvariable.unbdd returned");
	if (nround) {
		wb = "";
		if (nround < 0) {
			nround = -nround;
			wb = "would be ";
			}
		i += Sprintf(hbuf+i,
			"\n%d integer variables %srounded (maxerr = %g).",
			nround, wb, w[0]);
		if (w[0] > w[2]) {
			i += Sprintf(hbuf+i,
			  "\nAssigning integrality = %g might help.", w[2]);
			i += Sprintf(hbuf+i,
			  "\nCurrently integrality = %g.", w[1]);
			}
		}
	write_sol(hbuf, x, y, &Oinfo);
	}

 static double Times[4];

 static void
#ifdef KR_headers
show_times()
#else
show_times(void)
#endif
{
	int i;

	Times[3] = xectim_();
	for(i = 1; i <= 2; i++)
	    if (time_flag & i) {
		fprintf(i == 1 ? stdout : Stderr,
		"\nTimes (seconds):\nInput =  %g\nSolve =  %g\nOutput = %g\n",
			Times[1] - Times[0], Times[2] - Times[1],
			Times[3] - Times[2]);
		}
	}

 static void
#ifdef KR_headers
basread(asl, cpx) ASL *asl; cpxlp *cpx;
#else
basread(ASL *asl, cpxlp *cpx)
#endif
{
	if (CPXreadcopybase(Env, cpx, startbas)) {
		printf("Could not read starting basis file \"%s\"\n",
			startbas);
		need_nl = 0;
		}
	}

 static void CPXPUBLIC
#ifdef KR_headers
mymsgfunc(handle, msg) char *handle, *msg;
#else
mymsgfunc(void *handle, char *msg)
#endif
{
	static int wantnl;
	char *msg1;
	int n;

	handle = handle; /* shut up non-use warning */
	if (Logf)
		CPXfputs(msg, Logf);
	if (dispzap) {
		if (!strncmp(msg, "MIP start va", 12))
			goto print_msg;
		n = strlen(msg);
		if (n > 20 && !strcmp(msg+n-20, "variable selection.\n"))
			goto print_msg;
		if (*msg != '\n')
			return;
		dispzap = 0;
		}
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
3	upp	at upper bound\n";

 static SufDecl
suftab[] = {
	{ "bestnode", 0, ASL_Sufkind_obj  | ASL_Sufkind_outonly },
	{ "bestnode", 0, ASL_Sufkind_prob | ASL_Sufkind_outonly },
	{ "current", 0, ASL_Sufkind_con | ASL_Sufkind_outonly },
	{ "current", 0, ASL_Sufkind_var | ASL_Sufkind_outonly },
	{ "direction", 0, ASL_Sufkind_var },
	{ "down", 0, ASL_Sufkind_con | ASL_Sufkind_outonly },
	{ "down", 0, ASL_Sufkind_var | ASL_Sufkind_outonly },
	{ "iis", iis_table, ASL_Sufkind_var | ASL_Sufkind_outonly },
	{ "iis", 0, ASL_Sufkind_con | ASL_Sufkind_outonly },
	{ "priority", 0, ASL_Sufkind_var },
	{ "ref", 0, ASL_Sufkind_var | ASL_Sufkind_real },
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

 static char *
#ifdef KR_headers
strcpy1(t, s) char *t; char *s;
#else
strcpy1(char *t, char *s)
#endif
{
	while(*t = *s++)
		t++;
	return t;
	}

 static void
#ifdef KR_headers
adjust_version(asl) ASL *asl;
#else
adjust_version(ASL *asl)
#endif
{
	char *d, *de, *s, *t;
	int cvlen, dvlen;

	for(d = cplex_version; *d != '('; d++);
	for(de = d; *de++ != ')'; );
	*de = 0;
	dvlen = de - d;
	s = CPXversion(Env);
	cvlen = strlen(s);
	Oinfo.bsname = M1alloc(2*cvlen + dvlen + 24);
	t = strcpy1(Oinfo.bsname, "CPLEX ");
	t = strcpy1(t, s);
	Oinfo.version = version = ++t;
	t = strcpy1(t, Oinfo.bsname);
	t = strcpy1(t, ", driver");
	t = strcpy1(t, d);
	strcpy1(t, "\n");
	}

 static int
#ifdef KR_headers
Optimize1(e, c) CPXENVptr e; cpxlp *c;
#else
Optimize1(CPXENVptr e,  cpxlp *c)
#endif
{
	int rc = Optimize(e, c);
	if (rc == 1101) {
		CPXsetintparam(e, CPX_PARAM_AGGIND, 0);
		CPXsetintparam(e, CPX_PARAM_PREIND, 0);
		rc = Optimize(e, c);
		}
	return rc;
	}

 int
#ifdef KR_headers
main(argc, argv) int argc; char **argv;
#else
main(int argc, char **argv)
#endif
{
	int	nint1, nos, *nosp, rc, status, z;
	dims	d;
	char	*stub;
	CPXCHANNELptr cpxresults;
	 /* static values for use with longjmp */
	static cpxlp *cpx;
	static sig_func_type *oic;
	static FILE *nl;

	Times[0] = xectim_();

#ifdef BARRIER
	asl = ASL_alloc(ASL_read_fg);
#else
	asl = ASL_alloc(ASL_read_f);
#endif

	Env = 0;
	cpx = 0;
	breaking = 0;
	oic = 0;
	nl = 0;
	rc = setjmp(Jb);
	if (rc) {
		if (nl)
			fclose(nl);
		if (solve_result_num > 0) {
			rc = 0;
			goto ws_now;
			}
		--rc;
		goto done;
		}
	if (argc < 2)
#ifdef Stand_alone
		return cpxmain(argc,argv);
#else
		usage_ASL(&Oinfo, 1);
#endif
	if (Env = CPXopenCPLEXdevelop(&z))
		adjust_version(asl);

	/* prepare to return error messages in .sol file */
	if (!(stub = getstub(&argv, &Oinfo)))
		usage_ASL(&Oinfo, 1);
	nl = jac0dim(stub, (fint)0);
	oic = signal(SIGINT, intcatch);

	if (!Env) {
		switch(z) {
		  case 32020:
			badretfmt(572, "Expired CPLEX license.");
			break;
		  case 32027:
			badretfmt(571,
			"CPLEX licensing limit: maximum processes exceeded.");
			break;
		  default:
			badretfmt(570,
	"CPLEX licensing problem: error code %d from CPXopenCPLEXdevelop.",
				z);
		  }
		if (amplflag)
			goto ws_now;
		goto done;
		}

	if (z = CPXgetchannels(Env, &cpxresults, 0,0,0)) {
		badret("CPXgetchannels", z, 531);
		goto ws_now;
		}

	suf_declare(suftab, sizeof(suftab)/sizeof(SufDecl));
	amplin(asl, &cpx, &nl, &d, &nint1, argv);
	if (!cpx) {
		if (solve_result_num > 0) {
 ws_now:
			write_sol((char*)asl->i.uinfo, 0, 0, &Oinfo);
			}
		else {
			printf("Failed to load the problem!\n");
			rc = 1;
			}
		goto done;
		}
	if (wrtfname)
		CPXwriteprob(Env, cpx, wrtfname, NULL);
	if (startbas)
		basread(asl, cpx);
#ifdef BARRIER
	if (startvec && (status = CPXreadcopyvec(Env, cpx, startvec)))
		printf("\n*** return %d from CPXreadcopyvec.\n", status);
#endif
#ifdef CPLEX_MIP
	if (nint1 && starttree)
		treeio(cpx, starttree, "read", CPXreadcopytree);
#endif
	fflush(stdout);
	Times[1] = xectim_();
	nosp = 0;

	disconnectchannel(cpxresults);
	addfuncdest(cpxresults, 0, mymsgfunc);
	if (use_netopt) {
		nos = netopt(cpx,&d,&net_status,&net_nodes,&net_arcs,&netiters);
		if (nos)
			net_status = 0;
		nosp = &net_status;
		if (nos != -2)
			status = Optimize1(Env,cpx);
		}
	else
		status = Optimize1(Env,cpx);
	Times[2] = xectim_();
	if (endbas && !nint1)
		write_basis(cpx);
#ifdef CPLEX_MIP
	if (nint1 && endtree)
		treeio(cpx, endtree, "write", CPXtreewrite);
#endif
	amplout(asl, cpx, &d, status, nint1, nosp);
 done:
	if (oic)
		signal(SIGINT, oic);
	if (cpx)
		CPXfreeprob(Env, &cpx);
	if (Env)
		CPXcloseCPLEX(&Env);
	ASL_free(&asl);
	show_times();
	return rc;
	}

 void
#ifdef KR_headers
mainexit_ASL(rc) int rc;
#else
mainexit_ASL(int rc)
#endif
{ longjmp(Jb, rc+1); }
