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

static ASL *asl;

typedef struct cpxlp cpxlp;

 typedef struct
mint_values {
	int L;
	int U;
	int val;
	} mint_values;

 enum { /* sf_mint f values */
	set_crossover	= 0,
	set_dualthresh	= 1,
	set_netopt	= 2,
	set_objno	= 3,
	set_conpri	= 4,
	set_objpri	= 5,
	set_prestats	= 6,
	set_sos2	= 7,
	set_timing	= 8
	};

 static mint_values
mint_val[9] = {
	/* set_crossover */	0, 2, 1,
	/* set_dualthresh */	-0x7fffffff, 0x7fffffff, 32000,
	/* set_netopt */	0, 2, 1,
	/* set_objno */		0, 0/*n_obj*/,	1,
	/* set_conpri */	0, 0x7fffffff, 1,
	/* set_objpri */	0, 0x7fffffff, 2,
	/* set_prestats */	0, 1, 0,
	/* set_sos2 */		0, 1, 1,
	/* set_timing */	0, 3, 0
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

 static char hybmethod = 'o';
 static int netiters = -1;
 static FILE *Logf;
 extern char cplex_bsname[], cplex_version[];
 static char *baralgname, *endbas, *endtree, *endvec, *logfname,
	*startbas, *starttree, *startvec;
 static int dispzap, mbas, method, nbas, netopting, objadj, objsen,
	relax, zap_lpcbf;
 static int net_status, net_nodes, net_arcs, net_its;
 static double obj_adj;
 static unsigned cstorsz, rstorsz;
 typedef struct
Cbfinfo {
	int agg, disp, m, mipdisp, n, n0[4], nint, npr, pres;
	} Cbfinfo;
 static Cbfinfo cbi;

/* Stuff to cope with current CPLEX political correctness... */

 static CPXENVptr Env;

#define disconnectchannel(a)	CPXdisconnectchannel(Env,a)
#define addfuncdest(a,b,c)	CPXaddfuncdest(Env,a,b,c)
#define getitfoind(a)		CPXgetintparam(Env, CPX_PARAM_ITFOIND, a)
#define setitfoind(a,b,c)	CPXsetintparam(Env, CPX_PARAM_ITFOIND, a)
#define netopt(a,b,c,d,e)	CPXnetopt(Env, a,b,c,d,e)

#ifdef KR_headers
 static int (*Optimize)();
#else
 static int (*Optimize)(CPXENVptr,  cpxlp*);
#endif

#ifdef BARRIER

 static int
#ifdef KR_headers
Optimize2(e, c) CPXENVptr e; cpxlp *c;
#else
Optimize2(CPXENVptr e, cpxlp *c)
#endif
{
	return CPXhybbaropt(e, c, hybmethod);
	}

 static int
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
			baralgname = " crossover primal simplex";
			hybmethod = 'o';
			}
		else {
			baralgname = " crossover dual simplex";
			hybmethod = 't';
			}
		}
	else {
		Optimize = Optimizebar;
		baralgname = " barrier";
		}
	}

 static int
#ifdef KR_headers
startcomp(n0, m, nextra, ka, kal, ia, a, b, c, cdualp, cprimp, rdualp, rprimp)
	int n0, m, nextra, *ka, *kal, *ia;
	double *a, *b, *c, **cdualp, **cprimp, **rdualp, **rprimp;
#else
startcomp(int n0, int m, int nextra, int *ka, int *kal, int *ia,
	double *a, double *b, double *c,
	double **cdualp, double **cprimp, double **rdualp, double **rprimp)
#endif
{
	double *r, t, *x, *y;
	int i, j, k;
	int n = n0 - objadj;

	/* Why can't loadstart do this calculation? */

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

 typedef struct {
	int	m;
	int	n;
	int	shift;
	int	nshift;
	double	*x;
	double	*y;
	} dims;

#ifdef CPLEX_MIP

 static int
#ifdef KR_headers
plsos2(m_p, d, nzc_p, nsosnz_p, sostype_p,
	sospri_p, sosbeg_p, sosind_p, sosref_p)
 char **sostype_p;
 dims *d;
 int *m_p, *nzc_p, *nsosnz_p, **sospri_p, **sosbeg_p, **sosind_p;
 double **sosref_p;
#else
plsos2(int *m_p, dims *d, int *nzc_p, int *nsosnz_p, char **sostype_p,
	int **sospri_p, int **sosbeg_p, int **sosind_p, double **sosref_p)
#endif
{
	int i, i1, j, jlim, k, k0, k1, L, L1, mdec, n0, nint, nsos, nsosnz, pr;
	int *ja0, *ja00, *ja1, *sospri, *sosbeg, *sosind, *z;
	char *sostype;
	double *lb, *lb1, *sosref, *ub, *ub1;
	double *a = A_vals;
	int *ia = A_rownos;
	int *ja = A_colstarts;
	int m0 = *m_p;
	int n = d->n;
	double *clb = LUrhs;
	double *cub = Urhsx;

	ja00 = ja;
	ja0 = ja += n0 = ampl_options[3] - 1;
	nint = niv + nbv;
	d->shift = jlim = n - n0 - nint;
	k0 = ampl_options[4];

	for(nsos = 0;;) {
		k = ia[i = *ja];
		if (k < k0)
			k = ia[++i];
		for(;;) {
			if ((i1 = ia[i = *++ja]) < k0)
				i1 = ia[++i];
			if (i1 != k)
				break;
			}
		nsos++;
		if (ja - ja0 >= jlim)
			break;
		}
	nsosnz = *nsosnz_p = ja - ja0;
	sosref = *sosref_p = (double *)Malloc(nsos
				+ (2*nsos + 1)*sizeof(int)
				+ nsosnz*(sizeof(int) + sizeof(double)));
	sosind = *sosind_p = (int *)(sosref + nsosnz);
	sosbeg = *sosbeg_p = sosind + nsosnz;
	sospri = *sospri_p = sosbeg + nsos + 1;	/* + 1 for safety */
	sostype = *sostype_p = (char *)(sospri + nsos);
	z = (int *)Malloc(m0*sizeof(int));
	for(L = 0; L < k0; L++)
		z[L] = L;
	L1 = L;
	*sosbeg = mdec = 0;
	ja = ja0;
	do {
		k = ia[i = *ja];
		pr = objpri;
		if (k < k0) {
			k = ia[++i];
			pr = conpri;
			}
		j = *++ja;
		if (ia[j] < k0) {
			j++;
			pr = conpri;
			}
		ja1 = ja;
		k1 = k + 1;
		*sosref++ = ia[++i] == k1 ? -a[i] : 0.;
		*sosref++ = ia[++j] == k1 ? -a[j] : 0.;
		for(;;) {
			if ((i1 = ia[i = *++ja]) < k0) {
				i1 = ia[++i];
				pr = conpri;
				}
			if (i1 != k)
				break;
			*sosref++ = ia[++i] == k1 ? -a[i] : 0.;
			}
		z[L++] = L1++;
		z[L++] = L1++;
		k = ja - ja1;
		mdec += k += k + 1;
		while(--k >= 0)	/* omit introduced binary variables */
			z[L++] = -1;
		*sospri++ = pr;
		*sostype++ = '2';
		*++sosbeg = ja - ja0;
		}
		while (ja - ja0 < jlim);
	*m_p = m0 - mdec;

	L1 = 0;
	ja1 = ja00;
	i = *ja1;
	ja += nbv - jlim;
	while(ja1 < ja) {
		i = *ja1;
		*ja1++ = L1;
		for(j = *ja1; i < j; i++)
			if ((k = z[ia[i]]) >= 0) {
				ia[L1] = k;
				a[L1++] = a[i];
				}
		}

	/* copy integer variables */

	j = ja1 - ja00;
	lb = LUv + j;
	lb1 = lb + jlim;
	ub = Uvx + j;
	ub1 = ub + jlim;
	ja += jlim;
	i = *ja;
	ja0 = ja00 + n;
	d->nshift = ja0 - ja;
	while(ja < ja0) {
		*lb++ = *lb1++;
		*ub++ = *ub1++;
		*ja1++ = L1;
		for(j = *++ja; i < j; i++)
			if ((k = z[ia[i]]) >= 0) {
				ia[L1] = k;
				a[L1++] = a[i];
				}
		}
	*ja1 = *nzc_p = L1;
	for(i = 0; i < L; i++)
		if ((j = z[i]) >= 0) {
			clb[j] = clb[i];
			cub[j] = cub[i];
			}
	free((char *)z);
	while(--nsosnz >= 0)
		*sosind++ = n0++;
	return nsos;
	}

 static void
#ifdef KR_headers
treeio(cpx, fname, what, tio) cpxlp *cpx; char *fname, *what; int (*tio)();
#else
treeio(cpxlp *cpx, char *fname, char *what, int (*tio)(CPXENVptr,cpxlp*,char*))
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
	set_version	= 8
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
		Optimize = CPXoptimize;
		baralgname = "";
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
		printf(Fmt, oi->nnl ? "\n" : "", cplex_version+1);
		oi->option_echo &= ~ASL_OI_echothis;
#undef Fmt
		if (oi->option_echo & ASL_OI_clopt)
			exit(0);
		break;
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

 static char **file_name[7] = { &endbas, &endtree, &startbas, &starttree,
				&endvec, &startvec, &logfname };

enum {	/* sf_char f values */
	set_endbas	= 0,
	set_endtree	= 1,
	set_startbas	= 2,
	set_starttree	= 3,
	set_endvector	= 4,
	set_startvector	= 5,
	set_logname	= 6
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
		if (!(Logf = fopen(t,"w"))) {
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
	{ "barobjrange",sf_dbl,		VP CPX_PARAM_BAROBJRNG },
	{ "baropt",	sf_known,	VP set_barrier },
	{ "barthreads",	sf_int,		VP CPX_PARAM_BARTHREADS },
	{ "barvarup",sf_dbl,		VP CPX_PARAM_BARVARUP },
#endif
	{ "basisinterval",sf_int,	VP CPX_PARAM_BASINTERVAL },
#ifdef CPLEX_MIP
	{ "branch",	sf_int,		VP CPX_PARAM_BRDIR },
	{ "branchdir",	sf_int,		VP CPX_PARAM_BRDIR },
	{ "cliques",	sf_int,		VP CPX_PARAM_CLIQUES },
#endif
	{ "clocktype",	sf_int,		VP CPX_PARAM_CLOCKTYPE },
#ifdef CPLEX_MIP
	{ "coeffreduce",sf_int1,	VP CPX_PARAM_COEREDIND },
#endif
#ifdef BARRIER
	{ "comptol",	sf_dbl,		VP CPX_PARAM_BAREPCOMP },
#endif
#ifdef CPLEX_MIP
	{ "covers",	sf_int,		VP CPX_PARAM_COVERS },
#endif
	{ "crash",	sf_int,		VP CPX_PARAM_CRAIND },
#ifdef BARRIER
	{ "crossover",	sf_mint,	VP set_crossover },
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
	{ "feasibility",sf_dbl,		VP CPX_PARAM_EPRHS },
#ifdef BARRIER
	{ "growth",	sf_dbl,		VP CPX_PARAM_BARGROWTH }, /*== bargrowth*/
#endif
#ifdef CPLEX_MIP
	{ "heuristic",	sf_int,		VP CPX_PARAM_HEURISTIC },
	{ "integrality",sf_dbl,		VP CPX_PARAM_EPINT },
#endif
	{ "iterations",	sf_int,		VP CPX_PARAM_ITLIM },
	{ "iterlim",	sf_int,		VP CPX_PARAM_ITLIM },
	{ "logfile",	sf_char,	VP set_logname },
#ifdef CPLEX_MIP
	{ "lowercutoff",sf_dbl,		VP CPX_PARAM_CUTLO },
#endif
	{ "lowerobj",	sf_dbl,		VP CPX_PARAM_OBJLLIM },
	{ "lowerobjlim",sf_dbl,		VP CPX_PARAM_OBJLLIM },
	{ "lpdisplay",  sf_int,		VP CPX_PARAM_SIMDISPLAY },
	{ "lpiterlim",  sf_int,		VP CPX_PARAM_ITLIM },
	{ "lptimelim",	sf_dbl,		VP CPX_PARAM_TILIM },
	{ "markowitz",	sf_dbl,		VP CPX_PARAM_EPMRK },
	{ "maximize",	sf_known,	VP set_max },
	{ "minimize",	sf_known,	VP set_min },
#ifdef CPLEX_MIP
	{ "mipalg",	sf_int,		VP CPX_PARAM_SUBALG },
	{ "mipalgorithm",sf_int,	VP CPX_PARAM_SUBALG },
	{ "mipdisplay",	sf_int,		VP CPX_PARAM_MIPDISPLAY },
	{ "mipgap",	sf_dbl,		VP CPX_PARAM_EPGAP },
	{ "mipinterval",sf_int,		VP CPX_PARAM_MIPINTERVAL },
	{ "mipsolutions", sf_int,	VP CPX_PARAM_INTSOLLIM },
	{ "mipstartalg",sf_int,		VP CPX_PARAM_STARTALG },
	{ "mipsubalg",	sf_int,		VP CPX_PARAM_SUBALG },
	{ "mipthreads",	sf_int,		VP CPX_PARAM_MIPTHREADS },
#endif
	{ "netfind",	sf_int,		VP CPX_PARAM_NETFIND },
	{ "netfinder",	sf_int,		VP CPX_PARAM_NETFIND },
	{ "netopt",	sf_mint,	VP set_netopt },
#ifdef CPLEX_MIP
	{ "node",	sf_int,		VP CPX_PARAM_NDLIM },
	{ "nodefile",	sf_int,		VP CPX_PARAM_NDFILIND },
	{ "nodelim",	sf_int,		VP CPX_PARAM_NDLIM },
	{ "nodes",	sf_int,		VP CPX_PARAM_NDLIM },
	{ "nodesel",	sf_int,		VP CPX_PARAM_NDSEL },
	{ "nodeselect",	sf_int,		VP CPX_PARAM_NDSEL },
	{ "objdifference", sf_dbl,	VP CPX_PARAM_OBJDIF },
#endif
	{ "objno",	sf_mint,	VP set_objno },
	{ "optimality",	sf_dbl,		VP CPX_PARAM_EPOPT },
	{ "optimize",	sf_known,	VP set_primalopt },
#ifdef BARRIER
	{ "ordering",	sf_int,		VP CPX_PARAM_BARORDER },
#endif
	{ "outlev",	sf_int,		VP CPX_PARAM_SIMDISPLAY },
	{ "perturb",	sf_int1,	VP CPX_PARAM_PERIND },
	{ "perturbation", sf_dbl,	VP CPX_PARAM_EPPER },
	{ "perturbconst", sf_dbl,	VP CPX_PARAM_EPPER },
	{ "pgradient",	sf_int,		VP CPX_PARAM_PPRIIND },
#ifdef CPLEX_MIP
	{ "plconpri",	sf_mint,	VP set_conpri },
	{ "plobjpri",	sf_mint,	VP set_objpri },
#endif
	{ "presolve",	sf_int1,	VP CPX_PARAM_PREIND },
	{ "prestats",	sf_mint,	VP set_prestats },
	{ "pricing",	sf_int,		VP CPX_PARAM_PRICELIM },
	{ "primal",	sf_known,	VP set_primal },
	{ "primalopt",	sf_known,	VP set_primalopt },
#ifdef CPLEX_MIP
	{ "priorities",	sf_int1,	VP CPX_PARAM_ORDIND },
#endif
	{ "readbasis",	sf_char,	VP set_startbas },
#ifdef BARRIER
	{ "readvector",sf_char,		VP set_startvector },
#endif
#ifdef CPLEX_MIP
	{ "reducecostfix",sf_int1,	VP CPX_PARAM_RCFIXIND },
#endif
	{ "refactor",	sf_int,		VP CPX_PARAM_REINV },
#ifdef CPLEX_MIP
	{ "relax",	sf_known,	VP set_relax },
	{ "relobjdif",	sf_dbl,		VP CPX_PARAM_RELOBJDIF },
	{ "relobjdiff",	sf_dbl,		VP CPX_PARAM_RELOBJDIF },
#endif
	{ "scale",	sf_int,		VP CPX_PARAM_SCAIND },
	{ "simthreads",	sf_int,		VP CPX_PARAM_SIMTHREADS },
	{ "singular",	sf_int,		VP CPX_PARAM_SINGLIM },
	{ "singularlim",sf_int,		VP CPX_PARAM_SINGLIM },
#ifdef CPLEX_MIP
	{ "solutionlim",sf_int,		VP CPX_PARAM_INTSOLLIM },
	{ "sos2",	sf_mint,	VP set_sos2 },
	{ "sosmin",	sf_int,		VP CPX_PARAM_SOSMINSZ },
	{ "sosscan",	sf_int1,	VP CPX_PARAM_SOSIND },
	{ "startalg",	sf_int,		VP CPX_PARAM_STARTALG },
	{ "startalgorithm",sf_int,	VP CPX_PARAM_STARTALG },
	{ "startbasis",	sf_char,	VP set_startbas },
	{ "starttree",	sf_char,	VP set_starttree },
#ifdef BARRIER
	{ "startvector",sf_char,	VP set_startvector },
#endif
	{ "subalg",	sf_int,		VP CPX_PARAM_SUBALG },
	{ "subalgorithm",sf_int,	VP CPX_PARAM_SUBALG },
#endif
	{ "time",	sf_dbl,		VP CPX_PARAM_TILIM },
	{ "timelimit",	sf_dbl,		VP CPX_PARAM_TILIM },
	{ "timing",	sf_mint, 	VP set_timing },
	{ "tranopt",	sf_known,	VP set_dualopt },
#ifdef CPLEX_MIP
	{ "treelimit",	sf_dbl,		VP CPX_PARAM_TRELIM },
	{ "treememlim",	sf_dbl,		VP CPX_PARAM_TRELIM },
	{ "treememory",	sf_dbl,		VP CPX_PARAM_TRELIM },
	{ "uppercutoff",sf_dbl,		VP CPX_PARAM_CUTUP },
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
#ifdef BARRIER
	{ "writevector", sf_char,	VP set_endvector },
#endif
	{ "xxxstart",	sf_int,		VP CPX_PARAM_XXXIND }
	};

 static Option_Info Oinfo = { "cplex", cplex_bsname, "cplex_options",
				keywds, nkeywds, 0, cplex_version+1};

 static void
#ifdef KR_headers
nonlin(n, what) int n; char *what;
#else
nonlin(int n, char *what)
#endif
{
	if (n) {
		printf("%s contains %s.\n", filename, what);
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

 static char **cname, **rname;

 static void
#ifdef KR_headers
make_names()
#else
make_names(void)
#endif
{
	char *s;

	cstorsz = namemem(nbas);
	rstorsz = namemem(mbas);
	cname = (char **)M1alloc((mbas+nbas)*sizeof(char *)
				+ cstorsz + rstorsz);
	rname = cname + nbas;
	s = namegen("x", cname, (char *)(rname + mbas), nbas);
	namegen("c", rname, s, mbas);
	}

#ifdef CPLEX_MIP
 static void
#ifdef KR_headers
mip_priorities(cpx, d) cpxlp *cpx; dims *d;
#else
mip_priorities(cpxlp *cpx, dims *d)
#endif
{
	int i, k, listsize, nnames, pvk, shlim, sk;
	int *ci, *colindex, *p, *priority, *start, *num, *pri;

	if (nnames = mip_pri(&start, &num, &pri, 2147483647)) {
		for(i = listsize = 0; i < nnames; i++)
			listsize += num[i];
		ci = colindex = (int *)M1alloc(listsize*(2*sizeof(int)));
		priority = p = ci + listsize;
		shlim = d->n - d->nshift;
		for(k = 0; k < nnames; k++) {
			i = num[k];
			pvk = pri[k];
			if ((sk = start[k]) >= shlim)
				sk -= d->shift;
			while(--i >= 0) {
				*ci++ = sk++;
				*p++ = pvk;
				}
			}
		if (CPXloadorder(Env,cpx,listsize,colindex,priority, (int *)0))
			printf("error in loadorder!\n");
		}
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

 static int
#ifdef KR_headers
lpcbf(env, lp, wf) CPXENVptr env; CPXLPptr lp; int wf;
#else
lpcbf(CPXENVptr env, CPXLPptr lp, int wf)
#endif
{
	int i;
	static int n[4];
	static int wi[4] = {
		CPX_CBINFO_PRESOLVE_ROWSGONE,
		CPX_CBINFO_PRESOLVE_COLSGONE,
		CPX_CBINFO_PRESOLVE_AGGSUBST,
		CPX_CBINFO_PRESOLVE_COEFFS };
	static char pfmt[] =
		"%sLP Presolve eliminated %d rows and %d columns.\n";
	static int wf0;
	char *pfx, *prob;

	if (wf == CPX_CALLBACK_PRESOLVE) {
		for(i = 0; i < 4; i++)
			CPXgetcallbackinfo(env, lp, wf, wi[i], &n[i]);
		wf0 = wf;
		goto done;
		}
	if (wf == wf0)
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
	return 0;
	}

 static cpxlp *
#ifdef KR_headers
amplin(stub, d, nintp, av)
	char *stub; dims *d; int *nintp; char **av;
#else
amplin(char *stub, dims *d, int *nintp, char **av)
#endif
{
	int *ia, *ia1, *ia2, *ja, *ka, *ka1, *ka2, *kal, *kal1, *kal2;
	double *a, *a1, *b, *b1, *c, *c1, *l, *l1,
		*le, *lrhs, *lx, *rngvec, *rv1, *u, *u1, *urhs, *ux;
	double nos, os;
	ograd *og, *og1;
	int i, j, k, m, m0, n, n0, nbv1, nextra, nint, nint0, nint1,
		nr, nsos, nz, nz1, nzc0, nzextra, nzr;
	char *s, *senx;
	Char **LUvmem;
#ifdef BARRIER
	fint nelqf, *colqf, *rowqf;
	int havestart, nelq, *colq, *colqcnt, *rowq;
	double *cdual, *cprim, *qmat, *rdual, *rprim;
#endif
	static char probname[48];
	FILE *nl;
	char *cstore, *rstore;
	cpxlp *rv;
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
	Optimize = CPXoptimize;
	baralgname = "";

	nl = jac0dim(stub, (fint)strlen(stub));

	if (!(mint_val[set_objno].U = n_obj))
		objno = 0;
	d->m = m = m0 = n_con;
	d->n = n = n0 = nr = n_var;
	d->nshift = d->shift = 0;
	nzc0 = nzc;
	if (n0 <= 0) {
		fprintf(Stderr, "%s has no variables.\n", filename);
		exit(4);
		}
	nonlin(nlc, "nonlinear constraints");
	nonlin(nlvoi, "nonlinear integer variables");
	nonlin(plterms, "piecewise-linear terms");

	if (getopts(av, &Oinfo))
		exit(1);
	obj_no = objno - 1;

	nint = nint0 = nbv + niv;
	nsos = 0;
	if (ampl_options[3] && (i = n0 - nint - ampl_options[3] + 1) > 0)
		nsos = sos2;
	if (nint
#ifdef CPLEX_MIP
		&& relax
#endif
			  ) {
		if (nint1 = nint - ((nint0 = nsos) ? i : 0)) {
			printf("Ignoring integrality of %d variable%s.\n",
				nint1, nint1 > 1 ? "s" : "");
			need_nl = nint = 0;
			}
		}
	*nintp = nint0;	/* would make nint file static,
			   but the inane Sun math.h declares
			   a mysterious nint()... */
	if (use_netopt == 1 && !lnc)
		use_netopt = 0;
	if (nint0 || m <= 0) {
		method = 1;
		use_netopt = 0;
		}
	else if (!method)
		method = m - n > dual_thresh ? -1 : 1;

	nzr = nzc0;
	if (method > 0 || nlo) {
		/* allow room for adjusting the objective */
		m++;
		nr = n + nranges + 1;
		nzr += nranges + 1;
		}
	if (method < 0 && use_netopt == 1)
		use_netopt = 0;
	b = LUrhs = (double *)M1alloc(2*Sizeof(double)*m);
	Urhsx = b + m;
	LUv = (double *)Malloc((2L*nr + nzr)*Sizeof(double)
				+ ((fint)nzr + nr + 1)*sizeof(int));
	LUvmem = M1record(LUv);
	Uvx = LUv + nr;
	a = A_vals = Uvx + nr;
	ia = A_rownos = (int *)(a + nzr);
	ka = A_colstarts = ia + nzr;
#ifdef BARRIER
	if (Optimize == Optimize2 || Optimize == Optimizebar || nlo)
		want_xpi0 = 3;
	want_deriv = 0;
	qp_read(nl,0);
	if (nelqf = qpcheck(&rowqf, &colqf, &qmat)) {
		nonlin(nint, "integer variables and a quadratic objective");
		crossover = 0;
		set_baropt();
		method = 1;
		}
#else
	ed0read(nl);
#endif

	nbv1 = nbv;
#ifdef CPLEX_MIP
	if (nsos) {
		nsos = plsos2(&m0, d, &nzc0, &nsosnz, &sostype,
			&sospri, &sosbeg, &sosind, &sosref);
		n0 -= d->shift;
		nbv1 -= d->shift;
		}
#endif
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
	nonlin(nl_obj(obj_no), "nonlinear objective terms");
#endif
	*stub_end = 0;
	for(i = 0, og1 = og; og1; i++)
		og1 = og1->next;
	sprintf(probname, "c%dv%di%do%d", n_var, n_con, niv, i);
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
			free((char *)ja);
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
		nzr = nz1 = nz += nzextra;
		a = (double *)M1alloc((nz1+2L*n+m)*sizeof(double)
					+ (nz1+1L+n)*sizeof(int));
		l = a + nz1;
		u = l + n;
		b = u + n;
		ia = (int *)(b + m);
		ka = ia + nz1;
		nr = n;
		}
	else {
		m += objadj;
		nz1 = nz;
		l = LUv;
		u = Uvx;
		if (nint0) {
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
			fprintf(Stderr,
		"CPLEX's MIP option is needed to solve problems with\n%s.\n",
			"nonconvex (nonconcave) piecewise-linear terms");
			exit(1);
#endif
			}
		}

	d->x = (double *)M1alloc((m+2L*nr)*sizeof(double));
	d->y = d->x + nr;
	c = d->y + m;
	/* return zeros for infeasible problems... */
	memset((char *)d->x, 0, (m+2*nr)*sizeof(double));

	rngvec = (double *)M1alloc(m*(Sizeof(double) + 1) + nr*Sizeof(int));
	kal = (int *)(rngvec + m);
	senx = (char *)(kal + nr);
#ifdef CPLEX_MIP
	if (nint0) {
		Optimize = CPXmipoptimize;
		baralgname = "";
		}
#endif
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
		free((char *)LUv);
		*LUvmem = 0;
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
		i = d->shift;
		j = d->n - d->nshift;
		for(; og; og = og->next) {
			if ((k = og->varno) >= j)
				k -= i;
			c[k] = og->coef;
			}
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
	nbas = n; /* was method > 0 ? n + nranges : n; */
	cstore = rstore = 0;
	if (startbas || endbas || startvec || endvec
	 || Optimize == CPXoptimize || Optimize == CPXdualopt
	 || parval(CPX_PARAM_XXXIND)) {
		make_names();
		cstore = *cname;
		rstore = *rname;
		}
#ifdef CPLEX_MIP
	if (nint0) {
		rv = CPXloadmprob(Env,probname,n,m,0,objsen,c,b,senx, /*8*/
		 ka, kal, ia, a, l, u, rngvec, 0 /*norowind*/, 0,0,0,0,0, /*21*/
		 0 /*dataname*/, 0,0,0,0, cname, cstore, rname, rstore,
		 0,0, nr, m, nzr, 0, 0, cstorsz, rstorsz, 0, ctype);
#ifdef CPLEX_MIP
		if (nsos) {
			if (CPXloadsos(Env,rv,nsos,nsosnz,sostype,sospri,
					sosbeg,sosind,sosref)) {
				printf("Failure of loadsos: ");
				return 0;
				}
			free((char *)sosref);
			}
#endif
		if (nint)
			mip_priorities(rv, d);
		}
	else
#endif
		{
#ifdef BARRIER
		havestart = 0;
		if (Optimize == Optimize2 || Optimize == Optimizebar)
			havestart = startcomp(n, m, nextra, ka, kal, ia,
				a, b, c, &cdual, &cprim, &rdual, &rprim);
#endif
		rv = CPXloadprob(Env,probname,n,m,0,objsen,c,b,senx, /*8*/
		 ka, kal, ia, a, l, u, rngvec, 0 /*norowind*/, 0,0,0,0,0, /*21*/
		 0 /*dataname*/, 0,0,0,0, cname, cstore, rname, rstore,
		 0,0, nr, m, nzr, 0, 0, cstorsz, rstorsz, 0);
#ifdef BARRIER
		if (nelqf && rv) {
			nelq = (int)nelqf;
			colq = (int*)M1alloc(2*sizeof(int)*nr);
			colqcnt = colq + nr;
			if (sizeof(fint) != sizeof(int)) {
				/* should not happen */
				rowq = (int*)M1alloc(nelq*sizeof(int));
				for(i = 0; i < nelq; i++)
					rowq[i] = (int)rowqf[i];
				}
			else
				rowq = (int*)rowqf;
			k = n_var;
			for(i = j = 0; i < k; i++) {
				colq[i] = (int)colqf[i];
				colqcnt[i] = (int)(colqf[i+1] - colqf[i]);
				}
			for(;i < nr; i++)
				colq[i] = colqcnt[i] = 0;
			for(i = 0; ; i++) {
				if (i == k) {
					qmat = (double *)Realloc(qmat,
							nr*sizeof(double));
					while(i < nr)
						qmat[i++] = 0.;
					if (CPXloadqsep(Env,rv,qmat)) {
						printf("loadqsep failure!\n");
						exit(4);
						}
					baralgname = " separable QP barrier";
					break;
					}
				if (!(j = colqcnt[i]))
					continue;
				if (j > 1 || rowq[colq[i]] != i) {
					if (CPXloadquad(Env,rv,colq,colqcnt,
							rowq,qmat,nelq)) {
						printf("loadquad failure!\n");
						exit(4);
						}
					baralgname = " QP barrier";
					break;
					}
				}
			}
		if (havestart) {
			if (i = CPXloadstart(Env,rv,0,0,cprim,rprim,cdual,rdual))
				fprintf(Stderr,
					"Return %d from CPXloadstart.\n", i);
			}
#endif
		}
	cbi.disp = parval(CPX_PARAM_SIMDISPLAY);
#ifdef CPLEX_MIP
	cbi.mipdisp = parval(CPX_PARAM_MIPDISPLAY);
#endif
	if (prestats && (dispzap = !cbi.disp || nint && !cbi.mipdisp)) {
		zap_lpcbf = 1;
		cbi.agg = parval(CPX_PARAM_AGGIND);
		cbi.m = m;
		cbi.n = n;
		cbi.nint = nint;
		cbi.pres = parval(CPX_PARAM_PREIND);
		CPXsetlpcallbackfunc(Env, lpcbf);
		}
	return rv;
	}

 static void
#ifdef KR_headers
round(x, n) double *x; fint n;
#else
round(double *x, fint n)
#endif
{
	double *xe;
	for(xe = x + n; x < xe; x++)
		*x = floor(*x + 0.5);
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
	  case CPXERR_NO_ROWS:
		return "no rows in problem";
	  case CPXERR_NO_COLS:
		return "no columns in problem";
	  case CPXERR_NOT_LP:
		return "not a linear program";
	  case CPXERR_SB_INCOMPAT:
		return "primal values exist";
	  case CPXERR_PRIIND:
		return "invalid pricing setting";
	  case CPXERR_BOUNDS_INFEAS:
		return "contradictory bounds";
	  case CPXERR_PRESLV_INForUNBD:
		return "presolve: infeasible or unbounded problem";
	  case CPXERR_INDEX_RANGE:
		return "illegal list length";
	  case CPXERR_NET_SMALL:
		return "no network";
	  case CPXERR_NET_IMBALANCE:
		return "an infeasible or unbounded network";
	  case CPXERR_BAD_METHOD:
		return "unknown method after network optimization";
	  case CPXERR_LOCK_TAMPER:
		return "licensing problem during optimization";
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
	if (stat > 100 && stat < 117)
		stat -= 71;
	else if (stat > 18 && stat < 42)
		stat -= 13;
	else if (stat == 1101)
		stat = 29;
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
amplout(cpx, d, status, nint, nosp) cpxlp *cpx; dims *d; int status, nint, *nosp;
#else
amplout(cpxlp *cpx, dims *d, int status, int nint, int *nosp)
#endif
{
	char buf[32], hbuf[256];
	double obj, *l, *le, *u, *x, *x1, *x2, *y, *z, *z1;
	int ccpr, cl0, cl, cov, i, itc, itci, m, mipstat, n, nodelim, nos,
		stat, stat0, stat10;
	int ii = 0, stat1 = 1;
#ifdef CPLEX_MIP
	int nodecnt;
	double obj1;
#endif
	static char *statmsg[] = {
		"unrecoverable failure",
		"optimal solution",
		"infeasible problem",		/* 2 */
		"unbounded problem",
		"phase II objective limit exceeded",		/* 4 */
		"phase II iteration limit",
		"phase I time limit",				/* 6 */
		"phase II time limit",
		"phase I time limit",				/* 8 */
		"infeasible with phase II singularities",
		"infeasible with phase I singularities",
		"optimal with unscaled infeasibilities",	/* 11 */
		"aborted in phase II",
		"aborted in phase I",
		"aborted in barrier, dual infeasible",		/* 14* */
		"aborted in barrier, primal infeasible",
		"aborted in barrier, primal and dual infeasible",
		"aborted in barrier, primal and dual feasible",
		"aborted in crossover",				/* 18 */

		"converged, dual feasible, primal infeasible",	/*32 */
		"converged, primal feasible, dual infeasible",
		"converged, primal and dual infeasible",
		"primal objective limit reached",		/* 35 */
		"dual objective limit reached",
		"primal has unbounded optimal face",
		"best solution found, primal-dual feasible",	/* 38 */
		"best solution found, primal infeasible",
		"best solution found, dual infeasible",
		"best solution found, primal-dual infeasible",	/* 41 */

		"infeasible or unbounded in presolve",		/*1101*/
#ifdef CPLEX_MIP
		"optimal integer solution",			/* 101 */
		"optimal integer solution within mipgap or absmipgap",
		"integer infeasible",				/* 103 */
		"mixed-integer solutions limit",
		"node limit with integer solution",		/* 105 */
		"node limit with no integer solution",
		"time limit with integer solution",
		"time limit with no integer solution",		/* 108 */
		"unrecoverable failure with integer solution",
		"unrecoverable failure with no integer solution",
		"treememory limit with integer solution",
		"treememory limit with no integer solution",	/* 112 */
		"aborted, integer solution exists",
		"aborted, no integer solution",			/* 114 */
		"integer optimal with unscaled infeasibilities",
		"out of memory, no tree; solution may exist"	/* 116 */
#endif
		};
	static char *netmsg[] = {
		"memory failure",
		"an infeasible or unbounded network",
		"an unbounded network",
		"an infeasible network",
		/*"an optimal network solution"*/ "",
		"no network",
		"a bug" };
	static char wantobj[] = {0,1,0,0,1,1,1,1,1,0,0,1,	/*  0-11 */
				 0,0,0,0,0,0,0,			/* 12-18 */
				 1,1,1,1,1,0,1,1,1,1,		/* 32-41 */
				 0,				/* 1101 */
				 1,1,0,1,1,0,1,0,1,0,1,0,
				 1,0,1,1};

	mipstat = 0;
#ifdef CPLEX_MIP
	ccpr = 0;
	if (Optimize == CPXmipoptimize) {
		if (prestats && nint
		 && (parval(CPX_PARAM_CLIQUES) || parval(CPX_PARAM_COVERS))) {
			cl0 = CPXgetgclqc(Env, cpx);
			cl  =  CPXgetclqc(Env, cpx);
			cov =  CPXgetcovc(Env, cpx);
			ccpr = 1;
			}
		if (status > 0) {
			mipstat = statadjust(stat0 = status);
			status = 0;
			}
		}
#endif
	m = d->m;
	n = d->n;
	x = d->x;
	y = d->y;
	itc = CPXgetitc(Env, cpx);
	itci = CPXgetitci(Env, cpx);
	nos = nosp ? *nosp + 5 : 0;
	if (status) {
		stat0 = status;
		stat = 0;
		x = y = 0;
		goto have_stat;
		}
#ifdef CPLEX_MIP
	if (nint) {
		if (mipstat) {
			stat = mipstat;
			if (mipstat < 4 || mipstat > 8) {
				ii = 0;
				goto have_ii;
				}
#if 0
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
#endif
			}
		else
			stat = statadjust(stat0 = CPXgetstat(Env, cpx));
		ii = CPXgetitcm(Env, cpx);
		nodecnt = CPXgetndc(Env, cpx);
 have_ii:
		if (wantobj[stat] || endbas) {
			if (zap_lpcbf)
				CPXsetlpcallbackfunc(Env, NULL);
			CPXgetmobjval(Env, cpx, &obj);
			CPXchgprobtype(Env, cpx, 3);
			status = CPXoptimize(Env, cpx);
			if (status)
				fprintf(Stderr,
				"CPLEX status %d with fixed integers:\n\t%s\n",
					status, statmsg[statadjust(status)]);
			if (endbas)
				write_basis(cpx);
			if (!wantobj[stat])
				goto have_stat;
			CPXsolution(Env, cpx, &stat1, &obj1, x, y, 0, 0);
			stat1 = statadjust(stat10 = stat1);
			if (i = d->nshift) {
				x2 = x + n;
				x1 = x2 - d->shift;
				while(--i >= 0)
					*--x2 = *--x1;
				for(i = d->shift; --i >= 0; )
					*--x2 = 0;
				}
			/* round to integer */
			if (!relax) {
				x2 = x + n_var - nint;
				if (d->nshift) {
					round(x2, nbv - d->nshift);
					round(x2 + nbv, niv);
					}
				else
					round(x2, nint);
				}
			}
		}
	else
#endif
		{
		CPXsolution(Env, cpx, &stat, &obj, x, y, 0, 0);
		stat = nos == 3 ? 2 : statadjust(stat0 = stat);
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
	i = Sprintf(hbuf, "%s: %s", cplex_bsname, statmsg[stat]);
	if (!stat)
		i += Sprintf(hbuf+i, ": %s", failstat(stat0));
#ifdef CPLEX_MIP
	else if (stat == 16 || stat == 17) {
		CPXgetintparam(Env, CPX_PARAM_NDLIM, &nodelim);
		i += Sprintf(hbuf+i, ".\nCurrent node limit = %d", nodelim);
		}
#endif
	if (wantobj[stat]) {
		g_fmtop(buf, obj);
		i += Sprintf(hbuf+i, "; objective %s", buf);
		}
	if (nosp) {
		if (nos < 0 || nos > 6)
			nos = 6;
		i += nos == 4
			? Sprintf(hbuf+i, "\n%d network simplex iterations.",
				netiters)
			: Sprintf(hbuf+i, "\nnetopt found %s.", netmsg[nos]);
		}
#ifdef CPLEX_MIP
	if (nint && !mipstat)
		i += Sprintf(hbuf+i,
			"\n%d simplex iterations\n%d branch-and-bound nodes",
			ii, nodecnt);
	else
#endif
	     if (netiters < 0 || itc) {
		i += Sprintf(hbuf+i, "\n%d%s iterations", itc,
				ii ? " simplex" : baralgname);
		if (itc)
#ifdef BARRIER
		if (!*baralgname)
#endif
			i += Sprintf(hbuf+i, " (%d in phase I)", itci);
		if (ii)
			Sprintf(hbuf+i, "\n%d integer iterations", ii);
		}
#ifdef CPLEX_MIP
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
			 statmsg[stat1]);
		if (!stat1)
			i += Sprintf(hbuf+i, ": stat = %d", stat10);
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
basread(cpx) cpxlp *cpx;
#else
basread(cpxlp *cpx)
#endif
{
	int *cstat, *rstat;

	cstat = (int *)M1alloc((nbas+mbas)*sizeof(int));
	rstat = cstat + nbas;
	if (CPXmbaseread(Env,startbas,nbas,mbas,cname,rname,cstat,rstat)) {
		printf("Could not read starting basis file \"%s\"\n",
			startbas);
		need_nl = 0;
		}
	else if (CPXloadbase(Env, cpx, cstat, rstat)) {
		printf("Failed to load basis!\n");
		need_nl = 0;
		}
	}

 static void
#ifdef KR_headers
mymsgfunc(handle, msg) char *handle, *msg;
#else
mymsgfunc(void *handle, char *msg)
#endif
{
	static int wantnl;
	char *msg1;

	handle = handle; /* shut up non-use warning */
	if (Logf)
		fprintf(Logf, msg);
	if (dispzap) {
		if (*msg != '\n')
			return;
		dispzap = 0;
		}
	if (!strncmp(msg, "Total network iterations = ", 27))
		netiters = atoi(msg+27);
	else if (prestats) {
		if (!netopting || !strncmp(msg,"Extract",7))
			printf("%s", msg);
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
			else
				printf("%s", msg);
		  }
		}
	}

 static jmp_buf Jb;

 int
#ifdef KR_headers
main(argc, argv) int argc; char **argv;
#else
main(int argc, char **argv)
#endif
{
	int	itfoind, nint, nos, *nosp, rc, status, z;
	cpxlp	*cpx;
	dims	d;
	char	*stub;
	CPXCHANNELptr cpxresults;

	Times[0] = xectim_();

	Env = 0;
	cpx = 0;
	rc = setjmp(Jb);
	if (rc) {
		--rc;
		goto done;
		}
	if (argc < 2)
#ifdef Stand_alone
		return cpxmain(argc,argv);
#else
		usage_ASL(&Oinfo, 1);
#endif
	if (!(Env = CPXopenCPLEX(&z))) {
		fprintf(Stderr, "openCPLEX failed: error code %d\n", z);
		return 1;
		}
	if (CPXgetchannels(Env, &cpxresults, 0,0,0)) {
		printf("CPXgetchannels failure!\n");
		rc = 1;
		goto done;
		}

#ifdef BARRIER
	asl = ASL_alloc(ASL_read_fg);
#else
	asl = ASL_alloc(ASL_read_f);
#endif
	if (!(stub = getstub(&argv, &Oinfo)))
		usage_ASL(&Oinfo, 1);
	cpx = amplin(stub, &d, &nint, argv);
	if (!cpx) {
		printf("Failed to load the problem!\n");
		return 1;
		}
	if (startbas)
		basread(cpx);
#ifdef BARRIER
	if (startvec && (status = CPXvecread(Env, cpx, startvec)))
		printf("\n*** return %d from CPXvecread.\n", status);
#endif
#ifdef CPLEX_MIP
	if (nint && starttree)
		treeio(cpx, starttree, "read", CPXtreeread);
#endif
	fflush(stdout);
	Times[1] = xectim_();
	nosp = 0;

	disconnectchannel(cpxresults);
	addfuncdest(cpxresults, 0, mymsgfunc);
	if (use_netopt) {
		getitfoind(&itfoind);
		if (!itfoind) {
			netopting = 1; /* suppress netopt termination msg */
			setitfoind(1, &z, &z);
			}
		nos = netopt(cpx,&net_status,&net_nodes,&net_arcs,&net_its);
		if (nos)
			net_status = -5;
		nosp = &net_status;
		netopting = 0;
		if (!itfoind)
			setitfoind(0, &z, &z);
		if (nos != -2)
			status = Optimize(Env,cpx);
		}
	else
		status = Optimize(Env,cpx);
	Times[2] = xectim_();
	if (endbas && !nint)
		write_basis(cpx);
#ifdef CPLEX_MIP
	if (nint && endtree)
		treeio(cpx, endtree, "write", CPXtreewrite);
#endif
	amplout(cpx, &d, status, nint, nosp);
 done:
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
