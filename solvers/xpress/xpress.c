/* An interface from AMPL to Xpress-MP */
/* 1997-2008 Y. Colombani and others, Dash Optimization */
/* Requires Xpress-Optimizer libraries 18.00 or later */

/* Adjustments and additions to keywords and their descriptions, */
/* modifications for QPs, etc., by David M. Gay (2005). */

/* Remark: The options list MUST be in alphabetical order */

/* Some of the material in this file is derived from sample
   AMPL/solver interfaces that are available from netlib and
   which bear copyright notices of the following form... */

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

#include "xprs.h"
#include "nlp.h"
#include "getstub.h"
#include <signal.h>
#include <stdarg.h>

#define MININT (-XPRS_MAXINT - 1)

#ifndef XPRESS
#define XPRESS NULL
#endif

/* Basic status values returned by XPRESS function getbasis() */
#define XP_NBASLO 0  /* Vector non-basic at lower bound */
#define XP_BASIC  1  /* Vector basic */
#define XP_NBASUP 2  /* Vector non-basic at upper bound */


/* Define problem pointer */

XPRSprob prob;

typedef struct
 dims {
	double  *x;
	double  *y;
	int  *cstat;
	int  *rstat;
	SufDesc *csd;
	SufDesc *rsd;
	int miqp;
  } dims;

 static int mipststat = 1, Ray = 0, Round = 1, sos = 1, sos2 = 1;

static char *set_known(Option_Info *oi, keyword *kw, char *v);
static char *set_int(Option_Info *oi, keyword *kw, char *v);
static char *set_dbl(Option_Info *oi, keyword *kw, char *v);
static char *set_fln(Option_Info *oi, keyword *kw, char *v);
static void nonlin(int n, char *what);
static void amplin(char *stub, char *argv[], dims*);
static void amplout(dims*);
static void show_times(void);
static void xperror(const char *where, ...);
static void killtempprob(void);
static void mip_priorities(void);

static ASL *asl;
struct LU_bounds {real lower, upper;};
static double objadj;

enum glstat_e {
      GLSTAT_NOPROB           = 0,
      GLSTAT_LP_UNFINISHED    = 1,
      GLSTAT_LP_FINISHED      = 2,
      GLSTAT_UNFINISHED_NOSOL = 3,
      GLSTAT_UNFINISHED_SOL   = 4,
      GLSTAT_FINISHED_NOSOL   = 5,
      GLSTAT_FINISHED_SOL     = 6,
      GLSTAT_MIP_UNBOUNDED    = 7
     };

enum lpstat_e {
      LPSTAT_OPTIMAL      = 1,
      LPSTAT_INFEASIBLE   = 2,
      LPSTAT_CUTOFF       = 3,
      LPSTAT_UNFINISHED   = 4,
      LPSTAT_UNBOUNDED    = 5,
      LPSTAT_CUTOFFINDUAL = 6,
      LPSTAT_UNSOLVED	  = 7,
      LPSTAT_NONCONVEX    = 8
     };

enum known_parameters {
      set_primal,        /* Choose algorithm */
#ifdef RWA_DEBUG
      set_debug,
#endif
      set_dual,
      set_barrier,
      set_maxim,         /* What to do */
      set_minim,
      set_relax,
      set_timing,
      set_iis,
     };

 typedef struct Defer_setting Defer_setting;
 struct
Defer_setting {
	Defer_setting *next;
	int (*Setter)(Defer_setting*);
	keyword *kw;
	int ipar;
	union {
		int i;
		double d;
		} u;
	};

 static Defer_setting *Defer1, **Defer_next = &Defer1;

 static Defer_setting *
new_Defer_setting(int (*Dset)(Defer_setting*), keyword *kw, int ipar)
{
	Defer_setting *ds = (Defer_setting*)M1alloc(sizeof(Defer_setting));
	*Defer_next = ds;
	Defer_next = &ds->next;
	ds->Setter = Dset;
	ds->kw = kw;
	ds->ipar = -ipar;
	return ds;
	}

 void
Do_Defer(VOID)
{
	Defer_setting *ds;
	int nbad;

	*Defer_next = 0;
	nbad = 0;
	for(ds = Defer1; ds; ds = ds->next)
		nbad += (*ds->Setter)(ds);
	if (nbad)
		exit(2);
	}

char probname[L_tmpnam+4]; /* Use a temporary problem name */
char *endbasis=NULL, *startbasis=NULL, *logfile=NULL;
static double Times[4];      /* Timing stats */
int timing=0;
int iis_find=0;              /* find IIS */
#ifdef RWA_DEBUG
char debugopt[]="                   ";
#endif

static int nobj=1;           /* Which objective we optimise */
                             /* Which optimisation function */
static int  (XPRS_CC *Optimise)(XPRSprob ,const char *)=NULL;

static int prtmsg = 3;       /* Message output level */

char optimopt[3]={'X', '\0', '\0'};  /* Which algorithm */
/* [0] is b|d|X(meaning primal will be used - warning: non standard) */
/* [1] is l|g|blank */

static char
	autoperturb_desc[]	= "whether to introduce perturbations when the simplex\n\
			method encounters too many degnerate pivots:\n\
			1 = yes (default); 0 = no",
	backtrack_desc[]	= "choice of next node when solving MIP problems:\n\
			 1 = withdrawn; formerly choice 2 until a feasible\n\
				integer solution has been found, then\n\
				Forrest-Hirst-Tomlin choice\n\
			 2 = node with best estimated solution\n\
			 3 = node with best bound on the solution (default)\n\
			 4 = deepest node (depth-first search)\n\
			 5 = highest node (breadth-first search)\n\
			 6 = earliest-created node\n\
			 7 = most recently created node\n\
			 8 = random choice\n\
			 9 = node with fewest LP relaxation infeasibilities\n\
			10 = combination of 2 and 9\n\
			11 = combination of 2 and 4",
	backtracktie_desc[]	= "how to break ties for the next MIP node:\n\
			same choices as for \"backtrack\"",
	barcrash_desc[]		= "choice of crash procedure for crossover:\n\
			0 = no crash\n\
			1-6 = available strategies:\n\
			1 = most conservative, 6 = most agreessive",
	bardualstop_desc[]	= "barrier method convergence tolerance on\n\
			dual infeasibilities (default = 1e-8)",
	bargapstop_desc[]	= "barrier method convergence tolerance on\n\
			the relative duality gap (default = 0)",
	barindeflimit_desc[]	=
		"maximum indefinite factorizations to tolerate in the barrier\n\
			algorithm for solving a QP: stop when the limit is hit;\n\
			default = 15",
	bariterlimit_desc[]	= "maximum number of Newton Barrier iterations (default 200)",
	barorder_desc[]		= "Cholesky factorization pivot order for barrier algorithm:\n\
			0 = automatic choice (default)\n\
			1 = minimum degree\n\
			2 = minimum local fill\n\
			3 = nested disection",
	barout_desc[]		= "amount of output for the barrier method:\n\
			0 = no output\n\
			1 = each iteration (default)",
	barpresolve_desc[]	= "level of barrier-specific presolve effort:\n\
			0 = use standard presolve (default)",
	barprimalstop_desc[]	= "barrier method convergence tolerance on\n\
			primal infeasibilities (default = 1e-8)",
	barrier_desc[]		= "[no assignment] use the Newton Barrier algorithm",
	barstart_desc[]		= "choice of starting point for barrier method:\n\
			0 = automatic choice (default)\n\
			1 = heuristics based on magnitudes of matrix entries\n\
			2 = use pseudoinverse of constraint matrix",
	barstepstop_desc[]	= "barrier method convergence tolerance: stop when\n\
		step size <= barstepstop (default = 1e-10)",
	barthreads_desc[]	= "number of threads (default 1) used in the Newton Barrier\n\
		algorithm",
	bigmmethod_desc[]	= "0 = phase I/II, 1 = BigM method (default)",
	branchchoice_desc[]	= "whether to explore branch with min. or max.\n\
			estimate first:\n\
			0 = explore branch with min. estimate first (default)\n\
			1 = explore branch with max. estimate first\n\
			2 = if an incumbent solution exists, first explore\n\
				the branch satisfied by the incumbent;\n\
				otherwise use choice 0 (min. est. first)"
#if XPVERSION >= 22
		"\n\
			3 = explore the first branch that moves the branching\n\
				variable away from its value at the root node;\n\
				if the branching entity is not a simple\n\
				variable, assume branchchoice=0"
#endif
			,
	branchdisj_desc[]	= "whether to branch on general split disjunctions\n\
			while solving MIPs:\n\
			-1 = automatic choice (default)\n\
			 0 = disabled\n\
			 1 = cautious strategy: create branches only for\n\
				general integers with a wide range\n\
			 2 = moderate strategy\n\
			 3 = agressive strategy:  create disjunctive branches\n\
				for both binary and integer variables",
	branchstruct_desc[]	= "whether to search for special structure\n\
			during branch and bound:\n\
			-1 = automatic choice (default)\n\
			 0 = no\n\
			 1 = yes",
	breadthfirst_desc[]	=
		"number of MIP nodes included in best-first search\n\
			(default 10) before switching to local-first search",
	cachesize_desc[]	= "cache size in Kbytes -- relevant to Newton Barrier:\n\
			-1 = determined automatically for Intel\n\
			default = system-dependent (-1 for Intel)",
	choleskyalg_desc[]	= "type of Cholesky factorization used:\n\
			0 = Push (default), 1 = Pull",
	choleskytol_desc[]	= "zero tolerance for Cholesky pivots (default 1e-15)\n\
			in the Newton Barrier algorithm",
	convexitychk_desc[]	= "whether to check convexity before solving:\n\
			0 = no\n\
			1 = yes (default)",
#ifdef XPRS_CORESPERCPU
	corespercpu_desc[]	=
		"number of cores to assume per cpu; default = -1 ==> number\n\
			detected; the cachesize is divided by this number",
#endif
	covercuts_desc[]	=
		"for MIPS, the number of rounds of lifted-cover inequalities\n\
			at the top node (default = -1 = automatic choice)",
	cputime_desc[]		= "which times to report when logfile is speccified:\n\
			0 = elapsed time\n\
			1 = CPU time (default)",
	crash_desc[]		= "type of simplex crash:\n\
			0 = none\n\
			1 = one-pass search for singletons\n\
			2 = multi-pass search for singletons (default)\n\
			3 = multi-pass search including slacks\n\
			4 = at most 10 passes, only considering slacks\n\
			    at the end\n\
			n = (for n > 10) like 4, but at most n-10 passes",
	crossover_desc[]	= "whether to find a simplex basis after the barrier alg.:\n\
			1 = yes (default), 0 = no",
	cutdepth_desc[]		= "maximum MIP tree depth at which to generate cuts:\n\
			0  = no cuts\n\
			-1 = automatic choice (default)",
	cutfactor_desc[]	= "limit on number of cuts and cut coefficients\n\
			added while solving MIPs:\n\
			-1 = automatic choice (default)\n\
			 0 = do not add cuts\n\
			 > 0 ==> multiple of number of original constraints",
	cutfreq_desc[]		=
		"MIP cuts are only generated at tree depths that are integer\n\
			multiples of cutfreq; -1 = automatic choice (default)",
	cutselect_desc[]	= "detailed control of cuts at MIP root node:\n\
			sum of\n\
			   16 = clique cuts\n\
			   32 = mixed-integer founding (MIR) cuts\n\
			   64 = lifted cover cuts\n\
			 1024 = flow path cuts\n\
			 2048 = implication cuts\n\
			 4096 = automatic lift-and-project strategy\n\
			 8192 = disable cutting from cut rows\n\
			16384 = lifted GUB cover cuts\n\
			   -1 = all available cuts (default)",
	cutstrategy_desc[]	=
		"how aggressively to generate MIP cuts; more ==> fewer nodes\n\
			but more time per node:\n\
			-1 = automatic choice (default)\n\
			 0 = no cuts\n\
			 1 = conservative strategy\n\
			 2 = moderate strategy\n\
			 3 = aggressive strategy",
	defaultalg_desc[]	=
		"algorithm to use when none of \"barrier\", \"dual\", or \"primal\"\n\
			is specified:\n\
			1 = automatic choice (default)\n\
			2 = dual simplex\n\
			3 = primal simplex\n\
			4 = Newton Barrier",
#ifdef XPRS_DEGRADEFACTOR
	degradefactor_desc[]	=
		"factor to multiply estimated degradations in MIP objective\n\
			value from exploring an unexplored node; default = 1.0",
#endif
	densecollimit_desc[]	=
		"number of nonzeros above which a column is treated as dense\n\
			in the barrier algorithm's Cholesky factorization:\n\
			0 = automatic choice (default)",
	deterministic_desc[]	= "whether a MIP search should be deterministic:\n\
			0 = no\n\
			1 = yes (default)",
	dual_desc[]		= "[no assignment] use the dual simplex algorithm",
	dualgradient_desc[]	= "dual simplex pricing strategy:\n\
			-1 = automatic choice\n\
			 0 = Devex\n\
			 1 = steepest edge",
#ifdef XPRS_DUALIZE
	dualize_desc[]		= "whether the barrier algorithm should solve dual problems:\n\
			-1 = automatic choice (default)\n\
			 0 = solve primal problem\n\
			 1 = solve dual problem",
#endif /*XPRS_DUALIZE*/
	dualstrategy_desc[]	= "how to remove infeasibilities when re-optimizing\n\
			with the dual algorithm during MIP solves:\n\
			0 = use primal algorithm\n\
			1 = use dual algorithm (default)",
	eigenvaltol_desc[]	= "regard the matrix in a quadratic form as indefinite if its\n\
			smallest eigvenalue is < -eigevnaltol; default = 1e-6",
	elimtol_desc[]		= "Markowitz tolerance for the elimination phase of\n\
			XPRESS's presolve; default = 0.001",
	etatol_desc[] =		"zero tolerance on eta elements; default varies with XPRESS\n\
			version; default = 1e-12 or 1e-13 with some versions.\n\
			Use etatol=? to see the current value.",
	feaspump_desc[]		= "whether to run the Feasibility Pump heuristic at the top\n\
		node during branch-and-bound:  one of\n\
			0 = no (default)\n\
			1 = yes\n\
			2 = only if other heurstics found no integer solution",
	gomcuts_desc[]		= "gomory cuts at root: -1 = automatic choice (default)",
	heureffort_desc[]	= "factor (default 1.0) affecting how much work local search\n\
			heuristics should do.  Higher values cause more local\n\
			searches over larger neighborhoods",
	hdive_rand_desc[]	= "value between 0 and 1 inclusive affecting randomization\n\
			in the diving heuristic:  0 (default) ==> none;\n\
			1 ==> full;\n\
			intermediate values ==> intermediate behavior",
	hdive_speed_desc[]	= "controls tradeoff between speed and solution quality\n\
			in the diving heuristic:  an integer between -2 and 3:\n\
			-2 = automatic bias toward quality\n\
			-1 = automatic bias toward speed\n\
			 0 = emphasize quality\n\
			 4 = emphasize speed\n\
			 1-3 = intermediate emphasis",
	hdive_strategy_desc[]	= "strategy for diving heuristic:  integer between -1 and 10:\n\
			-1 = automatic choice (default)\n\
			 0 = do not use the diving heursistic\n\
			1-10 = preset strategies for diving",
	heurdepth_desc[]	=
#if XPVERSION >= 22
		"deprecated:  no longer has any effect:\n\t\t"
#endif
		"maximum depth of branch-and-bound tree search at which to apply\n\
			heuristics; 0 = no heuristics (default)",
	heurfreq_desc[]		=
		"during branch and bound, heuristics are applied at nodes\n\
			whose depth from the root is zero modulo heurfreq\n\
			(default 5)",
	heurmaxsol_desc[]	=
#if XPVERSION >= 22
		"deprecated:  no longer has any effect:\n\t\t"
#endif
		"maximum number of heuristic solutions to find during\n\
			branch-and-bound tree search (default 10)",
	heurnodes_desc[]	=
#if XPVERSION >= 22
		"deprecated:  no longer has any effect:\n\t\t"
#endif
		"maximum nodes at which to use heuristics during\n\
			branch-and-bound tree search (default 1000)",
	heurroot_desc[]		= "bit vector controlling local search heuristics to\n\
			apply at the root node:  sum of\n\
			  1 = large-neighborhood search: may be slow, but may\n\
				find solutions far from the incumbent\n\
			  2 = small-neighborhood search about node LP solution\n\
			  4 = small-neighborhood search about integer solutions\n\
			default = 2",
	heursearch_desc[]	= "how often the local search heurstic should be run\n\
			during branch-and-bound:\n\
			-1 = automatic choice (default)\n\
			 0 = never\n\
			 n > 0 ==> every n nodes",
	heurstrategy_desc[]	= "heuristic strategy for branch and bound: one of\n\
			-1 = automatic choice (default)\n\
			 0 = no heuristics\n\
			 1 = rounding heuristics (sometimes useful)",
	heurthreads_desc[]	= "number of threads for the root node\n\
			of branch-and-bound:\n\
			-1 = determined from \"threads\" keyword\n\
			 0 = no separate threads (default)\n\
			 n > 0 ==> use n threaeds",
	heurtree_desc[]		= "heuristics to apply during tree search:  sum of\n\
			the same values as for heurroot; default 2",
	iis_desc[]	=
		"[no assignment] if the problem is infeasible, find an\n\
			Irreducible Independent Set of infeasible constraints",
	indlinbigm_desc[]	= "largest \"big M\" value to use in converting indicator\n\
			constraints to regular constraints;\n\
			default = 1e5",
	indprelinbigm_desc[]	= "largest \"big M\" value to use in converting indicator\n\
			constraints to regular constraints during\n\
			XPRESS presolve; default = 100.0",
	invertfreq_desc[]	=
		"maximum simplex iterations before refactoring the basis:\n\
			-1 = automatic choice (default)",
	invertmin_desc[]	=
		"minimum simplex iterations before refactoring the basis:\n\
			default = 3",
	keepbasis_desc[]	= "basis choice for the next LP iteration:\n\
			0 = ignore previous basis\n\
			1 = use previous basis (default)\n\
			2 = use previous basis only if the number of basic\n\
				variables == number of constraints",
	keepnrows_desc[]	=
		"1 (default) if unconstrained rows are to be kept, else 0",
	lnpbest_desc[]		=
		"number of global infeasible entities for which to create\n\
			lift-and-project cuts during each round of Gomory cuts\n\
			at the top node (default 50)",
	lnpiterlimit_desc[]	=
		"maximum iterations for each lift-and-project cut\n\
			(default 10)",
	localchoice_desc[]	= "when to backtrack between two child nodes\n\
		during a \"dive\":\n\
			1 = (default) never backtrack from the first child\n\
				unless it is dropped (i.e., is infeasible\n\
				or cut off)\n\
			2 = always solve both nodes first\n\
			3 = automatic choice",
	logfile_desc[]		= "name of log file (default = no log file)",
	lpthreads_desc[]	= "number of threads in concurrent LP solves:\n\
			-1 = determine from \"threads\" keyword (default)\n\
			n > 0 ==> use n threads",
	lpiterlimit_desc[]	=
		"simplex iteration limit (default 2147483645)",
	lplog_desc[]		=
		"frequency of printing simplex iteration log (default 100)",
	markowitztol_desc[]	=
		"Markowitz tolerance used when factoring the basis matrix\n\
			(default 0.01)",
	matrixtol_desc[]	= "zero tolerance on matrix elements (default 1e-9)",
	maxcuttime_desc[]	= "maximum time (CPU seconds) to spend generating cuts\n\
			and reoptimizing (default = 0 ==> no limit)",
	maxiis_desc[]		= "maximum number of Irreducible Infeasible Sets to find:\n\
			-1 = no limit (default)\n\
			 0 = none",
	maxim_desc[]		= "[no assignment] force maximization of the objective",
	maxlocalbt_desc[]	= "max height above current node to look\n\
			for a local backtrack candidate node; default = 1",
	maxmipsol_desc[]	= "maximum number of integer solutions to find:\n\
			0 = no limit (default)",
	maxnode_desc[]		=
		"maximum number of MIP nodes to explore; default = 100000000",
	maxpagelines_desc[]	= "maximum output lines between page breaks in logfile\n\
			(default 23)",
	maxlogcale_desc[]	= "max log2 of factors used in scaling; must be\n\
			>= 0 and <= 64; default 64",
#ifdef XPRS_MAXSLAVE
	maxslaves_desc[]	= "how many processors to use in parallel when solving\n\
			mixed-integer problems:  must be at most the number of\n\
			processors available and licensed",
#endif
	maxtime_desc[]		= "maximum solution time allowed (default = 0 ==> no limit)",
	minim_desc[]		= "[no assignment] force minimization of the objective",
	mipabscutoff_desc[]	=
		"initial MIP cutoff:  ignore MIP nodes with objective values\n\
			worse than mipabscutoff; default = 1e40 for\n\
			minimization, -1e40 for maximization",
	mipabsstop_desc[]	=
		"stop MIP search if abs(MIPOBJVAL - BESTBOUND) <= mipabsstop\n\
			(default 0)",
	mipaddcutoff_desc[]	=
		"amount to add to the objective function of the best integer\n\
			solution found to give the new MIP cutoff\n\
			(default -1e-5)",
	miplog_desc[]		= "MIP printing level to logfile (default -100):\n\
			-n = print summary line every n MIP nodes\n\
			 0 = no MIP summary lines\n\
			 1 = only print a summary at the end\n\
			 2 = log each solution found\n\
			 3 = log each node",
	mipops_desc[]		= "MIP solver options:  one of\n\
			0 = traditional primal first phase (default)\n\
			1 = Big M primal first phase\n\
			2 = traditional dual first\n\
			3 = Big M dual first\n\
			4 = always use artificial bounds in dual\n\
			5 = use original basis only when warmstarting\n\
			6 = skip primal bound flips for ranged primals\n\
			7 = also do single-pivot crash\n\
			8 = suppress agressive dual perturbations",
	mippresolve_desc[]	= "MIP presolve done at each node: sum of\n\
			1 = reduced-cost fixing\n\
			2 = logical preprocessing of binary variables\n\
			4 = probing of binary variables\n\
		default determined from constraint-matrix properties",
	miprelcutoff_desc[]	=
		"fraction of best integer solution found to add to MIP cutoff\n\
			(default 1e-4)",
	miprelstop_desc[]	= "stop MIP search if\n\
		  abs(MIPOBJVAL - BESTBOUND) < miprelstop * abs(BESTBOUND)\n\
			(default = 0)",
	mipstartstatus_desc[]	= "use incoming statuses on MIP problems (default 1 = yes)",
#ifdef XPRS_MIPTARGET
	miptarget_desc[]	=
		"initial MIP target objective value (default 1e40),\n\
			used in some node-selection rules and updated as MIP\n\
			solutions are found",
#endif
	mipthreads_desc[]	= "number of threads to use solving mixed-integer\n\
		programming problems:\n\
			-1 = use \"threads\" keyword (default)\n\
			n > 0 ==> use n threads",
	miptol_desc[]		= "integer feasibility tolerance (default 5e-6)",
	nodefilebias_desc[]	= "a value between 0 and 1 (inclusive) that influences\n\
			operations when \"treememlimit\" (on how much of the\n\
			branch-and-bound tree should be kept in memory) has\n\
			been exceeded:\n\
			  0 ==> compress every node before writing anything to\n\
				the \"nodefile\";\n\
			  1 ==> write nodes to the \"nodefile\" immediately;\n\
			values between 0 and 1 give intermediate behavior.\n\
			Default = 0.5",
	nodeselection_desc[]	= "next MIP node control (default determined from\n\
			matrix characteristics):\n\
			1 = local first:  choose among descendant and sibling\n\
			    nodes if available, else from all outstanding nodes\n\
			2 = best first of all outstanding nodes\n\
			3 = local depth first:  choose among descendant and\n\
			    sibling nodes if available, else from deepest nodes\n\
			4 = best first for breadthfirst nodes, then local first\n\
			5 = pure depth first:  choose among deepest nodes",
	optimalitytol_desc[]	= "tolerance on reduced cost (default 1e-6)",
	outlev_desc[]		= "message level:\n\
			1 = all\n\
			2 = information\n\
			3 = warnings & errors only (default)\n\
			4 = errors\n\
			5 = none",
	outputtol_desc[]	= "zero tolerance on print values (default 1e-5)",
	penalty_desc[]		= "minimum absolute penalty variable coefficient;\n\
	default = automatic choice",
	perturb_desc[]		= "perturb factor if autoperturb is set to 1;\n\
		0 = default = automatic choice",
	pivottol_desc[]		= "zero tolerance for pivots; default = 1e-9",
	ppfactor_desc[]		= "partial-pricing candidate-list size factor (default 1.0)",
	precoefelim_desc[]	= "whether XPRESSMP's presolve should recombine constraints:\n\
			0 = no,\n\
			1 = yes, as many as possible\n\
			2 = yes, cautiously (default)",
	predomcol_desc[]	= "whether XPRESSMP's presolve should remove variables\n\
			when solving MIP problems:\n\
			-1 = automatic choice (default)\n\
			 0 = no\n\
			 1 = yes, cautiously\n\
			 2 = yes, check all candidates",
	predomrow_desc[]	= "whether XPRESSMP's presolve should remove constraints\n\
		when solving MIP problems:\n\
			-1 = automatic choice (default)\n\
			 0 = no\n\
			 1 = yes, cautiously\n\
			 2 = yes, medium strategy\n\
			 3 = yes, check all candidates",
	preprobing_desc[]	= "how much probing on binary variables to do during\n\
			XPRESSMP's presolve:\n\
			-1 = automatic choice (default)\n\
			 0 = none\n\
			 1 = light probing\n\
			 2 = full probing\n\
			 3 = repeated full probing",
	presolve_desc[]		= "whether to use XPRESS's presolver:\n\
			0 = no\n\
			1 = yes, removing redundant bounds (default)\n\
			2 = yes, retaining redundant bounds",
	presolveops_desc[]	= "reductions to use in XPRESSMP's presolve:\n\
			sum of\n\
			    1 = 2^0  = remove singleton columns\n\
			    2 = 2^1  = remove singleton constraints (rows)\n\
			    4 = 2^2  = forcing row removal (whatever that is)\n\
			    8 = 2^3  = dual reductions\n\
			   16 = 2^4  = redundant constraint (row) removal\n\
			   32 = 2^5  = duplicate variable removal\n\
			   64 = 2^6  = duplicate constraint removal\n\
			  128 = 2^7  = strong dual reductions\n\
			  256 = 2^8  = variable eliminations\n\
			  512 = 2^9  = no IP reductions\n\
			 1024 = 2^10 = no semicontinuous variable detection\n\
			 2048 = 2^11 = no advanced IP reductions\n\
			16384 = 2^14 = remove linearly dependent constraints\n\
			32768 = 2^15 = no integer variable and SOS detection\n\
		default = 511 (bits 0-8 set)",
	pricingalg_desc[]	= "primal simplex pricing method:\n\
			-1 = partial pricing\n\
			 0 = automatic choice (default)\n\
			 1 = Devex pricing",
	primal_desc[]		= "[no assignment] use the primal simplex algorithm",
#ifdef XPRS_PRIMALOPS
	primalops_desc[]	= "primal simplex options:  sum of\n\
			1 = 2^0 = agressive dj scaling\n\
			2 = 2^1 = conventional dj scaling\n\
			4 = 2^2 = reluctant switching back to partial pricing\n\
			8 = 2^3 = dynamic switching between cheap and expensive pricing\n\
			default = all of the above; if bits 0 and 1 are the same (both on or\n\
			both off), choose dj scaling automatically",
#endif
	primalunshift_desc[]	= "whether the primal alg. calls the dual to unshift:\n\
			0 = yes (default)\n\
			1 = no",
	pseudocost_desc[]	=
		"default pseudo-cost assumed for forcing an integer variable\n\
			to an integer value (default = 0.01)",
	pseudocost_ud_desc[] = 	"how to update pseudocosts during branch-and-bound:\n\
			-1 = automatic choice (default)\n\
			 0 = no updates\n\
			 1 = use only regular branches\n\
			 2 = use regular and strong branch results\n\
			 3 = use results from all nodes",
	quadunshift_desc[]	= "whether quadratic simplex should do an extra\n\
			purification after finding a solution:\n\
			-1 = automatic choice (default)\n\
			 0 = no\n\
			 1 = yes",
	ray_desc[]		=
		"whether to return a ray of unboundedness in suffix .unbdd:\n\
			0 ==> no (default)\n\
			1 ==> yes, after suppressing XPRESS's presolve\n\
			2 ==> yes, without suppressing XPRESS's presolve\n\
			The last setting (ray=2) may give wrong results when\n\
			XPRESS's presolve detects infeasibility.  Both ray=1\n\
			and ray=2 cause reoptimization with primal simplex if\n\
			some other algorithm was used.  No ray is returned for\n\
			MIP problems.",
	relaxtreemem_desc[]	= "fraction of memory limit by which to relax \"treememlimit\"\n\
			when too much structural data appears; default 0.1",
	relpivottol_desc[]	= "relative pivot tolerance (default 1e-6)",
	repairindefq_desc[]	= "whether to repair indefinite quadratic forms:\n\
			0 = yes\n\
			1 = no (default)",
	rootpresolve_desc[]	= "whether to presolve after root cutting and heuristics:\n\
			-1 = automatic choice (default)\n\
			 0 = no\n\
			 1 = yes",
	round_desc[]		=
		"whether to round integer variables to integral values before\n\
			returning the solution, and whether to report that\n\
			XPRESS returned noninteger values for integer values\n\
			(default 1):  sum of\n\
			 1 ==> round nonintegral integer variables\n\
			 2 ==> do not modify solve_result\n\
			 4 ==> do not modify solve_message\n\
			 8 ==> modify even if maxerr < 1e-9\n\
			Modifications take place only if XPRESS assigned\n\
			nonintegral values to one or more integer variables,\n\
			and (for round < 8) only if the maximum deviation from\n\
			integrality exceeded 1e-9.",
	sbbest_desc[]		= "For MIP problems, the number of infeasible\n\
		global entities on which to perform strong branching\n\
			(default -1 = automatic)",
	sbeffort_desc[]		= "multiplier on strong-branching controls that\n\
		are set to \"automatic\"; default = 1.0",
	sbestimate_desc[]	= "how to compute pseudo costs from the local node\n\
		when selecting an infeasible entity to branch on:\n\
			-1 = automatic choice (default)\n\
			1-6 = particular strategies (not described)",
	sbiterlimit_desc[]	=
		"Number of dual iterations to perform the strong branching\n\
			(default 0 for none)",
	sbselect_desc[]		= "size of candidate list for strong branching:\n\
			-2 = low-effort automatic choice (default)\n\
			-1 = high-effort automatic choice\n\
			n >= 0 ==> include max(n, sbbest) candidates",
	scaling_desc[]		=
		"how to scale the constraint matrix before optimizing: sum of\n\
			   1 = 2^0 = row scaling\n\
			   2 = 2^1 = column scaling\n\
			   4 = 2^2 = row scaling again\n\
			   8 = 2^3 = maximum scaling\n\
			  16 = 2^4 = Curtis-Reid\n\
			  32 = 2^5 = scale by maximum element (rather\n\
					than by geometric mean)\n\
			 128 = 2^7 = objective-function scaling\n\
			 256 = 2^8 = excluding quadratic part of constraint\n\
					when calculating scaling factors\n\
			 512 = 2^9 = scale before presolve\n\
			1024 = 2^10 = do not scale constraints (rows) up\n\
			2048 = 2^11 = do not scale variables up\n\
			default = 163",
#ifdef XPRS_SLEEPONTHREADWAIT
	sleeponthreadwait_desc[]	=
		"whether threads should sleep while awaiting work:\n\
			0 = no (busy-wait; default)\n\
			1 = yes (sleep; may add overhead)",
#endif
	sos_desc[]		= "whether to use explicit SOS information (default 1 = yes)",
	sos2_desc[]		= "whether to use implicit SOS information (default 1 = yes)",
	sosreftol_desc[]	= "minimum relative gap between reference row entries\n\
			default = 1e-6",
	tempbounds_desc[]	= "whether dual simplex should put temporary bounds on\n\
		unbounded variables:\n\
			-1 = automatic choice (default)\n\
			 0 = no\n\
			 1 = yes",
	threads_desc[]		= "default number of threads to use:\n\
			-1 = automatic choice (based on hardware)\n\
			 n > 0 ==> use n threads",
	trace_desc[]		= "whether to explain infeasibility:\n\
			0 = no (default)\n\
			1 = yes",
	treecovercuts_desc[]	=
		"number of rounds of lifted-cover inequalities at MIP nodes\n\
			other than the top node (cf covercuts); default = 1",
	treecompress_desc[]	= "level of effort at data compression when branch-and-bound\n\
			memory exceedss \"treememlimit\":  higher ==> greater\n\
			effort (taking more time); default = 2",
	treecuts_desc[] 	= "cuts to generate at nodes during tree search:\n\
		sum of\n\
			    32 = 2^5  = clique cuts\n\
			    64 = 2^6  = mixed-integer rounding (MIR) cuts\n\
			    64 = 2^7  = lifted-cover cuts\n\
			  2048 = 2^11 = flow-path cuts\n\
			  4096 = 2^12 = implication cuts\n\
			  8192 = 2^13 = lift-and-project cuts\n\
			 16384 = 2^14 = disable cutting from row cuts\n\
			 32768 = 2^15 = lifted GUB cover cuts\n\
			 65536 = 2^16 = zero-half cuts\n\
			131072 = 2^17 = indicator cuts\n\
		default = 521983",
	treegomcuts_desc[]	=
		"number of rounds of Gomory cuts to generate at MIP nodes\n\
			other than the top node (cf covercuts); default = 1",
	treeoutlev_desc[]	= "how much to report about branch-and-bound trees\n\
		(if allowed by outlev):  sum of\n\
			1 = regular summaries\n\
			2 = report tree compression and output to nodefile\n\
			default = 3",
	treememlimit_desc[]	= "an integer: soft limit in megabytes on memory to use for\n\
		branch-and-bound trees.  Default = 0 ==> automatic choice.",
	treememtarget_desc[] 	= "fraction of \"treememlimit\" to try to recover by compression\n\
			or writing to nodefile when  \"treememlimit\" is\n\
			exceeded.  Default = 0.1",
	varselection_desc[]	=
		"how to score the integer variables at a MIP node, for\n\
			branching on a variable with minimum score:\n\
			-1 = automatic choice (default)\n\
			 1 = minimum of the 'up' and 'down' pseudo-costs\n\
			 2 = 'up' pseudo-cost + 'down' pseudo-cost\n\
			 3 = maximum of the 'up' and 'down' pseudo-costs plus\n\
			     twice their minimum\n\
			 4 = maximum of the 'up' and 'down' pseudo-costs\n\
			 5 = the 'down' pseudo-cost\n\
			 6 = the 'up' pseudo-cost";

        /* The list of Xpress-MP options */
        /******MUST BE IN ALPHABETIC ORDER!****/

static keyword keywds[]={
  KW("autoperturb",	set_int, XPRS_AUTOPERTURB,	autoperturb_desc),
  KW("backtrack",	set_int, XPRS_BACKTRACK,	backtrack_desc),
  KW("backtracktie",	set_int, XPRS_BACKTRACKTIE,	backtracktie_desc),
  KW("barcrash",	set_int, XPRS_BARCRASH,		barcrash_desc),
  KW("bardualstop",	set_dbl, XPRS_BARDUALSTOP,	bardualstop_desc),
  KW("bargapstop",	set_dbl, XPRS_BARGAPSTOP,	bargapstop_desc),
  KW("barindeflimit",	set_int, XPRS_BARINDEFLIMIT,	barindeflimit_desc),
  KW("bariterlimit",	set_int, XPRS_BARITERLIMIT,	bariterlimit_desc),
  KW("barorder",	set_int, XPRS_BARORDER,		barorder_desc),
  KW("baroutput",	set_int, XPRS_BAROUTPUT,	barout_desc),
  KW("barpresolve",	set_int, XPRS_BARPRESOLVEOPS,	barpresolve_desc),
  KW("barprimalstop",	set_dbl, XPRS_BARPRIMALSTOP,	barprimalstop_desc),
  KW("barrier",		set_known, set_barrier,		barrier_desc),
  KW("barstart",	set_int, XPRS_BARSTART,		barstart_desc),
  KW("barstepstop",	set_dbl, XPRS_BARSTEPSTOP,	barstepstop_desc),
  KW("barthreads",	set_int, XPRS_BARTHREADS,	barthreads_desc),
  KW("basisin",		set_fln, &startbasis,		"load initial basis from specified file"),
  KW("basisout",	set_fln,  &endbasis,		"save final basis to specified file"),
  KW("bigm",		set_dbl, XPRS_BIGM,		"infeasibility penalty (default 1024)"),
  KW("bigmmethod",	set_int, XPRS_BIGMMETHOD,	bigmmethod_desc),
  KW("branchchoice",	set_int, XPRS_BRANCHCHOICE,	branchchoice_desc),
  KW("branchdisj",	set_int, XPRS_BRANCHDISJ,	branchdisj_desc),
  KW("branchstruct",	set_int, XPRS_BRANCHSTRUCTURAL, branchstruct_desc),
  KW("breadthfirst",	set_int, XPRS_BREADTHFIRST,	breadthfirst_desc),
  KW("cachesize",	set_int, XPRS_CACHESIZE,	cachesize_desc),
  KW("choleskyalg",	set_int, XPRS_CHOLESKYALG,	choleskyalg_desc),
  KW("choleskytol",	set_dbl, XPRS_CHOLESKYTOL,	choleskytol_desc),
  KW("convexitychk",	set_int, XPRS_IFCHECKCONVEXITY,	convexitychk_desc),
#ifdef XPRS_CORESPERCPU
  KW("corespercpu",	set_int, XPRS_CORESPERCPU,	corespercpu_desc),
#endif
  KW("covercuts",	set_int, XPRS_COVERCUTS,	covercuts_desc),
  KW("cputime",		set_int, XPRS_CPUTIME,		cputime_desc),
  KW("crash",		set_int, XPRS_CRASH,		crash_desc),
  KW("crossover",	set_int, XPRS_CROSSOVER,	crossover_desc),
  KW("cutdepth",	set_int, XPRS_CUTDEPTH,		cutdepth_desc),
  KW("cutfactor",	set_dbl, XPRS_CUTFACTOR,	cutfactor_desc),
  KW("cutfreq",		set_int, XPRS_CUTFREQ,		cutfreq_desc),
  KW("cutselect",	set_int, XPRS_CUTSELECT,	cutselect_desc),
  KW("cutstrategy",	set_int, XPRS_CUTSTRATEGY,	cutstrategy_desc),
#ifdef RWA_DEBUG
  KW("debug",		set_known, set_debug,		"RWA's debug switch [no assignment]"),
#endif
  KW("defaultalg",	set_int, XPRS_DEFAULTALG,	defaultalg_desc),
#ifdef XPRS_DEGRADEFACTOR
  KW("degradefactor",	set_dbl, XPRS_DEGRADEFACTOR,	degradefactor_desc),
#endif
  KW("densecollimit",	set_int, XPRS_DENSECOLLIMIT,	densecollimit_desc),
  KW("deterministic",	set_int, XPRS_DETERMINISTIC,	deterministic_desc),
  KW("dual",		set_known, set_dual,		dual_desc),
  KW("dualgradient",	set_int, XPRS_DUALGRADIENT,	dualgradient_desc),
#ifdef XPRS_DUALIZE
  KW("dualize",		set_int, XPRS_DUALIZE,		dualize_desc),
#endif
  KW("dualstrategy",	set_int, XPRS_DUALSTRATEGY,	dualstrategy_desc),
  KW("eigenvaltol",	set_dbl, XPRS_EIGENVALUETOL,	eigenvaltol_desc),
  KW("elimtol",		set_dbl, XPRS_ELIMTOL,		elimtol_desc),
  KW("etatol",		set_dbl, XPRS_ETATOL,		etatol_desc),
  KW("feaspump",	set_int, XPRS_FEASIBILITYPUMP,	feaspump_desc),
  KW("feastol",		set_dbl, XPRS_FEASTOL,		"zero tolerance on RHS; default = 1e-6"),
  KW("gomcuts",		set_int, XPRS_GOMCUTS,		gomcuts_desc),
  KW("hdive_rand",	set_dbl, XPRS_HEURDIVERANDOMIZE, hdive_rand_desc),
  KW("hdive_speed",	set_int, XPRS_HEURDIVESPEEDUP,	hdive_speed_desc),
  KW("hdive_strategy",	set_int, XPRS_HEURDIVESTRATEGY,	hdive_strategy_desc),
  KW("heurdepth",	set_int, XPRS_HEURDEPTH,	heurdepth_desc),
  KW("heureffort",	set_dbl, XPRS_HEURSEARCHEFFORT,	heureffort_desc),
  KW("heurfreq",	set_int, XPRS_HEURFREQ,		heurfreq_desc),
  KW("heurmaxsol",	set_int, XPRS_HEURMAXSOL,	heurmaxsol_desc),
  KW("heurnodes",	set_int, XPRS_HEURNODES,	heurnodes_desc),
  KW("heurroot",	set_int, XPRS_HEURSEARCHROOTSELECT, heurroot_desc),
  KW("heursearch",	set_int, XPRS_HEURSEARCHFREQ,	heursearch_desc),
  KW("heurstrategy",	set_int, XPRS_HEURSTRATEGY,	heurstrategy_desc),
  KW("heurthreads",	set_int, XPRS_HEURTHREADS,	heurthreads_desc),
  KW("heurtree",	set_int, XPRS_HEURSEARCHTREESELECT, heurtree_desc),
  KW("iis",		set_known, set_iis,		iis_desc),
  KW("indlinbigm",	set_dbl, XPRS_INDLINBIGM,	indlinbigm_desc),
  KW("indprelinbigm",	set_dbl, XPRS_INDPRELINBIGM,	indprelinbigm_desc),
  KW("invertfreq",	set_int, XPRS_INVERTFREQ,	invertfreq_desc),
  KW("invertmin",	set_int, XPRS_INVERTMIN,	invertmin_desc),
  KW("keepbasis",	set_int, XPRS_KEEPBASIS,	keepbasis_desc),
  KW("keepnrows", 	set_int, XPRS_KEEPNROWS,	keepnrows_desc),
  KW("lnpbest",		set_int, XPRS_LNPBEST,		lnpbest_desc),
  KW("lnpiterlimit",	set_int, XPRS_LNPITERLIMIT,	lnpiterlimit_desc),
  KW("localchoice",	set_int, XPRS_LOCALCHOICE,	localchoice_desc),
  KW("logfile",		set_fln, &logfile,		logfile_desc),
  KW("lpiterlimit",	set_int, XPRS_LPITERLIMIT,	lpiterlimit_desc),
  KW("lplog",		set_int, XPRS_LPLOG,		lplog_desc),
  KW("lpthreads",	set_int, XPRS_LPTHREADS,	lpthreads_desc),
  KW("markowitztol",	set_dbl, XPRS_MARKOWITZTOL, 	markowitztol_desc),
  KW("matrixtol",	set_dbl, XPRS_MATRIXTOL,	matrixtol_desc),
  KW("maxcuttime",	set_int, XPRS_MAXCUTTIME,	maxcuttime_desc),
  KW("maxiis",		set_int, XPRS_MAXIIS,		maxiis_desc),
  KW("maxim",		set_known, set_maxim,		maxim_desc),
  KW("maximise",	set_known, set_maxim,		maxim_desc),
  KW("maximize",	set_known, set_maxim,		maxim_desc),
  KW("maxlocalbt",	set_int, XPRS_MAXLOCALBACKTRACK, maxlocalbt_desc),
  KW("maxlogcale",	set_int, XPRS_MAXSCALEFACTOR,	maxlogcale_desc),
  KW("maxmipsol",	set_int, XPRS_MAXMIPSOL,	maxmipsol_desc),
  KW("maxnode",		set_int, XPRS_MAXNODE,		maxnode_desc),
  KW("maxpagelines",	set_int, XPRS_MAXPAGELINES,	maxpagelines_desc),
#ifdef XPRS_MAXSLAVE
  KW("maxslaves",	set_int, XPRS_MAXSLAVE,		maxslaves_desc),
#endif
  KW("maxtime",		set_int, XPRS_MAXTIME,		maxtime_desc),
  KW("minim",		set_known, set_minim,		minim_desc),
  KW("minimise",	set_known, set_minim,		minim_desc),
  KW("minimize",	set_known, set_minim,		minim_desc),
  KW("mipabscutoff",	set_dbl, XPRS_MIPABSCUTOFF,	mipabscutoff_desc),
  KW("mipabsstop",	set_dbl, XPRS_MIPABSSTOP,	mipabsstop_desc),
  KW("mipaddcutoff",	set_dbl, XPRS_MIPADDCUTOFF,	mipaddcutoff_desc),
  KW("miplog",		set_int, XPRS_MIPLOG,		miplog_desc),
  KW("mipops",		set_int, XPRS_QSIMPLEXOPS,	mipops_desc),
  KW("mippresolve",	set_int, XPRS_MIPPRESOLVE,	mippresolve_desc),
  KW("miprelcutoff",	set_dbl, XPRS_MIPRELCUTOFF,	miprelcutoff_desc),
  KW("miprelstop",	set_dbl, XPRS_MIPRELSTOP,	miprelstop_desc),
  KW("mipstartstatus",	I_val, &mipststat,		mipstartstatus_desc),
#ifdef XPRS_MIPTARGET
  KW("miptarget",	set_dbl, XPRS_MIPTARGET,	miptarget_desc),
#endif
  KW("mipthreads",	set_int, XPRS_MIPTHREADS,	mipthreads_desc),
  KW("miptol",		set_dbl, XPRS_MIPTOL,		miptol_desc),
  KW("nodefilebias",	set_dbl, XPRS_GLOBALFILEBIAS,	nodefilebias_desc),
  KW("nodeselection",	set_int, XPRS_NODESELECTION,	nodeselection_desc),
  KW("objno",		I_val, &nobj,			"objective number (0=none, 1=first...)"),
  KW("optimalitytol",	set_dbl, XPRS_OPTIMALITYTOL,	optimalitytol_desc),
  KW("outlev",		I_val,	 &prtmsg,		outlev_desc),
  KW("outputtol",	set_dbl, XPRS_OUTPUTTOL,	outputtol_desc),
  KW("penalty",		set_dbl, XPRS_PENALTY,		penalty_desc),
  KW("perturb",		set_dbl, XPRS_PERTURB,		perturb_desc),
  KW("pivottol",	set_dbl, XPRS_PIVOTTOL,		pivottol_desc),
  KW("ppfactor",	set_dbl, XPRS_PPFACTOR,		ppfactor_desc),
  KW("precoefelim",	set_int, XPRS_PRECOEFELIM,	precoefelim_desc),
  KW("predomcol",	set_int, XPRS_PREDOMCOL,	predomcol_desc),
  KW("predomrow",	set_int, XPRS_PREDOMROW,	predomrow_desc),
  KW("preprobing",	set_int, XPRS_PREPROBING,	preprobing_desc),
  KW("presolve",	set_int, XPRS_PRESOLVE,		presolve_desc),
  KW("presolveops",	set_int, XPRS_PRESOLVEOPS,	presolveops_desc),
  KW("pricingalg",	set_int, XPRS_PRICINGALG,	pricingalg_desc),
  KW("primal",		set_known, set_primal,		primal_desc),
#ifdef XPRS_PRIMALOPS
  KW("primalops",	set_int, XPRS_PRIMALOPS,	primalops_desc),
#endif
  KW("primalunshift",	set_int, XPRS_PRIMALUNSHIFT,	primalunshift_desc),
  KW("pseudocost",	set_dbl, XPRS_PSEUDOCOST,	pseudocost_desc),
  KW("pseudocost_ud",	set_int, XPRS_HISTORYCOSTS,	pseudocost_ud_desc),
  KW("quadunshift",	set_int, XPRS_QUADRATICUNSHIFT,	quadunshift_desc),
  KW("ray",		I_val, &Ray,			ray_desc),
  KW("relax",		set_known, set_relax,		"[no assignment] ignore integrality"),
  KW("relaxtreemem",	set_dbl, XPRS_RELAXTREEMEMORYLIMIT, relaxtreemem_desc),
  KW("relpivottol",	set_dbl, XPRS_RELPIVOTTOL,	relpivottol_desc),
  KW("repairindefq",	set_int, XPRS_REPAIRINDEFINITEQ, repairindefq_desc),
  KW("rootpresolve",	set_int, XPRS_ROOTPRESOLVE,	rootpresolve_desc),
  KW("round",		I_val, &Round,			round_desc),
  KW("sbbest",		set_int, XPRS_SBBEST,		sbbest_desc),
  KW("sbeffort",	set_dbl, XPRS_SBEFFORT,		sbeffort_desc),
  KW("sbestimate",	set_int, XPRS_SBESTIMATE,	sbestimate_desc),
  KW("sbiterlimit",	set_int, XPRS_SBITERLIMIT,	sbiterlimit_desc),
  KW("sbselect",	set_int, XPRS_SBSELECT,		sbselect_desc),
  KW("scaling",		set_int, XPRS_SCALING,		scaling_desc),
#ifdef XPRS_SLEEPONTHREADWAIT
  KW("sleeponthreadwait", set_int, XPRS_SLEEPONTHREADWAIT, sleeponthreadwait_desc),
#endif
  KW("sos",		I_val, &sos,			sos_desc),
  KW("sos2",		I_val, &sos2,			sos2_desc),
  KW("sosreftol",	set_dbl, XPRS_SOSREFTOL,	sosreftol_desc),
  KW("tempbounds",	set_int, XPRS_TEMPBOUNDS,	tempbounds_desc),
  KW("threads",		set_int, XPRS_THREADS,		threads_desc),
  KW("timing",		set_known, set_timing,		"[no assignment] give timing statistics"),
  KW("trace",		set_int, XPRS_TRACE,		trace_desc),
  KW("treecompress",	set_int, XPRS_TREECOMPRESSION,	treecompress_desc),
  KW("treecovercuts",	set_int, XPRS_TREECOVERCUTS,	treecovercuts_desc),
  KW("treecuts",	set_int, XPRS_TREECUTSELECT,	treecuts_desc),
  KW("treegomcuts",	set_int, XPRS_TREEGOMCUTS,	treegomcuts_desc),
  KW("treememlimit",	set_int, XPRS_TREEMEMORYLIMIT,	treememlimit_desc),
  KW("treememtarget",	set_dbl, XPRS_TREEMEMORYSAVINGTARGET, treememtarget_desc),
  KW("treeoutlev",	set_int, XPRS_TREEDIAGNOSTICS,	treeoutlev_desc),
  KW("varselection",	set_int, XPRS_VARSELECTION,	varselection_desc),
  KW("wantsol",		WS_val, 0,			WS_desc_ASL),
     };

static Option_Info Oinfo = { "xpress", NULL, "xpress_options",
           keywds,nkeywds,0,"XPRESS", 0,0,0,0,0, 20120417 };

 static int breaking;
 static jmp_buf Jb;

 static void
intcatch(int n)
{
	printf("\n<BREAK> (xpress)\n");
	fflush(stdout);
	if (!breaking++)
		XPRSinterrupt(prob, XPRS_STOP_CTRLC);
	else if (breaking > 3)
		longjmp(Jb, 2);
	}

/* Xpress-MP callback in case the user wants some output from Optimizer */
void XPRS_CC xpdisplay(XPRSprob prob, void *data, const char *ch, int n, int msglvl)
{
  /*
   msglvl gives the message level as follows:
   * 1 dialogue
   * 2 information
   * 3 warnings
   * 4 errors
   * a negative value indicates the XPRESS is about to finish and
   * buffers should be flushed.
   */

  /* You could divert the messages to your own log file if you wanted */

  if (msglvl < 0)
    fflush(NULL);
  else if (msglvl >= prtmsg && (msglvl != 4 || strncmp(ch, "?899 ", 5)))
    printf("%s\n", ch);
}

/***********************/
/* Print abort message */
/***********************/
static void xperror(const char *fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  fprintf(stderr, progname ? "%s: Error " : "Error ", progname);
  vfprintf(stderr, fmt, ap);
  fprintf(stderr, ".\n");
  va_end(ap);
  exit(1);
}

/******************************/
/* Delete .sol and .glb files */
/******************************/
static void killtempprob(void)
{
  int len = strlen(probname);
#if 0 /* would be needed if a call on XPRSwritebinsol were added */
  strcat(probname,".sol");
  remove(probname);
  probname[len]='\0';
#endif
  strcat(probname,".glb");
  remove(probname);
  probname[len]='\0';
}

static SufDecl
suftab[] = {
  { "direction", 0, ASL_Sufkind_var },
  { "priority", 0, ASL_Sufkind_var },
  { "ref", 0, ASL_Sufkind_var | ASL_Sufkind_real },
  { "sos", 0, ASL_Sufkind_var },
  { "sos", 0, ASL_Sufkind_con },
  { "sosno", 0, ASL_Sufkind_var | ASL_Sufkind_real },
  { "sosref", 0, ASL_Sufkind_var | ASL_Sufkind_real },
  { "sstatus", 0, ASL_Sufkind_var, 1 },
  { "sstatus", 0, ASL_Sufkind_con, 1 },
  { "unbdd", 0, ASL_Sufkind_var | ASL_Sufkind_real },
  };

/**********************/
/* The main procedure */
/**********************/

int main(int argc, char *argv[])
{
 char *stub;
 static char xpprompt[20];
 int iret, version;
 dims d;
 void (*oic)(int);

 Times[0] = xectim_();

 /* for debugging, allow -=, -v, -? to work */
 if (argc == 2 && *(stub = argv[1]) == '-' && stub[1] && !stub[2]) {
	asl = (ASL*)ASL_alloc(ASL_read_fg);
	getstub(&argv, &Oinfo);
	return 0;
	}

 iret = XPRSinit(XPRESS);
 if (iret)
   xperror("initialising Xpress-MP (return code %d)\n%s",
	iret, "Have you set XPRESS?/Is Flex running?");

 iret = XPRScreateprob(&prob);
 if (iret)
   xperror("Creating Xpress-MP problem(return code %d)\n%s",
	iret, "Have you set XPRESS?");

#ifdef UNIX
  XPRSsetlogfile(prob,"/dev/null");  /* Kill output from XPRSinit() */
#endif

 XPRSsetcbmessage(prob, xpdisplay, NULL);

 XPRSgetintcontrol(prob,XPRS_VERSION, &version);  /* Get the version number */
 sprintf(xpprompt,"XPRESS %2.2f",(double)version/100.0);/* set the banner */
 Oinfo.bsname = Oinfo.version = xpprompt;

 asl = (ASL*)ASL_alloc(ASL_read_fg);		/* Allocate a structure */
 if(!(stub = getstub(&argv, &Oinfo)))		/* Get the 'stub' name */
   usage_ASL(&Oinfo,1);

 if(tmpnam(probname)==NULL) xperror("tmpnam() failed: cannot obtain a problem name");
                                                 /* Get temporary problem name */

 suf_declare(suftab, sizeof(suftab)/sizeof(SufDecl));

 amplin(stub,argv,&d);        /* Read and load the problem */

#ifdef XPRS_SOLUTIONFILE
 if(optimopt[1]!='g')
   XPRSsetintcontrol(prob,XPRS_SOLUTIONFILE,0);  /* Don't save solution */
#endif


 if(startbasis!=NULL)                            /* Load an initial basis */
   if(XPRSreadbasis(prob,startbasis,""))
     xperror("loading an initial basis");

 Times[1] = xectim_();

 breaking = 0;
 iret = setjmp(Jb);
 oic = signal(SIGINT, intcatch);
 if (iret) {
	signal(SIGINT, oic);
	--iret;
	Times[2] = xectim_();
	if (amplflag | Oinfo.wantsol & 1) {
		if (solve_result_num <= 0)
			solve_result_num = 520;
		iret = 0;
		write_sol("Break!", 0, 0, &Oinfo);
		}
	else
		printf("Break received.\n");
	}
 else {
	iret = Optimise(prob,optimopt);                /* Optimise the problem */
	if (iret < 0 || (iret > 8 && iret != 32))
		xperror("optimising the problem:  surprise return %d", iret);
	Times[2] = xectim_();

	if((endbasis!=NULL)&&(optimopt[1]!='g')        /* Save the final basis */
	  && XPRSwritebasis(prob,endbasis,""))
		xperror("saving the current basis");

	amplout(&d);                                    /* Export the solution */
	}
 show_times();
 XPRSfree();
 return iret;
}

/*************************/
/* Set a known parameter */
/*************************/
static char *set_known(Option_Info *oi, keyword *kw, char *v)
{
#ifdef RWA_DEBUG
 char *d;
#endif
 switch((size_t) kw->info)
 {
  case set_primal:   optimopt[0]='p';  break;
#ifdef RWA_DEBUG
  case set_debug:  d=debugopt;
                        while( *v != '\0' && *v != ' ' ) *d++ = *v++;
      *d++ = '\0';
                        break;
#endif
  case set_dual:    optimopt[0]='d';       break;
  case set_barrier: optimopt[0]='b';       break;
  case set_relax:   optimopt[1]='l';       break;
  case set_maxim:   Optimise=XPRSmaxim;    break;
  case set_minim:   Optimise=XPRSminim;    break;
  case set_timing:  timing=1;              break;
  case set_iis:     iis_find=1;            break;
  default:          printf("Unknown option %s\n",kw->name);
                    badopt_ASL(oi);
 }
 return v;
}

/**************************************/
/* Set an Xpress-MP integer parameter */
/**************************************/

 static int
IntDset(Defer_setting *ds)
{
	if (XPRSsetintcontrol(prob,ds->ipar,ds->u.i)) {
		printf("The value %d could not be assigned belatedly to %s\n",
			ds->u.i, ds->kw->name);
		return 1;
		}
	return 0;
	}

static char *set_int(Option_Info *oi, keyword *kw, char *v)
{
  Defer_setting *ds;
  int ipar,isval;
  char *rv;

  ipar = (int) (unsigned long) kw->info;

  if ((*v=='?') && (v[1]<=' ')) {
    if (ipar < 0)
	ipar = -ipar;
    XPRSgetintcontrol(prob,ipar,&isval);
    printf("%s=%d\n",kw->name,isval);
    oi->option_echo &= ~ASL_OI_echothis;
    return v+1;
  }
  isval = (int)strtol(v,&rv,10);
  if (v==rv) {
    printf("Expected a numeric value for %s, not \"%s\"\n",kw->name,v);
    badopt_ASL(oi);
  } else if (ipar < 0) {
	ds = new_Defer_setting(IntDset, kw, ipar);
	ds->u.i = isval;
  } else if (XPRSsetintcontrol(prob,ipar,isval)) {
    printf("The value %d is not allowed for %s\n",isval,kw->name);
    badopt_ASL(oi);
  }

  return rv;
}

/*************************************/
/* Set an Xpress-MP double parameter */
/*************************************/

 static int
DblDset(Defer_setting *ds)
{
	if (XPRSsetdblcontrol(prob,ds->ipar,ds->u.d)) {
		printf("The value %g could not be assigned belatedly to %s\n",
			ds->u.d, ds->kw->name);
		return 1;
		}
	return 0;
	}

static char *set_dbl(Option_Info *oi, keyword *kw, char *v)
{
 Defer_setting *ds;
 int ipar;
 double dgval;
 char *rv;

 ipar= (int) (unsigned long) kw->info;
 if((*v=='?') && (v[1]<=' '))
 {
  if (ipar < 0)
	ipar = -ipar;
  XPRSgetdblcontrol(prob,ipar,&dgval);
  printf("%s=%g\n",kw->name,dgval);
  oi->option_echo &= ~ASL_OI_echothis;
  return v+1;
 }
 dgval= strtod(v,&rv);
 if(v==rv)
 {
  printf("Expected a numeric value for %s, not \"%s\"\n",kw->name,v);
  badopt_ASL(oi);
 }
 else if (ipar < 0) {
	ds = new_Defer_setting(DblDset, kw, ipar);
	ds->u.d = dgval;
	}
 else
  if(XPRSsetdblcontrol(prob,ipar,dgval))
  {
   printf("The value %g is not allowed for %s\n",dgval,kw->name);
   badopt_ASL(oi);
  }
 return rv;
}

/***********************************/
/* Set a filename option/parameter */
/***********************************/
static char *set_fln(Option_Info *oi, keyword *kw, char *v)
{
 char *rv, *t,**f;
 int n, q;

 if(!*v)
 {
  printf("rejecting %s: no following file name\n", kw->name);
  badopt_ASL(oi);
  return v;
 }
 f = (char**)kw->info;
 if (*v == '?' && v[1] <= ' ')
 {
  if (!(t = *f))
	t = "";
  printf("%s=\"%s\"\n", kw->name, t);
  oi->option_echo &= ~ASL_OI_echothis;
  return v + 1;
 }
 /* Allow optional quoting of file name with ' or " */
 /* Quoted file names may contain  blanks. */
 if ((q = *v) == '"' || q == '\'')
	for(rv = ++v; *rv != q && *rv; rv++);
 else
	for(rv = v; *++rv > ' ';);
 if (!(n = rv - v))
	t = 0;
 else {
	t = M1alloc(n + 1);
	strncpy(t, v, n);
	t[n] = '\0';
	}
 *f = t;
 if (q && *rv == q)
	rv++;
 return rv;
}

/*********************************/
/* Exit if problem is non linear */
/*********************************/

 static int
nlify(int k, const char *msg)
{
#ifndef NLIFY_LL
#define NLIFY_LL 72
#endif
	/* break error message into lines <= NLIFY_LL characters long */
	const char *s;
	int j;

	for(;;) {
		while(*msg <= ' ')
			if (!*msg++)
				return k;
		for(s = msg; *++s > ' ';);
		j = (int)(s - msg);
		if (j + k > NLIFY_LL)
			k = fprintf(Stderr, "\n%.*s", j, msg) - 1;
		else
			k += fprintf(Stderr, " %.*s", j, msg);
		msg += j;
		}
	}

 static void
nonlin(int n, char *what)
{
	if (n) {
		n = fprintf(Stderr, "Sorry, %s", Oinfo.bsname);
		n = nlify(n, "cannot handle");
		nlify(n, what);
		putc('\n', Stderr);
		exit(4);
		}
	}

/*****************************************/
/* Read priorities for integral entities */
/*****************************************/

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


 static int
priority_suf(void)
{
  int baddir, badpri, i, j, k, listsize, nnames;
  int *colindex, *dir, *p, *priority;
  char *direction;
  SufDesc *dp, *dd;

  baddir = badpri = i = listsize = 0;
  dd = suf_get("direction", ASL_Sufkind_var);
  dir = dd->u.i;
  dp = suf_get("priority", ASL_Sufkind_var);
  if (!(p = dp->u.i) && !dir)
    return 0;
  direction = 0;
  nnames = n_var;
  k = p ? 2 : 1;
  if (p && dir) {
    for(; i < nnames; i++)
      if (p[i] || dir[i])
        listsize++;
    }
  else if (p)
    listsize = nzeros(p, nnames);
  else
    listsize = nzeros(dir, nnames);
  if (!listsize)
    return 1;
  priority = 0;
  colindex = (int*)M1alloc(k*listsize*sizeof(int));
  if (dir)
    direction = (char*)M1alloc(listsize);
  i = k = 0;
  if (dir && p) {
    priority = colindex + listsize;
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
          j = 'D';
          break;
          case 1:
          j= 'U';
          break;
          default:
          baddir++;
          /* no break */
          case 0:
          j = 'N';
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
    for(; i < nnames; i++)
      if (dir[i]) {
        switch(dir[i]) {
          case -1:  j = 'D';      break;
          case 1:   j = 'U';      break;
          default:  baddir++;   /* no break */
          case 0:   j = 'N';
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
  if (XPRSloaddirs(prob,listsize,colindex,priority,direction,NULL,NULL))
    xperror("loading priorities");
  return 1;
  }

 static void
 mip_priorities(void)
{
 int nbpri, *start, *pri, *num;
 int ndir, *mcols, *mpri;
 int p, c, curr;

 if (priority_suf())
  return;
 if((nbpri=mip_pri(&start, &num, &pri, 2147483647))>0)
 {
  ndir=0;
  for(p=0;p<nbpri;p++) ndir+=num[p];  /* Number of priorities */
  mcols=(int *)M1alloc(ndir*sizeof(int));
  mpri=(int *)M1alloc(ndir*sizeof(int));

  curr=0;
  for(p=0;p<nbpri;p++)        /* Create mcols & mpri... */
   for(c=0;c<num[p];c++)
   {
    mcols[curr]=start[p]+c;
    mpri[curr]=pri[p];
    curr++;
   }
  if(XPRSloaddirs(prob,ndir,mcols,mpri,NULL,NULL,NULL))
    xperror("loading priorities");
 }
}

 static void
stat_map(int *stat, int n, int *map, int mx, char *what)
{
  int bad, i, i1=0, j, j1=0;
  static char badfmt[] = "XPRESS driver: %s[%d] = %d\n";

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
    "Xpress-MP driver: %d messages about bad %s values suppressed.\n",
        bad-1, what);
    }
  }

 static void
get_statuses(dims *d)
{
  static int map[] = {0, 1, 3, 0, 2, 0, 2};

  if ((!mipststat && niv + nbv)
   || (!(d->csd->kind & ASL_Sufkind_input)
   && !(d->rsd->kind & ASL_Sufkind_input)))
    return;
  stat_map(d->cstat, n_var, map, 7, "incoming cstat");
  stat_map(d->rstat, n_con, map, 7, "incoming rstat");
  if (XPRSloadbasis(prob,d->rstat, d->cstat))
    xperror("loading statuses");
  }

 static void
qcadj(int **pnelq, int **pqcrows, int **qcol1, int **qcol2, real **pqv)
{
	/* Adjust for quadratic constraints. */

	cgrad **cgp, *cg, *cga, *cg1, **cgt;
	char errbuf[256];
	double *a, *a1, *qv, *x, *y;
	fint *colqf, *cqf, *rowqf;
	int i, i1, i2, j, j1, je, k, kl, nl, nqc, m, n, nnz, nq, nqcnl, nqv;
	int *col0, *col1, *col2, *ia, *ia1, *ka, *ka1, *nelq, *qcrows;

	n = n_var;
	nqc = nlc;
	nqv = nlvc;
	ia = A_rownos;
	ka = A_colstarts;
	a = A_vals;

	j = nzc;
	for(i = kl = nl = 0; i < j; i++) {
		if (ia[i] < nqc)
			++nl;
		else
			++kl;
		}

	cg = cga = (cgrad*)Malloc((nqc+n)*sizeof(cgrad*) + nl*sizeof(cgrad));
	Cgrad = cgp = (cgrad**)(cga + nl);
	cgt = cgp + nqc;
	memset(cgp, 0, (nqc+n)*sizeof(cgrad*));

	for(i = j = 0; i < n; i++) {
		for(k = ka[i+1]; j < k; j++)
			if ((i1 = ia[j]) < nqc) {
				cg->varno = i;
				cg->coef = a[j];
				cg->next = cgp[i1];
				cgp[i1] = cg++;
				}
		}
	x = LUrhs;
	for(i = nqcnl = 0; i < nqc; i++) {
		i2 = i << 1;
		if (x[i2] > negInfinity && x[i2+1] < Infinity) {
			snprintf(errbuf, sizeof(errbuf),
			 "constraint %s, which is not convex quadratic since it is %s constraint.",
				con_name(i), x[i2] == x[i2+1] ? "an equality"
						: "a two-sided");
			nonlin(1, errbuf);
			}
		j = mqpcheck(-(i+1), 0, 0, 0);
		if (j < 0) {
			nonlin(j == -2,
			 "a quadratic constraint involving division by 0");
			nonlin(1, "a nonquadratic nonlinear constraint");
			}
		nonlin(j == 0,
		 "a driver bug: no quadratic terms in a \"nonlinear\" constraint");
		nqcnl += j;
		}

	x = *pqv = qv = (double*)Malloc(nqc*sizeof(int)
					+ nqcnl*(2*sizeof(int) + sizeof(double)));
	*pnelq = nelq = (int*)(qv + nqcnl);
	*qcol1 = col1 = nelq + nqc;
	*qcol2 = col2 = col1 + nqcnl;

	for(i = m = nl = 0; i < nqc; i++) {
		nq = mqpcheck(-(i+1), &rowqf, &colqf, &y);
		cqf = colqf;
		i1 = *++cqf;
		col0 = col1;
		/* Discard elements above the diagonal, and account for XPRESS's */
		/* expected scaling of quadratic constraints. */
		for(j = k = 0; j < nq; j++) {
			while (j >= i1) {
				k++;
				i1 = *++cqf;
				}
			if (rowqf[j] <= k) {
				*x++ = 0.5 * y[j];
				*col1++ = rowqf[j];
				*col2++ = k;
				}
			}
		nelq[i] = col1 - col0;
		free(colqf);
		free(rowqf);
		free(y);
		}

	nl = 0;
	while(i > 0) {
		for(cg = cgp[--i]; cg; cg = cg1) {
			cg1 = cg->next;
			if (cg->coef != 0.) {
				++nl;
				cg->next = cgt[j = cg->varno];
				cgt[j] = cg;
				cg->varno = i;
				}
			}
		}

	nnz = kl + nl + ka[n] - ka[nqv] + 1;
	A_vals = a1 = (real*)Malloc(nnz*sizeof(real) + (nqc + nnz + n + 2)*sizeof(int));
	A_rownos = ia1 = (int*)(a1 + nnz);
	A_colstarts = ka1 = ia1 + nnz;
	*pqcrows = qcrows = ka1 + n + 2;
	*ka1 = 0;
	for(i = i1 = j = j1 = 0; i < n; ++i) {
		for(cg = cgt[i]; cg; cg = cg->next) {
			ia1[j1] = cg->varno;
			a1[j1++] = cg->coef;
			}
		for(je = ka[i+1]; j < je; ++j)
			if (ia[j] >= nqc) {
				ia1[j1] = ia[j];
				a1[j1++] = a[j];
				}
		*++ka1 = j1;
		}
	for(i = 0; i < nqc; ++i)
		qcrows[i] = i;
	free(cga);
	free(a);
	}

 typedef struct
IndicInfo {
	int nic, nr;
	int *rn;
	int *indv;
	int *comps;
	} IndicInfo;

 static int
add_indic(void *v, int iv, int compl, int sense, int nz, int *ig, real *g, real rhs)
{
	IndicInfo *II = (IndicInfo*)v;
	int i, j, k, mstart[2];

	mstart[0] = 0;
	mstart[1] = nz;
	if (sense <= 1) {
		i = XPRSaddrows(prob, 1, nz, "LGE" + sense, &rhs, NULL, mstart, ig, g);
		k = 1;
		}
	else {
		/* Since XPRESS disallows equalities here, we supply two inequalities. */
		i = XPRSaddrows(prob, 1, nz, "L", &rhs, NULL, mstart, ig, g)
		    ||
		    XPRSaddrows(prob, 1, nz, "G", &rhs, NULL, mstart, ig, g);
		k = 2;
		}
	if (i)
		return i;
	for(j = 0; j < k; ++j) {
		i = II->nic++;
		II->rn[i] = II->nr++;
		II->indv[i] = iv;
		II->comps[i] = compl ? -1 : 1;
		}
	return 0;
	}

 static void
indicator_constrs(void)
{
	IndicInfo II;
	int errinfo[2], i, nlogc;

	nlogc = 4*n_lcon;	/* Allow space for "else" constraints and two */
				/* inequalities instead of an equality. */
	II.nic = 0;
	if (XPRSgetintattrib(prob, XPRS_ROWS, &II.nr))
		xperror("XPRSgetintattrib(XPRS_ROWS)");
	II.rn = (int*)Malloc(3*nlogc*sizeof(int));
	II.indv = II.rn + nlogc;
	II.comps  = II.indv + nlogc;
	if ((i = indicator_constrs_ASL(asl, &II, add_indic, errinfo))) {
		switch(i) {
		  case 1:
			xperror("logical constraint %s is not an indicator constraint.\n",
				lcon_name(errinfo[0]));
			break;
		  case 2:
			xperror("logical constraint %s is not an indicator constraint\n\
	due to bad comparison with %s.\n", lcon_name(errinfo[0]), var_name(errinfo[1]));
			break;
		  case 3:
			xperror("adding indicator row");
			break;
		  default:
			xperror("indicator_constrs_ASL");
		  }
		}
	if (XPRSsetindicators(prob, II.nic, II.rn, II.indv, II.comps))
		xperror("XPRSsetindicators");
	free(II.rn);
	}

/************************************************/
/* Read the matrix and call loadprob/loadglobal */
/************************************************/
 static void
amplin(char *stub, char *argv[], dims *d)
{
 FILE *nl;
 char *qgtype, *qrtype, *qstype;
 double *L, *U, *a, *dref, *obj, *q, *q0, *q1, *qe, *qv, *rhs, *rng;
 fint *colq, nelq, *rowq, *rq, *rq1;
 int *ia, *ka, *mgcols, *mqc1=NULL, *mqc2=NULL, *mscols, *msstart;
 int *qcol1, *qcol2, *qcrows, *qmn;
 int i, iret, j, m, m1, n, n0, n_bv, ngents, nq, nqc, nsets, nz, row;
 ograd *og;
 struct LU_bounds *rhs_bounds;
#ifndef NO_XPR_CON_INDICATOR
#define ALLOW_CLP ASL_allow_CLP
 int nlogc;
#else
#define ALLOW_CLP 0
#define nlogc 0
#endif

 nl = jac0dim(stub, (fint)strlen(stub));

 /* allow elbow room for objadj */
 m = n_con;
#ifndef NO_XPR_CON_INDICATOR
 nlogc = n_lcon;
#endif
 m1 = m + 4*nlogc + 1;	/* "+ 1" due to nextra == 1 for sstatus */
 n = n_var + 1;
 if (!m)
	m = 1;	/* we'll add a constraint to bypass a defect in XPRESS */
 nz = nzc + 1;
 d->cstat = (int*)M1zapalloc((m1 + n)*sizeof(int));
 d->rstat = d->cstat + n;
 d->csd = suf_iput("sstatus", ASL_Sufkind_var, d->cstat);
 d->rsd = suf_iput("sstatus", ASL_Sufkind_con, d->rstat);
 d->miqp = 0;

 obj = (real*)Malloc((5*m1+4*n)*sizeof(real) + m*sizeof(char));
 rng = obj + n;
 memset(obj, 0, (n+m1)*sizeof(real));
 LUv = rng + m1;
 Uvx = LUv + n;
 LUrhs = Uvx + n;
 rhs = LUrhs + 2*m1;
 d->x = rhs + m1;
 d->y = d->x + n; memset(d->y, 0, sizeof(real) * m1);
 qrtype = (char*)(d->y + m1);

 A_vals = a = (real*)Malloc(nz*sizeof(real)+(nz+n+1)*sizeof(int));
 A_rownos = ia = (int*)(a + nz);
 A_colstarts = (int*)(ia + nz);

 want_deriv = 0;
 qp_read(nl,ALLOW_CLP);
 if (!n_obj)
  nobj = 0;
 if(getopts(argv, &Oinfo)) exit(1);     /* Set options */
 if(logfile != NULL)
 {
  if(XPRSsetlogfile(prob,logfile))
    xperror("opening the logfile");
  XPRSsetintcontrol(prob,XPRS_OUTPUTLOG,1);        /* Chat mode */
 }
 qstype=NULL;
 msstart=NULL;
 mscols=NULL;
 qv = dref = NULL;

 i = sos ? 0 : ASL_suf_sos_ignore_sosno;
 if (!sos2)
  i |= ASL_suf_sos_ignore_amplsos;
 nsets = suf_sos(i, 0, &qstype, 0,0, &msstart, &mscols, &dref);
 ngents = niv + nbv + nlvbi + nlvci + nlvoi;
 if ((ngents>0) && (optimopt[1]=='l')) {
	printf("Ignoring integrality of %d variable%s.\n",
	    ngents, ngents > 1 ? "s" : "");
	nlvoi = nbv = niv = ngents = 0;
	}
 m = n_con;
 n = n_var;

                        /* Preparing the objective function */
 obj_no = --nobj;
 og = 0;
 nq = nelq = 0;
 q0 = 0;
 if(nobj < n_obj && nobj >= 0) {
	if ((nelq = mqpcheck(nobj, &rowq, &colq, &q0))) {
		if (nelq < 0) {
			nonlin(nelq == -2,
			 "a quadratic objective involving division by 0");
			nonlin(1, "a non-quadratic nonlinear objective");
			}
		d->miqp = ngents;
		/* discard quadratic terms below the diagonal */
		q = q1 = q0;
		rq = rq1 = rowq;
		for(i = 1; i <= n; i++) {
			for(qe = q0 + colq[i]; q < qe; q++)
				if ((*rq1 = *rq++) < i) {
					rq1++;
					*q1++ = *q;
					}
			colq[i] = q1 - q0;
			}
		nq = colq[n];
		mqc2 = (int*)Malloc(nq*sizeof(int));
		for(rq = rowq, i = j = 0; i < n; i++)
			for(rq1 = rowq + colq[i+1]; rq < rq1; rq++)
				mqc2[j++] = i;
		free(colq);
		if (sizeof(int) == sizeof(fint)) /* most likely */
			mqc1 = (int*)rowq;
		else {
			mqc1 = (int*)Malloc(nq*sizeof(int));
			for(i = 0; i < nq; i++)
				mqc1[i] = rowq[i];
			free(rowq);
			}
		}
	og = Ograd[nobj];
	objadj = objconst(nobj);
	if(Optimise==NULL)
	Optimise = objtype[nobj] == 0 ? XPRSminim : XPRSmaxim;
	}
  else if (nobj != -1) {
	fprintf(Stderr,"Objective %d does not exist.\n",nobj+1);
	exit(1);
	}
  else if (!Optimise)
	Optimise = XPRSminim;

  qmn = qcrows = qcol1 = qcol2 = 0;
  if ((nqc = nlc))
	qcadj(&qmn, &qcrows, &qcol1, &qcol2, &qv);
  a = A_vals; /* qcadj may have changed A_vals, A_rownos and A_colstarts */
  ia = A_rownos;
  ka = A_colstarts;

  for(; og; og = og->next)  /* The coefficients */
    obj[og->varno] = og->coef;

  n0 = n;
  if(objadj) {      /* Getting adjustment value */
                    /* If necessary, add a new variable */
	LUv[n]=objadj;
	Uvx[n]=objadj;
	obj[n]=1.0;
	n++;                    /* pretend 1 more column */
	A_colstarts[n]=A_colstarts[n0];
	}

      /* Make RHS, QRTYPE and RNG */

 rhs_bounds=(struct LU_bounds *)LUrhs;

 for(row=0;row<m;row++)
  if(rhs_bounds[row].upper==Infinity)
  {
   if(rhs_bounds[row].lower==negInfinity)
   {
    qrtype[row]='N';        /* Non-binding constraint */
    rhs[row]=0;
   }
   else
    {
     qrtype[row]='G';       /* >= constraint */
     rhs[row]=rhs_bounds[row].lower;
    }
  }
  else
   if(rhs_bounds[row].lower==rhs_bounds[row].upper)
    {
     qrtype[row]='E';       /* == constraint */
     rhs[row]=rhs_bounds[row].lower;
    }
   else
    if(rhs_bounds[row].lower==negInfinity)
     {
      qrtype[row]='L';      /* <= constraint */
      rhs[row]=rhs_bounds[row].upper;
     }
    else
     {
      qrtype[row]='R';      /* Range constraint */
      rhs[row]=rhs_bounds[row].upper;
      rng[row]=rhs_bounds[row].upper-rhs_bounds[row].lower;
     }

 mgcols = 0;
 if(ngents + nsets == 0)    /* Is it just LP, or MIP ? */
 {
  iret = nqc
	? XPRSloadqcqp(prob,probname,n,m,qrtype,rhs,rng,obj,
		ka,NULL,ia,a,LUv,Uvx,
		nq, mqc1, mqc2, q0,
		nqc,qcrows,qmn,qcol1,qcol2,qv)
	: nq
	? XPRSloadqp(prob,probname,n,m,qrtype,rhs,rng,obj,
		ka,NULL,ia,a,LUv,Uvx,
		nq, mqc1, mqc2, q0)
	: XPRSloadlp(prob,probname,n,m,qrtype,rhs,rng,obj,
		ka,NULL,ia,a,LUv,Uvx);
  if(iret)
    xperror("loading the problem");
 }
 else
 {
  qgtype = 0;
  if (ngents) {
	mgcols = (int *) Malloc(ngents*(sizeof(int)+1));
	qgtype = (char *)(mgcols + ngents);

	L = LUv;
	U = Uvx;
	j = 0;
	n_bv = nlvb;
	for(i = nlvb - nlvbi; i < n_bv; ++i, ++j) {	/* nonlinear integer or binary variables */
						/* in both constraints and objectives */
		qgtype[j] = L[i] == 0. && U[i] == 1. ? 'B' : 'I';
		mgcols[j] = i;
		if (L[i] < MININT)
			L[i] = MININT;
		if (U[i] > XPRS_MAXINT)
			U[i] = XPRS_MAXINT;
		}
	n_bv = nlvc;
	for(i = n_bv - nlvci; i < n_bv; ++i, ++j) {	/* nonlinear integer or binary variables */
							/* just in constraints */
		qgtype[j] = L[i] == 0. && U[i] == 1. ? 'B' : 'I';
		mgcols[j] = i;
		if (L[i] < MININT)
			L[i] = MININT;
		if (U[i] > XPRS_MAXINT)
			U[i] = XPRS_MAXINT;
		}
	n_bv += nlvo - nlvc;
	for(i = n_bv - nlvoi; i < n_bv; ++i, ++j) {	/* nonlinear integer or binary variables */
							/* just in objectives */
		qgtype[j] = L[i] == 0. && U[i] == 1. ? 'B' : 'I';
		mgcols[j] = i;
		if (L[i] < MININT)
			L[i] = MININT;
		if (U[i] > XPRS_MAXINT)
			U[i] = XPRS_MAXINT;
		}

	n_bv = n_var;
	for(i = n_bv - (nbv + niv); i < n_bv; ++i, ++j) { /* linear integer or binary variables */
		qgtype[j] = L[i] == 0. && U[i] == 1. ? 'B' : 'I';
		mgcols[j] = i;
		if (L[i] < MININT)
			L[i] = MININT;
		if (U[i] > XPRS_MAXINT)
			U[i] = XPRS_MAXINT;
		}
	}

  if( !(m) ) { /* code around XPRESS defect: add a nonbinding constraint */
	m = 1;
	A_vals[0] = 1.;
	qrtype[0] = 'N';
	rhs[0] = 0.;
	}

  iret = nqc
	? XPRSloadqcqpglobal(prob,probname,n,m,qrtype,rhs,rng,obj,
		ka,NULL,ia,a,LUv,Uvx,
		nq, mqc1, mqc2, q0,
		nqc,qcrows,qmn,qcol1,qcol2,qv,
		ngents,nsets,qgtype,mgcols,NULL/*mplim*/,
		qstype,msstart,mscols,dref)
	: nq
	? XPRSloadqglobal(prob,probname,n,m,qrtype,rhs,rng,obj,
         	ka,NULL,ia,a,LUv,Uvx,
		nq, mqc1, mqc2, q0,
		ngents,nsets,qgtype,mgcols,NULL/*mplim*/,
		qstype,msstart,mscols,dref)
	: XPRSloadglobal(prob,probname,n,m,qrtype,rhs,rng,obj,
		ka,NULL,ia,a,LUv,Uvx,
		ngents,nsets,qgtype,mgcols,NULL/*mplim*/,
		qstype,msstart,mscols,dref);
  if(iret) xperror("loading the problem");
#ifndef NO_XPR_CON_INDICATOR
	if (nlogc)
		indicator_constrs();
#endif /*NO_XPR_CON_INDICATOR*/

  if(optimopt[1]=='\0')
  {
   optimopt[1]='g';       /* Search will be global */
   mip_priorities();      /* using provided priorities */
  }
 }
 if (mgcols)
	free(mgcols);
 if (qv)
	free(qv);
 if (q0) {
	free(mqc2);
	free(mqc1);
	free(q0);
	}
 free(a);

 Do_Defer();

#ifdef RWA_DEBUG
 if( strstr(debugopt,"save") ) XPRSsave(prob); /* save matrix to internal file */
#endif
 atexit(killtempprob);    /* Ensure temp files are removed on exit */
 get_statuses(d);
}

 static int
send_statuses(dims *d)
{
  int *cstat, *rstat;
  static int map[] = {3, 1, 4, 2};

  if (!(asl->i.flags & 1) && !amplflag)
    return 0;
  if (d->miqp)
	return 1; /* avoid unsuppressable message from XPRSgetbasis */
  cstat = d->cstat;
  rstat = d->rstat;
  memset(cstat, 0, n_var*sizeof(int));
  memset(rstat, 0, n_con*sizeof(int));
  if (XPRSgetbasis(prob,rstat, cstat))
	return 1;
  stat_map(cstat, n_var, map, 3, "outgoing cstat");
  stat_map(rstat, n_con, map, 3, "outgoing rstat");
  return 0;
}

 static int
getvec(int *mrow, double *dmat, int mxelt, int *pnelt, int jvec)
{
  int nrow;
  char qrtype;          /* Row type for slack variable */
  static int mbeg[2];

  XPRSgetintattrib(prob,XPRS_ROWS, &nrow);
  if (jvec < nrow) {   /* Vector is a slack/surplus */
    *pnelt = 1;
    if (mxelt < 1) return 0;
  if (XPRSgetrowtype(prob,&qrtype, jvec, jvec)) return 1;

    *mrow = jvec; *dmat = (qrtype == 'G') ? -1.0 : 1.0;
  }
  else {                /* Vector is a structural */
    if (XPRSgetcols(prob,mbeg, mrow, dmat, mxelt, pnelt, jvec - nrow, jvec - nrow))
      return 1;
  }

  return 0;
}

 void
unpack(int *mind, double *dnz, int size, int nnz, double *dvec)
{
  double *d = dvec;

  while (size-- > 0)
    *(d++) = 0.0;
  while (nnz-- > 0)
    dvec[*(mind++)] = *(dnz++);
}

 static int
send_ray(dims *d, char *hbuf, int *lenp)
{
  int *cstat, *mrow, *pivrow, *rstat, i, j, junb, nb, ncol, nelt, nrow, nrseq, ns;
  double *dmat, *dvec, *unbdd, dscale;

  int pstat,lstat;

  if (optimopt[1] == 'g') {
	if (!Ray)
		Ray = 3; /* don't say "not requested" */
	return 1;
	}
  switch(Ray) {
	default: return 1;
	case 2:
		if (optimopt[0] != 'p')
			goto use_primal;
		break;
	case 1:
		XPRSgetintcontrol(prob, XPRS_PRESOLVE, &i);
		if (!i)
			break;
		XPRSsetintcontrol(prob, XPRS_PRESOLVE, 0);
 use_primal:
		j = optimopt[0];
		optimopt[0] = 'p';
		if((*Optimise)(prob,optimopt))
			xperror("optimising the problem in send_ray");
		optimopt[0] = j;
		XPRSgetintattrib(prob, XPRS_LPSTATUS, &i);
		if (i != LPSTAT_UNBOUNDED)
			*lenp += Sprintf(hbuf + *lenp,
				"\nSurprise LPSTATUS = %d computing .unbdd", i);
		XPRSgetintattrib(prob, XPRS_BARITER, &nb);
		XPRSgetintattrib(prob, XPRS_SIMPLEXITER, &ns);
		if (nb)
			*lenp += Sprintf(hbuf + *lenp,
			  "\n%d extra barrier iteations computing .unbdd", nb);
		if (ns)
			*lenp += Sprintf(hbuf + *lenp,
			  "\n%d extra simplex iteations computing .unbdd", ns);
	}
  XPRSgetintattrib(prob,XPRS_PRESOLVESTATE, &pstat);
  XPRSgetintattrib(prob,XPRS_LPSTATUS, &lstat);

  XPRSgetintattrib(prob,XPRS_ROWS, &nrow);
  XPRSgetintattrib(prob,XPRS_COLS, &ncol);
  XPRSgetintattrib(prob,XPRS_SPAREROWS, &nrseq);
  nrseq += nrow;

  dmat  = (double*) M1alloc((nrseq+ncol+2*nrow)*sizeof(int) + 2*nrow*sizeof(double));
  dvec  = dmat + nrow;
  rstat = (int*)(dvec + nrow);
  cstat = rstat + nrseq;
  pivrow= cstat + ncol;
  mrow  = pivrow + nrow;

  i = (n_var > ncol) ? n_var : ncol;
  unbdd = (double*)M1zapalloc(i*sizeof(double));

  if (XPRSgetunbvec(prob,&junb))
    return 1;

  if (getvec(mrow, dmat, nrow, &nelt, junb))
    return 1;

  unpack(mrow, dmat, nrow, nelt, dvec);

  if (XPRSftran(prob,dvec))
    return 1;

  if (XPRSgetbasis(prob,rstat, cstat))
    return 1;

  if (XPRSgetpivotorder(prob,pivrow))
    return 1;

  dscale = (rstat[junb] == XP_NBASUP) ? -1 : 1;

  if (junb >= nrow) unbdd[junb-nrow] = dscale;

  dscale = -dscale;

  for (i = 0; i < nrow; i++){
    j = pivrow[i];
    if (j < nrseq)
  continue; /* it's a row that's basic */
    j -=nrseq;           /* get col seq number starting at 0 */
    if (cstat[j] == XP_BASIC)
  unbdd[j] = dscale * dvec[i];
  }

  suf_rput("unbdd", ASL_Sufkind_var, unbdd);

  return 0;
}


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

/***************************************************************/
/* Check the Xpress-MP status parameter and write the solution */
/***************************************************************/
static void amplout(dims *d)
{
 char hbuf[640], buf[32], *wb;
 double objvalue, w;
 int *cstatus, *rstatus;
 int didbarrier, i, ipstat, len, lpstat=-1, m, n;
 int nbit = 0, nbs, ncol, nint, nround, nrow, nsit = 0;
 real *x, *x1, *y;
 typedef struct { char *msg; int code; } Sol_info;
 static Sol_info report[]={
    { "Problem has not been loaded", 500 },
    { "Optimal solution found", 000 },
    { "Infeasible problem", 200 },
    { "Objective is worse than cutoff", 100 },
    { "Unfinished optimisation", 400 },
    { "Unbounded problem", 300 },
    { "Cutoff in dual", 101 },
    { "Problem unsolved", 502 }, /* should not happen */
    { "Problem is not convex", 510 }
    };
 static Sol_info repglb[]={
  { "Problem has not been loaded", 500},
  { "LP has not been optimized (probably LP Infeasible)", 501},
  { "LP has been optimised", 001},
  { "Global search incomplete - no integer solution found", 401},
  { "Global search incomplete", 102},
  { "Global search complete - no integer solution found", 201},
  { "Global search complete", 002},
  { "Unbounded problem with some integer variables", 301}
  };

 m = n_con;
 n = n_var;
 x = d->x;
 y = d->y;

 len=sprintf(hbuf,"%s: ",Oinfo.bsname);
 didbarrier = (optimopt[0]=='b');
 if(optimopt[1]=='g')    /* We did a Global search */
 {
  XPRSgetintattrib(prob,XPRS_MIPSTATUS, &ipstat);
  if (ipstat <= 2)  /* ...but never started global, so .sol file not created */
  {

#ifdef XPRS_SOLUTIONFILE
   XPRSsetintcontrol(prob,XPRS_SOLUTIONFILE,0);/*  In case we try to get the soln anyway:
                             prevent looking for absent .sol file */
               /* changed from seticv(N_IFMEM,1|4|8|16)*/
#endif
  }
 }
 if(optimopt[1]=='g' && ipstat > GLSTAT_LP_FINISHED)  /* We have a valid IP sol on the .sol file */
 {
  len += Sprintf(hbuf+len,"%s",repglb[ipstat].msg);
  solve_result_num = repglb[ipstat].code;
  switch (ipstat)
  {
   case GLSTAT_UNFINISHED_NOSOL:
   case GLSTAT_FINISHED_NOSOL:
           XPRSgetdblattrib(prob,XPRS_BESTBOUND,&objvalue);
           g_fmtop(buf,objvalue);
           len+=Sprintf(hbuf+len,"\nBest bound determined so far %s",buf);
	   /* no break */
   case GLSTAT_MIP_UNBOUNDED:
           x=y=NULL;
           break;
   case GLSTAT_UNFINISHED_SOL:
   case GLSTAT_FINISHED_SOL:
           XPRSgetdblattrib(prob,XPRS_MIPOBJVAL,&objvalue);
           g_fmtop(buf,objvalue);
           len += Sprintf(hbuf+len,"\nBest integer solution found %s",buf);
           XPRSgetintattrib(prob,XPRS_MIPSOLS,&nbs);
           if(nbs>1)      /* At least 2 solutions here */
              len+=Sprintf(hbuf+len,
                            "\n%d integer solutions have been found",nbs);
           if(XPRSgetmipsol(prob,x,NULL)) /* there are no dual variables for mip solutions */
           xperror("preparing solution file");
           break;
   default:
           xperror("unrecognised global status");
  }
  XPRSgetintattrib(prob,XPRS_NODES,&nbit);
  len += Sprintf(hbuf+len,"\n%d branch and bound node%s", nbit, (nbit!=1) ? "s" : "");
 }
 else        /* Just LP minim or maxim - even if prob was global */
 {
  XPRSgetintattrib(prob,XPRS_LPSTATUS, &lpstat);

  len += Sprintf(hbuf+len,"%s",report[lpstat].msg);
  solve_result_num = report[lpstat].code;
  if(lpstat==LPSTAT_OPTIMAL)
  {
   XPRSgetdblattrib(prob,XPRS_LPOBJVAL, &objvalue);
   g_fmtop(buf,objvalue);
   len += Sprintf(hbuf+len,"\nObjective %s",buf);
  }
  XPRSgetintattrib(prob,XPRS_BARITER,&nbit);
  XPRSgetintattrib(prob,XPRS_SIMPLEXITER,&nsit);

  if (didbarrier || nbit > 0) {
	i = -1;
	XPRSgetintcontrol(prob,XPRS_BARTHREADS,&i);
	if (i > 1)
		len += Sprintf(hbuf+len, "\n%d processors used.", i);
	}

  /* Get IIS */
  if(lpstat==LPSTAT_INFEASIBLE && iis_find==1) {
#if XPVERSION >= 21
    XPRSiisfirst(prob,1,&i);
#else
    XPRSiis(prob, "");
#endif
  }


/* If (lpstat==LPSTAT_INFEASIBLE or lpstat==LPSTAT_UNBOUNDED) and 'its'==0
   then no solution is available (presolve has proven infeasibility or
   unboundedness). */

  if (lpstat==LPSTAT_NONCONVEX)
	x = y = NULL;
  else if (nbit > 0
   || nsit > 0
   || !(lpstat==LPSTAT_INFEASIBLE
   || lpstat==LPSTAT_UNBOUNDED
   || lpstat==LPSTAT_NONCONVEX)) {
    if(XPRSgetlpsol(prob,x,NULL,y,NULL))
      xperror("preparing solution file");
  }
  else {

    /* there will be no basis, so set an all slack one */

    rstatus = (int *) d->y;
    cstatus = (int *) d->x;
    XPRSgetintattrib(prob,XPRS_ROWS,&nrow);
    XPRSgetintattrib(prob,XPRS_COLS,&ncol);
    for(i=0;i<nrow;i++) rstatus[i]=1;
    for(i=0;i<ncol;i++) cstatus[i]=0;
    XPRSloadbasis(prob,rstatus,cstatus);
    x = y = NULL;
  }
  if (nsit > 0)
	len += Sprintf(hbuf + len, "\n%d simplex iteration%s", nsit,
			"s" + (nsit == 1));
  if (nbit > 0)
	len += Sprintf(hbuf + len, "\n%d barrier iteration%s", nbit,
			"s" + (nbit == 1));
 }
 if ((nbit > 0 && !nsit) || send_statuses(d)) {
	d->csd->kind &= ~ASL_Sufkind_output;
	d->rsd->kind &= ~ASL_Sufkind_output;
	len += Sprintf(hbuf+len, "\nNo basis.");
	}
 if (lpstat == LPSTAT_UNBOUNDED && send_ray(d,hbuf,&len))
  len += Sprintf(hbuf+len, "\nNo unbounded vector%s.", Ray ? "" : " requested");


 /* make sure integer variables near to integer values are integers */

 if (niv + nbv + nlvoi && x && optimopt[1] != 'l' && Round >= 0) {
	nround = 0;
	w = 0;
	if ((nint = niv + nbv)) {
		x1 = x + n - nint;
		nround = xround(x1, nint, Round & 1, &w);
		}
	if ((nint = nlvoi)) {
		x1 = x + (nlvo - nint);
		nround += xround(x1, nint, Round & 1, &w);
		}
	if (w < 1e-9 && !(Round & 8))
		nround = 0;
	else if (nround) {
		if (solve_result_num < 200 && !(Round & 2))
			solve_result_num += 10;
		if (Round & 4)
			nround = 0;
		else if (!(Round & 1))
			nround = -nround;
		}
	if (nround) {
		wb = "";
		if (nround < 0) {
			nround = -nround;
			wb = "would be ";
			}
		len += Sprintf(hbuf+len,
			"\n%d integer variables %srounded (maxerr = %g).\n%s\n",
			nround, wb, w,
			"Reducing miptol to something < maxerr might help.");
		}
	}

 write_sol(hbuf,x,y,&Oinfo);
}

/***************************/
/* Some timing information */
/***************************/
static void show_times(void)
{
 Times[3] = xectim_();
 if(timing)
  printf("\nTimes (seconds):\nInput =  %g\nSolve =  %g\nOutput = %g\n",
      Times[1] - Times[0], Times[2] - Times[1],
      Times[3] - Times[2]);
}
