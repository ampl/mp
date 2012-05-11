#include "MCP_Interface.h"
#include "Path.h"
#include "PathOptions.h"
#include "Macros.h"
#include "Output.h"
#include "Options.h"

#undef Char
#include "getstub.h"

#define CRIPPLE_SIZE 200
#define CRIPPLE_NNZ 10000

extern char *Path_Version(void);
static int debug, functimes, quitnow, wantfuncs;
static int sideineq = 1, sqwarn = 1;
#ifndef LOGFILE
#define LOGFILE 0
#endif
#ifndef STATUSFILE
#define STATUSFILE 0
#endif
#ifndef OPTFILE
#define OPTFILE 0
#endif
static char *logfile = LOGFILE, *statusfile = STATUSFILE, *optfile = OPTFILE;
#ifdef MDEBUG /* debugging malloc */
extern int zapvalue;
#endif

 typedef struct
FuncInfo {
	real	*c_scratch;
	char	*ctype;
	char	*vtype;
	int	*vb;	/* variables with bounds */
	int	*sba;	/* constraints tweaked by jadj */
	int	*sbae;	/* end of sba */
	int	*cb;	/* constraints involving inequalities */
	int	*cbe;	/* end of cb */
	int	n;	/* adjusted problem size */
	int	ineq;	/* number of (simple) side inequalities */
	int	jadj;	/* whether to convert side inequalities from */
			/* expr >= 0 complements newvar >= 0  to */
			/* expr + newvar >= 0 complements newvar >= 0 */
	int	nvb;	/* number of bounded variables */
	int	nvr;	/* number of ranged variables */
	int	nnz;	/* adjusted number of Jacobian nonzeros */
	int	ncsq;	/* number of constraints added to make square */
	int	nvsq;	/* number of variables added to make square */
	} FuncInfo;

 static FuncInfo FI;

 static int numJacEvals;
 static Options_Interface *opt;
 static void install_interface(MCP *m);

 static fint
mkey(char *phrase, ftnlen len)
{
	Not_Used(len);

	Options_Set(opt, phrase);
	return 0;
	}

 static char *
Optfile(Option_Info *oi, keyword *kw, char *value)
{
	char *rv = C_val(oi, kw, value);
	if (optfile)
		Options_Read(opt, optfile);
	return rv;
	}

 static char
debug_desc[] = "debug level (default 0): sum of\n\
		1 = show initial z and bounds\n\
		2 = show z and f(z)\n\
		4 = show z and J = f'(z)\n\
		8 = show sparsity of J";

 static char
sideineq_text[] = "handling of side inequalities:\n\
		0 = no warning\n\
		1 = warn (default)\n\
		2 = quit, permitting AMPL loops to test solve_result\n\
		3 = quit, terminating AMPL loops and scripts\n\
		4 = do not make Jacobian nonsingular (avoid nonuniqueness)\n\
		5 = warn and do not make Jacobian nonsingular";

 static char
sqwarn_text[] = "whether to warn of nonsquare systems:\n\
		0 = no\n\
		1 = yes (default)\n\
		2 = quit, permitting AMPL loops to test solve_result\n\
		3 = quit, terminating AMPL loops and scripts";

 static keyword
keywds[] = {
	KW("debug", I_val, &debug, debug_desc),
#ifndef NEOS
	KW("logfile", C_val, &logfile, "name of log file"),
#endif
	KW("maxfwd",  IA_val, voffset_of(ASL, p.maxfwd_), "max vars in fwd AD of common exprs (default 5)"),
#ifndef NEOS
	KW("optfile", Optfile, &optfile, "name of options file"),
#endif
	KW("sideineq", I_val, &sideineq, sideineq_text),
	KW("sqwarn", I_val, &sqwarn, sqwarn_text),
#ifndef NEOS
	KW("statusfile", C_val, &statusfile, "name of status file"),
#endif
	KW("version", Ver_val, 0, "report version"),
	KW("wantsol", WS_val, 0, WS_desc_ASL + 5)
#ifdef MDEBUG
	,KW("zapvalue", I_val, &zapvalue, "zapvalue (in debugging malloc)")
#endif
	};

 static
keyword cloptions[] = {
	KW("f", IK1_val, &wantfuncs, "list available user-defined functions"),
	KW("q", IK1_val, &quitnow,   "quit now"),
	KW("t", IK1_val, &functimes, "time function evaluations")
	};

 static char xxxvers[] = "\nAMPL/PATH Driver Version 20020506\n";

 static
Option_Info Oinfo = {
	"path", 0, "path_options", keywds, nkeywds, 3, xxxvers,
	0, mkey, 0, cloptions, sizeof(cloptions) / sizeof(keyword),
	20020506 };

 static void
show_ver(FILE *f)
{
	fprintf(f, "%s", Oinfo.version);
	fflush(f);
	}

 static int
cvarcomp(const void *a, const void *b, void *v)
{
	return ((int*)v)[*(int*)a] - ((int*)v)[*(int*)b];
	}

enum {
	C_COMPL = 0,
	C_LB = 1,
	C_UB = 2,
	C_RANGE = 3,
	C_EQN = 4,
	C_FREE = 5,
	C_ZERO = 6
	};

 static int
mcp_adj(ASL *asl)
{
	int *Cvar, *z, *zc;
	int ct, i, j, k, *ka, L, nc, ncs, nv, nvb, nvr, nvs;
	real *lu;
	cgrad *cg;
	size_t st, st1;

	/* Match variables to equations and inequalities, adding one
	 * new variable for each simple inequality, including simple
	 * variable bounds, and two new variables for ranges:
	 * range constraints, -Infinity < L <= f(x) <= U < Infinity,
	 * with L < U, are represented by
	 *	f(x) = y complements v
	 *	v	 complements L <= y <= U
	 * in which v and y are new variables.  Ranged variables are
	 * handled similarly, with f(x) := x[i] for variable i if it
	 * has distinct finite lower and upper bounds.   This adds
	 * a new nonzero to column i of the Jacobian matrix.
	 */

	nc = n_con;
	ncs = nc + FI.ncsq;
	nv = n_var;
	nvs = nv + FI.nvsq;
	if (!(Cvar = cvar) || FI.ncsq) {
		Cvar = (int * )M1zapalloc(ncs * sizeof(int));
		if (cvar)
			memcpy(Cvar, cvar, nc*sizeof(int));
		cvar = Cvar;
		}
	st = ncs * sizeof(real);
	st1 = (ncs + nvs) * sizeof(int);
	if (st < st1)
		st = st1;
	FI.c_scratch = (real *)M1zapalloc(st);
	z = (int *)FI.c_scratch;
	zc = z + nvs;
	for (i = 0; i < nc; )
		if (j = Cvar[i++])
			z[j-1] = i;

	lu = LUv;

	/* nvb := number of variables with one or two finite bounds. */
	/* nvr := number of ranged variables: two finite bounds */

	for (i = nvb = nvr = 0; i < nv; i++, lu += 2)
		if (!z[i]) {
			if (lu[0] > negInfinity) {
				nvb++;
				if (lu[1] < Infinity)
					nvr++;
				}
			else if (lu[1] < Infinity)
				nvb++;
			}
	FI.nvb = nvb;
	FI.nvr = nvr;
	FI.vb = (int *)M1zapalloc((sizeof(int)+1)*nvb + ncs);
	FI.vtype = (char * )(FI.vb + nvb);
	FI.ctype = FI.vtype + nvb;
	if (nvb) {
		lu = LUv;
		for (i = k = 0; i < nv; i++, lu += 2)
			if (!z[i]) {
				ct = 0;
				if (lu[0] > negInfinity)
					ct = lu[1] < Infinity ? C_RANGE : C_LB;
				else if (lu[1] < Infinity)
					ct = C_UB;
				else
					continue;
				FI.vtype[k] = ct;
				FI.vb[k++] = i;
				}
		}

	/* match equations */

	j = nvs;
	for (i = k = 0; i < nc; i++) {
		if (!Cvar[i]) {
			lu = LUrhs + (i << 1);
			if (lu[0] < lu[1]) {
				Cvar[i] = ++j;

				/* For range constraints, */
				/* increment j twice... */

				if (lu[0] > negInfinity)
					FI.ctype[i] = lu[1] < Infinity
					     ? (j++,C_RANGE) : C_LB;
				else
					FI.ctype[i] = lu[1] < Infinity
					     ? C_UB : C_FREE;
				}
			else {
				FI.ctype[i] = C_EQN;
				while(z[k])
					k++;
				z[k++] = i + 1;
				Cvar[i] = k;
				}
			}
		}
	while(i < ncs) {
		FI.ctype[i] = C_ZERO;
		while(z[k])
			k++;
		Cvar[i++] = ++k;
		}

	/* Prepare zc = permutation of constraints. */

	for (i = 0; i < ncs; i++)
		--Cvar[zc[i] = i];
	if (nc > 1)
		qsortv(zc, nc, sizeof(int), cvarcomp, Cvar);

	/* Adjust A_colstarts to insert new singleton Jacobian rows */
	/* for ranged variables. */

	ka = A_colstarts + 1;
	j = k = 0;
	L = FI.nvb ? FI.vb[0] : nv;
	for (i = 0; i < nv; i++) {
		z[i] = ka[i] + j - 1;
		if (i == L) {
			j++;
			L = ++k < FI.nvb ? FI.vb[k] : nv;
			}
		ka[i] += j;
		}

	/* Adjust cg->goff for use in Jacobian computations */
	/* -- see sparseJ() below. */

	for (i = nc; i-- > 0; ) {
		j = zc[i];
		for (cg = Cgrad[j]; cg; cg = cg->next)
			cg->goff = z[cg->varno]--;
		}

	FI.nnz = nzc + 2*(nranges + nvr) + nvb;
	return ncs + nranges + nvb + nvr;
	}

 static void
prepare_jadj(ASL *asl)
{
	int *Cvar, i, j, n, ncb, nvb, *x, *y;

/* Adjust single side inequalities to make their contribution
 * to the Jacobian matrix nonsingular -- suggested by
 * Todd Munson (May 2000):  change, e.g., from
 *
 *	expr >= 0 complements newvar >= 0
 * to
 *	expr - newvar >= 0 complements newvar >= 0
 *
 * which makes newvar's Jacobian column nonzero but introduces
 * nonuniqueness (for expr > 0, we can have either newvar = 0
 * or newvar = expr).  Preliminary testing suggests that this
 * change (now the default) often saves iterations and takes
 * fewer pivots than the alternative of treating sinqle
 * inequalities the same as double inequalities, i.e.,
 * converting expr >= 0 to
 *
 *	newvar1 - expr	complements  -Infinity <= newvar2 <= Infinity
 *	newvar2		complements  0 <= newvar1 <= Infinity
 *
 * which also avoids adding zero columns to the Jacobian matrix
 * and contributes only the new nonuniqueness that newvar2 can
 * assume any nonnegative value.
 */

	ncb = FI.ineq - nranges;
	nvb = FI.nvb - FI.nvr;
	FI.nnz += ncb + nvb;
	n = FI.ineq + nvb + ncb;
	x = FI.sba = (int*)M1alloc(n*sizeof(int));
	y = FI.cb = FI.sbae = x + nvb + ncb;
	FI.cbe = y + FI.ineq;
	if (FI.ineq) {
		n = n_con;
		Cvar = cvar;
		for(i = 0; i < n; i++)
			switch(FI.ctype[i]) {
			case C_LB:
			case C_UB:
			case C_FREE:
				*x++ = Cvar[i];
				/* no break */
			case C_RANGE:
				*y++ = i;
			}
		}
	if (n = FI.nvb) {
		j = FI.n - (n + FI.nvr);
		for(i = 0; i < n; i++, j++)
			switch(FI.vtype[i]) {
			case C_RANGE:
				j++;
				break;
			case C_LB:
			case C_UB:
				*x++ = j;
				break;
			}
		}
	}

 static FILE *
Fopen(char *name, char *kind)
{
	FILE *f;
	if (!(f = fopen(name, "w"))) {
		fprintf(Stderr, "Cannot open %s file \"%s\".\n", kind, name);
		exit(2);
		}
	return f;
	}

 int
main(int argc, char **argv)
{
	ASL *asl;
	FILE *nl, *logfl = NULL, *statusfl = NULL;
	Information info;
	MCP *mcp;
	char buf[1024], *stub;
	int i, n, tc;
	typedef struct { char *msg; int pathcode, code; } Sol_info;
	Sol_info *si, *sie;
	static Sol_info solinfo[] = {
		{ "Solution found", MCP_Solved, 0 },
		{ "Not enough progress", MCP_NoProgress, 510 },
		{ "Major iteration limit", MCP_MajorIterationLimit, 400 },
		{ "Cumulative iteration limit", MCP_MinorIterationLimit, 401 },
		{ "Time limit", MCP_TimeLimit, 402 },
		{ "User interrupt", MCP_UserInterrupt, 403 },
		{ "Inconsistent bounds", MCP_BoundError, 200 },
		{ "Domain error", MCP_DomainError, 503 },
		{ "Infeasible", MCP_Infeasible, 201 },
		{ "Solver error", MCP_Error, 501 }
		};

	Not_Used(argc);

	asl = ASL_alloc(ASL_read_fg);

	/* Get version of PATH */
	Oinfo.bsname = stub = Path_Version();
	Oinfo.version = (char*)M1alloc(strlen(stub) + 8);
	sprintf(Oinfo.version, "AMPL/%s", stub);

	asl->i.congrd_mode = 2;	/* sparse Jacobians */
	want_xpi0 = 1;	/* want X0 */
	stub = getstub(&argv, &Oinfo);
	if (wantfuncs)
		show_funcs();
	if (quitnow | wantfuncs)
		exit(0);
	if (!stub)
		usage_ASL(&Oinfo, 1);

	nl = jac0dim(stub, (fint)strlen(stub));

	opt = Path_GetOptions();
	Options_Default(opt);
	Options_SetDouble(opt, "infinity", Infinity);

	if (getopts(argv, &Oinfo))
		exit(2);
	Options_Display(opt);

	FI.nvsq = FI.ncsq = 0;
	if ((n = n_eqn + n_cc) != n_var) {
		if (sqwarn & 3) {
			sprintf(buf, "%s %s %d x %d.\n", Oinfo.bsname,
				"requires a square system; this one is",
				n, n_var);
			if (sqwarn & 2) {
				if (sqwarn & 1) {
					fprintf(Stderr, "%s", buf);
					exit(1);
					}
				solve_result_num = 504;
				goto quit_now;
				}
			}
		if ((n -= n_var) > 0)
			FI.nvsq = n;
		else
			FI.ncsq = -n;
		}
	if (statusfile)
		Output_SetStatus(statusfl = Fopen(statusfile, "statusfile"));
	if (logfile)
		Output_SetLog(logfl = Fopen(logfile, "logfile"));
	fg_read(nl, ASL_no_linear_cc_rhs_adjust);
	FI.n = mcp_adj(asl);

	FI.ineq = n_con - n_eqn - n_cc;
	if (sideineq & 3 && (n = FI.ineq + FI.nvb) > 0) {
		i = sideineq & 2 ? sprintf(buf, "%s: ", Oinfo.bsname) : 0;
		i += sprintf(buf+i,"%d%s%d%s%d%s%d%s%d%s",
			n + FI.nvr + nranges, " side inequalities:\n\t",
			FI.nvb - FI.nvr, " simple bounds,\n\t",
			FI.nvr, " range bounds,\n\t", FI.ineq - nranges,
			" single inequality constraints, and\n\t",
			nranges, " range constraints.");
		if (sideineq & 2) {
			if (sideineq & 1) {
				fprintf(Stderr, "%s\n", buf);
				exit(1);
				}
			solve_result_num = 505;
 quit_now:
			write_sol(buf, 0, 0, &Oinfo);
			goto bailout;
			}
		printf("%s\n", buf);
		need_nl = 0;
		}
	if ((FI.ineq + FI.nvb - FI.nvr - nranges) > 0
	 && (FI.jadj = !(sideineq & 4)))
		prepare_jadj(asl);

	info.generate_output = Output_Log | Output_Status | Output_Listing;
	info.use_start = True;
	info.use_basics = True;

	mcp = MCP_Create(FI.n, FI.nnz);
	MCP_Jacobian_Structure_Constant(mcp, 1);
	MCP_Jacobian_Data_Contiguous(mcp, 1);
	install_interface(mcp);

	tc = Path_Solve(mcp, &info);

	i = sprintf(buf, "%s: ", Oinfo.bsname);
	sie = solinfo + sizeof(solinfo) / sizeof(Sol_info);
	for(si = solinfo;; si++) {
		if (si >= sie) {
			i += sprintf(buf + i, "Unexpected return code %d", tc);
			solve_result_num = 500;
			break;
			}
		if (si->pathcode == tc) {
			i += sprintf(buf + i, si->msg);
			solve_result_num = si->code;
			break;
			}
		}
	i += sprintf(buf + i, ".\n%d iterations (%d for crash); %d pivots",
		info.major_iterations + info.crash_iterations,
		info.crash_iterations, info.minor_iterations);
	i += sprintf(buf + i, ".\n%d function, %d gradient evaluations.",
		info.function_evaluations, info.jacobian_evaluations);
	write_sol(buf, DenseVector_Array(MCP_GetModX(mcp)), 0, &Oinfo);
	MCP_Destroy(mcp);
 bailout:
	if (logfl != NULL)
		fclose(logfl);
	if (statusfl != NULL)
		fclose(statusfl);
	return 0;
	}

 static void
bad_n(char *who, char *what, int n)
{
	fprintf(Stderr,
	    "Bad dimension argument '%s' to %s: got %d, expected %d\nexiting . . .\n",
	    what, who, n, FI.n);
	exit(9);
	}

 static void
bad_nnz(char *who, char *what, int n)
{
	fprintf(Stderr,
	    "Bad dimension argument '%s' to %s: got %d, expected %d\nexiting . . .\n",
	    what, who, n, FI.nnz);
	exit(9);
	}

 static void
zshow(char *what, real *z, int n)
{
	int i;

	printf("\n%s:\n", what);
	for (i = 0; i < n; i++)
		printf("%d\t%.g\n", i, z[i]);
	}

 int
function_evaluation(void *v, int n, real *z, real *f)
{
	ASL *asl = cur_ASL;
	int *Cvar, i, j, nc, nr, *sb, *sbe;
	real *c, *lu;
	Not_Used(v);

	if (n != FI.n)
		bad_n("F", "n", n);
	nc = n_con + FI.ncsq;
	lu = LUrhs;
	conval(z, c = FI.c_scratch, 0);
	Cvar = cvar;
	for (i = 0; i < nc; i++, lu += 2) {
		j = Cvar[i];
		switch (FI.ctype[i]) {
		case C_COMPL:
			f[j] = c[i];
			break;
		case C_LB:
		case C_EQN:
			f[j] = c[i] - lu[0];
			break;
		case C_UB:
			f[j] = c[i] - lu[1];
			break;
		case C_RANGE:
			f[j] = c[i] - z[j+1];
			f[j+1] = z[j];
			break;
		case C_FREE:
		case C_ZERO:
			f[j] = 0;
			}
		}
	nr = n - (FI.nvb + FI.nvr);
	for (i = 0; i < FI.nvb; i++, nr++) {
		j = FI.vb[i];
		lu = LUv + 2*j;
		switch(FI.vtype[i]) {
		  case C_LB:
			f[nr] = z[j] - lu[0];
			break;
		  case C_UB:
			f[nr] = z[j] - lu[1];
			break;
		  case C_RANGE:
			f[nr] = z[j] - z[nr+1];
			++nr;
			f[nr] = z[nr-1];
			}
		}
	if (sb = FI.sba)
		for(sbe = FI.sbae; sb < sbe; ) {
			j = *sb++;
			f[j] -= z[j];
			}
	if (debug & 2) {
		printf("\nf = F(z):\ni\tz[i]\t\tf[i]\n");
		for (i = 0; i < n; i++)
			printf("%d\t%.g\t%.g\n", i, z[i], f[i]);
		}
	return 0;
	}

 void
sparse_struct(int n, int maxnnz, int *row, int *col, int *len)
{
	ASL *asl = cur_ASL;
	int *Cvar, i, j, ja, k, *ka, nc, nv, *vb, *vbe;
	cgrad *cg;
	char *vt;

	if (n != FI.n)
		bad_n("sparse_struct", "n", n);
	if (maxnnz < FI.nnz)
		bad_nnz("sparse_struct", "nnz", maxnnz);
	nv = n_var;
	ka = A_colstarts;
	for (i = 0; i < nv; i++) {
		col[i] = ka[i] + 1;
		len[i] = ka[i+1] - ka[i];
		}
	j = ka[nv] + 1;
	for(nv += FI.nvsq; i < nv; i++) {
		col[i] = j;
		len[i] = 0;
		}
	nc = n_con;
	Cvar = cvar;
	for (i = 0; i < nc; ) {
		j = Cvar[i] + 1;
		for (cg = Cgrad[i++]; cg; cg = cg->next)
			row[cg->goff] = j;
		}
	j = nzc + FI.nvb;
	nc += FI.ncsq;
	ja = FI.jadj;
	for (i = 0; i < nc; i++)
		switch (FI.ctype[i]) {
		case C_LB:
		case C_UB:
		case C_FREE:
			col[nv] = j+1;
			if (len[nv++] = ja)
				row[j++] = Cvar[i] + 1;
			break;

		case C_RANGE:
			row[j++] = k = Cvar[i] + 2;
			col[nv] = j;
			len[nv++] = 1;

			row[j++] = k - 1;
			col[nv] = j;
			len[nv++] = 1;
		}
	if (FI.nvb) {
		k = nc + nranges;
		vb = FI.vb;
		vbe = vb + FI.nvb;
		vt = FI.vtype;
		++ka;
		do {
			i = *vb++;
			if (*vt++ == C_RANGE) {
				row[j++] = k += 2;
				col[nv] = j;
				len[nv++] = 1;

				row[ka[i]-1] = row[j++] = k-1;
				col[nv] = j;
				len[nv++] = 1;
				}
			else {
				col[nv] = j+1;
				row[ka[i]-1] = ++k;
				if (len[nv++] = ja)
					row[j++] = k;
				}
			} while (vb < vbe);
		}
	if (debug & 8) {
		printf("\nJacobian nonzeros:\ni\tcol[i]\trow[i]\n");
		for(i = j = 0; j < n; j++)
			for(k = i + len[j]; i < k; i++)
				printf("%d\t%d\t%d\n", i, j, row[i]);
		}
	}

 int
sparseJ(int n, int maxnnz, double *z, double *J)
{
	ASL *asl = cur_ASL;
	int *cb, *cbe, i, j, *ka;

	if (n != FI.n)
		bad_n("sparseJ", "n", n);
	if (maxnnz < FI.nnz)
		bad_nnz("sparseJ", "nnz", maxnnz);
	if (debug & 4)
		zshow("z in sparseJ", z, n);
	jacval(z, J, 0);
	j = nzc + FI.nvb;
	if (FI.sba) {
		cb = FI.cb;
		cbe = FI.cbe;
		while(cb < cbe)
			switch(FI.ctype[*cb++]) {
			case C_LB:
			case C_UB:
			case C_FREE:
				J[j++] = -1;
				break;
			case C_RANGE:
				J[j++] = 1;
				J[j++] = -1;
			}
		for(i = 0; i < FI.nvb; i++)
			switch(FI.vtype[i]) {
			case C_LB:
			case C_UB:
			case C_FREE:
				J[j++] = -1;
				break;
			case C_RANGE:
				J[j++] = 1;
				J[j++] = -1;
			}
		}
	else
		for(i = FI.nvr + nranges; i > 0; --i) {
			J[j++] = 1;
			J[j++] = -1;
			}
	ka = A_colstarts + 1;
	--J;
	for(i = 0; i < FI.nvb; i++)
		J[ka[FI.vb[i]]] = 1;
	if (debug & 4)
		zshow("J in sparseJ", J+1, FI.nnz);
	return 0;
	}

 void
problem_size(void *v, int *n, int *nnz)
{
#ifdef Student_Edition
	ASL *asl = cur_ASL;
	Not_Used(v);

	if (FI.n  > CRIPPLE_SIZE || FI.nnz > CRIPPLE_NNZ) {
		fflush(stdout);
		fprintf(Stderr,
"\nSorry, the student edition is limited to\n\
%d variables and %d Jacobian nonzeros.\n\
%sou have %d variables and %d nonzeros.\n",
			CRIPPLE_SIZE, CRIPPLE_NNZ,
			n_var != FI.n || nzc != FI.nnz
			? "After problem adjustments,\ny" : "Y",
			FI.n, FI.nnz);
		exit(1);
		}
#else
	Not_Used(v);
#endif
	*n = FI.n;
	*nnz = FI.nnz;
	}

 void
bounds(void *v, int n, real *z, real *lower, real *upper)
{
	ASL *asl = cur_ASL;
	int i, j, k, nc, nv;
	real Lb, Ub, *lu;
	Not_Used(v);

	if (n != FI.n)
		bad_n("bounds", "n", n);
	nv = n_var;
	lu = LUv;
	if (X0)
		memcpy(z, X0, nv * sizeof(real));
	else
		memset(z, 0, nv * sizeof(real));
	for (i = 0; i < nv; i++, lu += 2) {
		lower[i] = lu[0];
		upper[i] = lu[1];
		}
	k = nc = n_con + FI.ncsq;
	for (i = 0; i < nc; i++)
		if ((j = cvar[i]) >= nv) {
			switch (FI.ctype[i]) {
			case C_LB:
				Lb = 0;
				Ub = Infinity;
				break;
			case C_UB:
				Lb = negInfinity;
				Ub = 0;
				break;
			case C_RANGE:
				lu = LUrhs + (i << 1);
				k = j + 1;
				z[k] = lu[0] + 0.5 * (lu[1] - lu[0]);
				lower[k] = lu[0];
				upper[k] = lu[1];
				/* no break; */
			case C_EQN:
			case C_ZERO:
				Lb = negInfinity;
				Ub = Infinity;
				break;
			case C_FREE:
				Lb = 0;
				Ub = 0;
				break;
			default:
				fprintf(Stderr,
				    "bounds: FI.ctype[%d] = %d\n",
				    i, FI.ctype[i]);
				exit(1);
				}
			z[j] = 0;
			lower[j] = Lb;
			upper[j] = Ub;
			}
	for(i = 0, k = n - FI.nvb - FI.nvr; k < n; i++, k++) {
		z[k] = 0;
		switch(FI.vtype[i]) {
		  case C_LB:
			lower[k] = 0;
			upper[k] = Infinity;
			break;
		  case C_UB:
			lower[k] = negInfinity;
			upper[k] = 0;
			break;
		  case C_RANGE:
			lower[k] = negInfinity;
			upper[k] = Infinity;
			z[++k] = 0;
			j = FI.vb[i];
			lu = LUv + 2*j;
			lower[k] = lu[0];
			upper[k] = lu[1];
			break;
		  default:
			fprintf(Stderr, "bounds: FI.vtype[%d] = %d\n",
				i, FI.vtype[i]);
			exit(1);
		  }
		}
	for(i = FI.nvb; i-- > 0; ) {
		j = FI.vb[i];
		lower[j] = negInfinity;
		upper[j] = Infinity;
		}
	if (debug & 1) {
		printf("i\tz[i]\tlower[i]\tupper[i]\n");
		for(i = 0; i < n; i++)
			printf("%d\t%.g\t%.g\t%.g\n", i, z[i],
				lower[i], upper[i]);
		}
	return;
	}

 void
variable_name(int variable, char *buffer, int buffer_size)
{
	ASL * asl = cur_ASL;

	strncpy(buffer, var_name(variable - 1), buffer_size - 1);
	buffer[buffer_size-1] = '\0';
	return;
	}

 void
constraint_name(int constraint, char *buffer, int buffer_size)
{
	ASL * asl = cur_ASL;

	strncpy(buffer, con_name(constraint - 1), buffer_size - 1);
	buffer[buffer_size-1] = '\0';
	return;
	}

 int
jacobian_evaluation(void *v, int n, double *z, int  wantf, double *f,
	int *nnz, int *col_start, int *col_len, int *row, double *data)
{
	int	err;

	if (wantf)
		err = function_evaluation(v, n, z, f);
	else
		err = 0;

	if (numJacEvals == 0)
		sparse_struct(n, *nnz, row, col_start, col_len);
	numJacEvals++;
	err += sparseJ(n, *nnz, z, data);

	*nnz = col_start[n - 1] + col_len[n - 1] - 1;
	return err;
	}

 static void
jac_typ(void *v, int nnz, int *typ)
{
	ASL *asl = cur_ASL;
	cgrad *cg;
	int i, n;
	Not_Used(v);

	for (i = 0, n = nnz; i < n; i++)
		typ[i] = PRESOLVE_LINEAR;

	for (i = 0, n = nlc; i < n; i++)
		for (cg = Cgrad[i]; cg; cg = cg->next)
			if (cg->varno < nlvc)
				typ[cg->goff] = PRESOLVE_NONLINEAR;
	}

 static MCP_Interface
interface = {
	NULL,
	problem_size, bounds,
	function_evaluation, jacobian_evaluation,
	NULL, NULL,
	NULL, NULL,
	NULL
	};

 static Presolve_Interface
presolve_interface = {
	NULL,
	NULL, NULL,
	NULL, NULL,
	jac_typ
	};

 static void
install_interface(MCP *m)
{
	MCP_SetInterface(m, &interface);
	MCP_SetPresolveInterface(m, &presolve_interface);
	}
