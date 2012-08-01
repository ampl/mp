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

/* This driver works with lp_solve 5.x (at least for x >= 5). */

#include "lpkit.h"
#define STDIO_H_included
#include "getstub.h"

 static int bb_first = -1, debug, degen, nobj = 1, pivoptions, pivoting = -1,
		presolve, prlp, psols, scalemode, simplextype, trace,
		verbose;
 static lprec *lpx;
 static I_Known cfirst	= { 0, &bb_first },
		ffirst	= { 1, &bb_first },
		piva	= { PRICE_ADAPTIVE, &pivoptions },
		pivf	= { PRICE_PRIMALFALLBACK, &pivoptions },
		pivla	= { PRICE_LOOPALTERNATE, &pivoptions },
		pivll	= { PRICE_LOOPLEFT, &pivoptions },
		pivm	= { PRICE_MULTIPLE, &pivoptions },
		pivp	= { PRICE_AUTOPARTIAL, &pivoptions },
		pivr	= { PRICE_RANDOMIZE, &pivoptions },
		pres0	= { PRESOLVE_ROWS | PRESOLVE_COLS, &presolve },
		presl	= { PRESOLVE_LINDEP, &presolve },
		presr	= { PRESOLVE_REDUCEMIP, &presolve },
		psolsa	= { 2, &psols },
		simpdd	= { SIMPLEX_DUAL_DUAL, &simplextype },
		simpdp	= { SIMPLEX_DUAL_PRIMAL, &simplextype },
		simppd	= { SIMPLEX_PRIMAL_DUAL, &simplextype },
		simppp	= { SIMPLEX_PRIMAL_PRIMAL, &simplextype };

 typedef struct
Isg_val {
	int (__WINAPI *Get)(lprec*);
	void (__WINAPI *Set)(lprec*, int);
	} Isg_val;

 typedef struct
Rsg_val {
	real (__WINAPI *Get)(lprec*);
	void (__WINAPI *Set)(lprec*, real);
	} Rsg_val;

 static Isg_val
		improve = { get_improve, set_improve };

 static Rsg_val
		epsel	= { get_epsel,	set_epsel },
		epsint	= { get_epsint,	set_epsint },
		epspiv	= { get_epspivot,set_epspivot },
		epsrc	= { get_epsd,	set_epsd },
		epsrhs	= { get_epsb,	set_epsb };
		/* obound = { get_objbound,set_objbound }; */

 static char *
Igs_val(Option_Info *oi, keyword *kw, char *v)
{
	Isg_val *iv;
	char *rv;
	int b, n;

	iv = (Isg_val*)kw->info;
	if (*v == '?') {
		n = iv->Get(lpx);
		printf("%s%s%d\n", kw->name, oi->eqsign, n);
		return v+1;
		}
	b = 10;
	if (*v == '0' && (v[1] == 'x' || v[1] == 'X')) {
		v += 2;
		b = 16;
		}
	n = (int)strtol(v, &rv, b);
	if (*(unsigned char *)rv <= ' ')
		iv->Set(lpx, n);
	else
		rv = badval_ASL(oi, kw, v, rv);
	return rv;
	}

 static char *
Rgs_val(Option_Info *oi, keyword *kw, char *v)
{
	Rsg_val *g;
	char buf[32], *rv;
	real t;

	g = (Rsg_val*)kw->info;
	if (*v == '?') {
		t = g->Get(lpx);
		g_fmt(buf, t);
		printf("%s%s%s\n", kw->name, oi->eqsign, buf);
		return v+1;
		}
	t = strtod(v, &rv);
	if (*(unsigned char *)rv <= ' ')
		g->Set(lpx, t);
	else
		rv = badval_ASL(oi, kw, v, rv);
	return rv;
	}

 static char *
SK_val(Option_Info *oi, keyword *kw, char *v)
{
	char *rv;
	int i;
	rv = Ival_ASL(oi, kw, v, &i);
	*(MYBOOL*)(oi->uinfo + Intcast kw->info) = i;
	return rv;
	}

 static char
antidegen_msg[] = "... anti-degeneracy scheme: sum of\n\
		 0 (default) nothing\n\
		 1 driving equality slacks from the basis in Phase 1\n\
		 2 ColumnCheck\n\
		 4 Stalling\n\
		 8 NumFailure\n\
		16 LostFeas\n\
		32 Infeasible\n\
		64 Dynamic\n\
	       128 Applying above during BB",
bb_msg[] = "=... branch-and-bound rule: based on one of\n\
			0 (default) lowest indexted non-integer variable\n\
			1 distance from the current bounds\n\
			2 largest current bound\n\
			3 largest fractional value\n\
			4 unweighted pseudo-costs of variables\n\
			5 minimizing number of integer infeasibilities\n\
			6 maximizing the normal pseudo-cost\n\
	    plus\n\
		   8 for weight reverse mode\n\
		  16 to reverse the direction choice when cauto is used\n\
		  32 for greedy mode\n\
		  64 for pseudo-cost mode\n\
		 128 to select the node most often previously selected\n\
		 256 for \"randomize\" mode\n\
		1024 to turn mode 128 off at the first feasible solution\n\
		2048 for restart mode\n\
		4096 to selected a node least-frequently selected so far\n",
improve_msg[] = "=... the iterative improvement level: one of\n\
			0 improve none.\n\
			1 Running accuracy measurement of solved equations based on\n\
				Bx=r (primal simplex), remedy is refactorization.\n\
			2 Improve initial dual feasibility by bound flips\n\
				(default, highly recommended)\n\
			4 Low-cost accuracy monitoring in the dual,\n\
				remedy is refactorization.\n\
			8 By default there is a check for primal/dual feasibility\n\
				at an optimum only for the relaxed problem; this\n\
				also activates the test at the node level.",
pivoting_msg[] = "=... simplex pivot rule: one of\n\
			0 first available pivot\n\
			1 Dantzig\n\
			2 Devex (default = Devix + piva)\n\
			3 steepest edge\n\
		possibly modified by one or more of the following piv* keywords,\n\
		or by adding to pivoting the number mentioned with the piv* keyword:",
pres0_msg[] = "use lp_solve's presolver (without presolvel or presolver\n\
			unless one or both is explicitly given",
presl_msg[] = "use lp_solve's presolver and eliminate linearly-dependent rows",
presr_msg[] = "use lp_solve's presolver and delete constraints found redundant\n\
			in phase I",
scalemode_msg[] = "=... scale mode:  one of\n\
			 1 for scale to convergence using largest absolute value\n\
			 2 for scale based on the simple numerical range\n\
			 3 for numerical range-based scaling\n\
			 4 for geometric scaling\n\
			 7 for Curtis & Reid scaling\n\
		plus\n\
			16 for scaling to logarithmic mean\n\
			32 to also do Power scaling\n\
			64 to make sure that no scaled number is above 1\n\
		       128 to scale integer variables\n\
		       256 to allow dynamic update\n\
		       512 to only scale constraints (rows)\n\
		      1024 to only scale variables (columns)";

#define v_off(x) voffset_of(lprec,x)

 static keyword
keywds[] = {
 KW("bb",	SK_val, v_off(bb_rule), bb_msg),
 KW("cfirst",	IK_val, &cfirst, "in IPs, take ceiling branch first"),
 KW("debug",	IK1_val, &debug, "debug mode"),
 KW("degen",	IK1_val, &degen, "same as degenx=1"),
 KW("degenx",	I_val, &degen, antidegen_msg),
 KW("eps",	Rgs_val, &epsint, "=... tolerance for rounding to integer"),
 KW("epsel",	Rgs_val, &epsel, "=... zero tolerance for matrix elements"),
 KW("epspiv",	Rgs_val, &epspiv, "=... minimum pivot"),
 KW("epsrc",	Rgs_val, &epsrc, "=... zero tolerance for reduced costs"),
 KW("epsrhs",	Rgs_val, &epsrhs, "=... zero tolerance for the RHS"),
 KW("ffirst",	IK_val, &ffirst, "in IPs, take floor branch first"),
 KW("improve",	Igs_val, &improve, improve_msg),
 KW("objno",	I_val, &nobj, "=... objective number: 0 = none, 1 (default) = first"),
/* KW("obound",	Rgs_val, &obound, "=... initial objective bound (for branch & bound)"),*/
 KW("piv", I_val, &pivoting, pivoting_msg),
 KW("piva", IK_val, &piva, "(32) temporarily use first index if cycling is detected"),
 KW("pivf", IK_val, &pivf, "(4) with steepest-edge, use Devex in primal"),
 KW("pivla", IK_val, &pivla, "(2048) alternate left and right scanning"),
 KW("pivll", IK_val, &pivll, "(1024) scan to the left rather than right"),
 KW("pivm", IK_val, &pivm, "(8) do multiple pricing"),
 KW("pivp", IK_val, &pivp, "(256) enable partial pricing"),
 KW("pivr", IK_val, &pivr, "(128) randomly perturb the selected pivot"),
 KW("presolve",  IK_val, &pres0, pres0_msg),
 KW("presolvel", IK_val, &presl, presl_msg),
 KW("presolver", IK_val, &presr, presr_msg),
 KW("prlp",  IK1_val, &prlp, "print the LP"),
 KW("psols", IK1_val, &psols, "print (intermediate) feasible solutions"),
 KW("psolsa", IK_val, &psolsa, "print (intermediate) feasible solutions (non-zeros)"),
 KW("scale", I_val, &scalemode, scalemode_msg),
 KW("simplexdd", IK_val, &simpdd, "dual Phase I and II"),
 KW("simplexdp", IK_val, &simpdp, "dual Phase I, primal Phase II (default)"),
 KW("simplexpd", IK_val, &simppd, "primal Phase I, dual Phase II"),
 KW("simplexpp", IK_val, &simppp, "primal Phase I and II"),
 KW("trace", IK1_val, &trace, "trace pivot selections"),
 KW("verbose", IK1_val, &verbose, "verbose mode"),
 KW("version", Ver_val, 0, "report version details"),
 KW("wantsol", WS_val, 0, WS_desc_ASL)
 };

#ifndef YYYYMMDD
#include "lpsolve5_date.h"
#endif
 static char lp_solve_Version[48] = "\nAMPL/LP_SOLVE Driver Version " qYYYYMMDD "\n";
 static Option_Info Oinfo = { "lpsolve", lp_solve_Version+6, "lpsolve_options",
				keywds, nkeywds, 0, lp_solve_Version+1,
				0,0,0,0,0, YYYYMMDD };

 static SufDecl
suftab[] = {
	{ "sstatus", 0, ASL_Sufkind_var | ASL_Sufkind_outonly },
	{ "sstatus", 0, ASL_Sufkind_con | ASL_Sufkind_outonly }
	};

 int
main(int argc, char **argv)
{
	ASL *asl;
	FILE *nl;
	char buf[256], *stub;
	int *bi, *cs, ct, i, ib, intmin, j, j0, j1, k, m, m1, maxlevel, n, n1, nalt, rc, *vs;
	lprec *lp;
	ograd *og;
	real *LU, *c, lb, objadj, *rshift, *shift, t, ub, *x, *x0, *x1, *y;
	unsigned long totalnodes;
	void *SI;
	typedef struct { char *msg; int code; } Sol_info;
	static Sol_info solinfo[] = {
		{ "optimal", 0 },
		{ "integer programming failure", 502 },
		{ "infeasible", 200 },
		{ "unbounded", 300 },
		{ "failure", 501 },
		{ "bug", 500 }
		};

	lp_solve_version(&i, &j, &j0, &j1);
	snprintf(lp_solve_Version+15, sizeof(lp_solve_Version)-15,
		"%d.%d.%d.%d", i, j, j0, j1);
	asl = ASL_alloc(ASL_read_f);
	stub = getstub(&argv, &Oinfo);
	if (!stub)
		usage_ASL(&Oinfo, 1);
	nl = jac0dim(stub, 0);
	SI = sos_add(nl,0);
	suf_declare(suftab, sizeof(suftab)/sizeof(SufDecl));
	n = n_var;
	m = n_con;

	/* set A_vals to get the constraints column-wise */
	A_vals = (real *)M1alloc(nzc*sizeof(real));

	f_read(nl,0);

	sos_finish(&SI, 0, &j, 0, 0, 0, 0, 0);

	lpx = lp = make_lp(m, 0);

	Oinfo.uinfo = (char *)lp;
	if (getopts(argv, &Oinfo))
		return 1;

	if (bb_first != -1)
		set_bb_floorfirst(lp, bb_first);
	if (debug)
		set_debug(lp, 1);
	if (degen)
		set_anti_degen(lp, degen);
	if (pivoting >= 0)
		set_pivoting(lp, pivoting | pivoptions);
	if (presolve)
		set_presolve(lp, presolve, get_presolveloops(lp));
	if (psols)
		set_print_sol(lp, psols);
	if (simplextype)
		set_simplextype(lp, simplextype);
	if (trace)
		set_trace(lp, 1);
	set_verbose(lp, verbose ? NORMAL : SEVERE);

	i = n + m + 1;
	x = (real*)M1alloc(i*sizeof(real));	/* scratch vector */
	memset(x, 0, i*sizeof(real));
	x0 = x++;
	c = x + m;

	/* supply objective */

	objadj = 0;
	if (--nobj >= 0 && nobj < n_obj) {
		for(og = Ograd[nobj]; og; og = og->next)
			c[og->varno] = og->coef;
		if (objtype[nobj])
			set_maxim(lp);
		objadj = objconst(nobj);
		}

	/* supply columns and variable bounds */

	LU = LUv;
	intmin = n - (nbv + niv);
	j1 = nalt = 0;
	rshift = shift = 0;
	for(i = 1; i <= n; i++, LU += 2) {
		lb = LU[0];
		ub = LU[1];
		j0 = j1;
		j1 = A_colstarts[i];
		*x0 = *c++;	/* cost coefficient */
		if (lb <= negInfinity && ub < Infinity) {
			/* negate this variable */
			nalt++;
			lb = -ub;
			ub = -LU[0];
			for(j = j0; j < j1; j++)
				x[A_rownos[j]] = -A_vals[j];
			*x0 = -*x0;
			add_column(lp, x0);
			if (lb)
				goto shift_check;
			}
		else {
			for(j = j0; j < j1; j++)
				x[A_rownos[j]] = A_vals[j];
			add_column(lp, x0);
			if (lb <= negInfinity) {
				nalt++;
				if (i > intmin)
					set_int(lp, lp->columns, TRUE);
				/* split free variable */
				*x0 = -*x0;
				for(j = j0; j < j1; j++)
					x[A_rownos[j]] *= -1.;
				add_column(lp,x0);
				}
			else if (lb) {
 shift_check:
				if (lb > 0)
					set_lowbo(lp, lp->columns, lb);
				else {
					if (!rshift) {
						rshift = (real*)M1zapalloc(
						  (m+n)*sizeof(real));
						shift = rshift + m - 1;
						}
					shift[i] = lb;
					for(j = j0; j < j1; j++) {
						k = A_rownos[j];
						rshift[k] += lb*x[k];
						}
					if (ub < Infinity)
						ub -= lb;
					objadj += lb**x0;
					}
				}
			if (ub < Infinity)
				set_upbo(lp, lp->columns, ub);
			}
		for(j = j0; j < j1; j++)
			x[A_rownos[j]] = 0;
		if (i > intmin)
			set_int(lp, lp->columns, TRUE);
		}

	if (objadj) {
		/* add a fixed variable to adjust the objective value */
		*x0 = objadj;
		add_column(lp, x0);
		set_lowbo(lp, i, 1.);
		set_upbo(lp, i, 1.);
		}

	/* supply constraint rhs */

	LU = LUrhs;
	for(i = 1; i <= m; i++, LU += 2) {
		t = LU[0];
		if (t == LU[1])
			ct = EQ;
		else if (t <= negInfinity) {
			t = LU[1];
			if (t >= Infinity) {
				/* This is possible only with effort: */
				/* one must turn presolve off and */
				/* explicitly specify a constraint */
				/* with infinite bounds. */
				fprintf(Stderr,
					"Sorry, can't handle free rows.\n");
				exit(1);
				}
			ct = LE;
			}
		else
			ct = GE;
		set_constr_type(lp, i, ct);
		set_rh(lp, i, rshift ? t - *rshift++ : t);
		if (ct == GE && LU[1] < Infinity)
			lp->orig_upbo[i] = LU[1] - t;
		}

	if (prlp)
		print_lp(lp);
	if (scalemode)
		set_scaling(lp, scalemode);

	/* Unfortunately, there seems to be no way to suggest */
	/* a starting basis to lp_solve; thus we must ignore  */
	/* any incoming .sstatus values. */

	rc = solve(lp);
	if (rc < 0 || rc > 5)
		rc = 5;
	solve_result_num = solinfo[rc].code;
	ib = sprintf(buf, "%s: %s", Oinfo.bsname, solinfo[rc].msg);
	if (rc == OPTIMAL)
		ib += sprintf(buf+ib, ", objective %.*g", obj_prec(),
			lp->best_solution[0]);
	ib += sprintf(buf+ib,"\n%d simplex iterations", lp->total_iter);
	maxlevel = get_max_level(lp);
	totalnodes = get_total_nodes(lp);
	if (maxlevel > 1 || totalnodes > 1)
		ib += sprintf(buf+ib, "\n%lu branch & bound nodes: depth %d",
			totalnodes, maxlevel);

	/* Prepare to report solution: deal with split free variables. */

	if (!get_ptr_variables(lp,&x1))
		x = 0;
	else if (nalt || shift) {
		x = x0;
		LU = LUv;
		for(i = 0; i < n; i++, LU += 2) {
			if (LU[0] > negInfinity)
				x[i] = *x1++;
			else if (LU[1] < Infinity)
				x[i] = -*x1++;
			else {
				x[i] = x1[0] - x1[1];
				x1 += 2;
				}
			if (shift)
				x[i] += *++shift;
			}
		}
	else
		x = x1;
	if (get_ptr_dual_solution(lp, &y))
		++y;
	else
		y = 0;

	if (solinfo[rc].code < 500 && !(nbv + niv)) {

		/* return .sstatus values */

		n1 = get_Ncolumns(lp);
		m1 = get_Nrows(lp);
		k = m1 + n1;
		cs = (int*)M1alloc((m + n + k + 1)*sizeof(int));
		vs = cs + m;
		bi = vs + n;
		if (!get_basis(lp, bi, 1)) {
			ib += sprintf(buf+ib, "\nNo basis returned.");
			goto no_basis;
			}
		memset(cs, 0, (m+n)*sizeof(int));
		for(i = 1; i <= m1; ++i) {
			if ((j = bi[i]) < 0)
				j = -j;
			if (--j < n)
				vs[j] = 1;
			else
				cs[j-n] = 1;
			}
		for(; i <= k; ++i) {
			k = 4;
			if ((j = bi[i]) < 0) {
				j = -j;
				k = 3;
				}
			if (--j < n)
				vs[j] = k;
			else
				cs[j-n] = k;
			}
		suf_iput("sstatus", ASL_Sufkind_con, cs);
		suf_iput("sstatus", ASL_Sufkind_var, vs);
		}

 no_basis:
	write_sol(buf, x, y, &Oinfo);
	/* The following calls would only be needed */
	/* if execution were to continue... */
	delete_lp(lp);
	ASL_free(&asl);
	return 0;
	}
