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

/* This driver works with lp_solve 3.x and 4.x. */

#include "lpkit.h"
#include "patchlevel.h"
#define STDIO_H_included
#include "getstub.h"

 static int nobj = 1, prlp, scaling;

#ifdef __PATCHLEVEL_H__ /*lp_solve >= 4.0*/
 static int scalemode;

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
scalemode_msg[] = "=... scale mode:  one of\n\
			0 for default scaling,\n\
			1 for geometric scaling, or\n\
			4 for Curtis & Reid scaling\n\
		plus\n\
			2 for powerscaling (not with 4)\n\
			8 to scale \"Lagrange\" variables\n\
			16 to scale integer variables",
bb_msg[] = "=... branch-and-bound rule: one of\n\
			0 (default) for lowest indexted non-integer variable\n\
			1 for a random non-integer variable\n\
			2 for largest deviation from integer\n\
			4 for median deviation from integer";
#else
#define MYBOOL unsigned short
#define set_scalemode(x,y) /*nothing*/
#endif /*__PATCHLEVEL_H__*/

 static char *
SK0_val(Option_Info *oi, keyword *kw, char *v)
{
	*(MYBOOL*)(oi->uinfo + Intcast kw->info) = 0;
	return v;
	}

 static char *
SK1_val(Option_Info *oi, keyword *kw, char *v)
{
	*(MYBOOL*)(oi->uinfo + Intcast kw->info) = 1;
	return v;
	}

#define v_off(x) voffset_of(lprec,x)

 static keyword
keywds[] = {
#ifdef __PATCHLEVEL_H__
 KW("bb", SK_val, v_off(bb_rule), bb_msg),
#endif
 KW("cfirst", SK0_val, v_off(floor_first), "in IPs, take ceiling branch first"),
 KW("debug", SK1_val, v_off(debug), "debug mode"),
 KW("degen", IK1_val, v_off(anti_degen), "perturb degeneracies"),
 KW("eps",   DU_val, v_off(epsilon), "=... tolerance for rounding to integer"),
#ifdef __PATCHLEVEL_H__
 KW("improve", SK_val, v_off(improve), "=... what to: sum of 1 for FTRAN, 2 for BTRAN; default = 0"),
#endif
 KW("objno", I_val, &nobj, "=... objective number: 0 = none, 1 (default) = first"),
 KW("obound", DU_val, v_off(obj_bound), "=... initial objective bound (for branch & bound)"),
 KW("prinv", SK1_val, v_off(print_at_invert), "print after inverting"),
 KW("prlp",  IK1_val, &prlp, "print the LP"),
 KW("psols", SK1_val, v_off(print_sol), "print (intermediate) feasible solutions"),
 KW("scale", IK1_val, &scaling, "scale the problem"),
#ifdef __PATCHLEVEL_H__
 KW("scalemode", I_val, &scalemode, scalemode_msg),
#endif
 KW("trace", SK1_val, v_off(trace), "trace pivot selections"),
 KW("verbose", SK1_val, v_off(verbose), "verbose mode"),
 KW("version", Ver_val, 0, "report version details"),
 KW("wantsol", WS_val, 0, WS_desc_ASL)
 };

 static char lp_solve_Version[48] = "\nAMPL/LP_SOLVE Driver Version 20030612\n";
 static Option_Info Oinfo = { "lpsolve", lp_solve_Version+6, "lpsolve_options",
				keywds, nkeywds, 0, lp_solve_Version+1,
				0,0,0,0,0, 20030612 };

 static SufDecl
suftab[] = {
	{ "sstatus", 0, ASL_Sufkind_var | ASL_Sufkind_outonly },
	{ "sstatus", 0, ASL_Sufkind_con | ASL_Sufkind_outonly }
	};

 int
main(int argc, char **argv)
{
	char *stub;
	ASL *asl;
	FILE *nl;
	lprec *lp;
	ograd *og;
	int ct, i, intmin, *is, j, j0, j1, k, nalt, rc;
	MYBOOL *basis, *lower;
	real *LU, *c, lb, objadj, *rshift, *shift, t, ub, *x, *x0, *x1;
	char buf[256];
	typedef struct { char *msg; int code; } Sol_info;
	void *SI;
	static Sol_info solinfo[] = {
		{ "optimal", 0 },
		{ "integer programming failure", 502 },
		{ "infeasible", 200 },
		{ "unbounded", 300 },
		{ "failure", 501 },
		{ "bug", 500 }
		};

#ifdef PATCHLEVEL
	snprintf(lp_solve_Version+15, sizeof(lp_solve_Version)-15,
		"%.*s", (int)sizeof(lp_solve_Version)-16, PATCHLEVEL);
#else
	lp_solve_version(&i, &j, &j0, &j1);
	snprintf(lp_solve_Version+15, sizeof(lp_solve_Version)-15,
		"%d.%d.%d.%d", i, j, j0, j1);
#endif
	asl = ASL_alloc(ASL_read_f);
	stub = getstub(&argv, &Oinfo);
	if (!stub)
		usage_ASL(&Oinfo, 1);
	nl = jac0dim(stub, 0);
	SI = sos_add(nl,0);
	suf_declare(suftab, sizeof(suftab)/sizeof(SufDecl));

	/* set A_vals to get the constraints column-wise */
	A_vals = (real *)M1alloc(nzc*sizeof(real));

	f_read(nl,0);

	sos_finish(&SI, 0, &j, 0, 0, 0, 0, 0);

	lp = make_lp(n_con, 0);

	Oinfo.uinfo = (char *)lp;
	if (getopts(argv, &Oinfo))
		return 1;

	i = n_var + n_con + 1;
	x = (real*)M1alloc(i*sizeof(real));	/* scratch vector */
	memset(x, 0, i*sizeof(real));
	x0 = x++;
	c = x + n_con;

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
	intmin = n_var - (nbv + niv);
	j1 = nalt = 0;
	rshift = shift = 0;
	for(i = 1; i <= n_var; i++, LU += 2) {
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
						  (n_var+n_con)*sizeof(real));
						shift = rshift + n_con - 1;
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
	for(i = 1; i <= n_con; i++, LU += 2) {
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
	if (scaling) {
		set_scalemode(lp, scalemode);
		auto_scale(lp);
		}

	/* Unfortunately, there seems to be no way to suggest */
	/* a starting basis to lp_solve; thus we must ignore  */
	/* any incoming .sstatus values. */

	rc = solve(lp);
	if (rc < 0 || rc > 5)
		rc = 5;
	solve_result_num = solinfo[rc].code;
	i = sprintf(buf, "%s: %s", Oinfo.bsname, solinfo[rc].msg);
	if (rc == OPTIMAL)
		i += sprintf(buf+i, ", objective %.*g", obj_prec(),
			lp->best_solution[0]);
	i += sprintf(buf+i,"\n%d simplex iterations", lp->total_iter);
	if (lp->max_level > 1 || lp->total_nodes > 1)
		sprintf(buf+i, "\n%d branch & bound nodes: depth %d",
			lp->total_nodes, lp->max_level);

	/* Prepare to report solution: deal with split free variables. */

	x1 = lp->best_solution+lp->rows+1;
	if (nalt || shift) {
		x = x0;
		LU = LUv;
		for(i = 0; i < n_var; i++, LU += 2) {
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

	if (solinfo[rc].code < 500 && !(nbv + niv)) {

		/* return .sstatus values */

		basis = lp->basis;
		lower = lp->lower;
		is = M1alloc((n_var + n_con)*sizeof(int));
		suf_iput("sstatus", ASL_Sufkind_con, is);
		for(i = 0; i < n_con; i++) {
			j = *++lower;
			*is++ = *++basis ? 1 : j ? 3 : 4;
			}
		suf_iput("sstatus", ASL_Sufkind_var, is);
		LU = LUv;
		for(i = 0; i < n_var; i++, LU += 2) {
			j0 = *++basis;
			j1 = *++lower;
			if (LU[0] > negInfinity)
				j = j0 ? 1 : j1 ? 3 : 4;
			else if (LU[1] < Infinity)
				j = j0 ? 1 : j1 ? 4 : 3;
			else {
				++lower;
				j = *++basis || j0;
				}
			*is++ = j;
			}
		}

	write_sol(buf, x, lp->duals+1, &Oinfo);
	/* The following calls would only be needed */
	/* if execution were to continue... */
	delete_lp(lp);
	ASL_free(&asl);
	return 0;
	}
