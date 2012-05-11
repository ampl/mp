/****************************************************************
Copyright (C) 1997-9 Lucent Technologies
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

static int ma, nobj = 1;
fint minmax, outlev;
static real bigb = 1.e15, mf = 6.;

Cextern int readpar_ ANSI((char *nam1, fint *outlev, fint *minmax,
	fint *maxmn, char *namstr, char *objnam, char *rhsnam,
	char *bndnam, char *rngnam, real *bigb, fint *inslck,
	ftnlen nam1_len, ftnlen namstr_len, ftnlen objnam_len,
	ftnlen rhsnam_len, ftnlen bndnam_len, ftnlen rngnam_len));

Cextern int solver_ ANSI((real *obj, real *rhs, real *lbound,
	real *ubound, real *diag, real *odiag,
	real *xs, real *dxs, real *dxsn, real *up,
	real *dspr, real *ddspr, real *ddsprn, real *dsup,
	real *ddsup, real *ddsupn, real *dv, real *ddv,
	real *ddvn, real *prinf, real *upinf, real *duinf,
	real *scale, real *nonzeros,
	fint *vartyp, fint *slktyp, fint *colpnt, fint *ecolpnt,
	fint *count, fint *vcstat, fint *pivots, fint *invprm,
	fint *snhead, fint *nodtyp, fint *inta1, fint *prehis,
	fint *rowidx, fint *rindex, fint *code, real *opt,
	fint *iter, fint *corect, fint *fixn, fint *dropn,
	fint *fnzmax, fint *fnzmin, real *addobj,
	real *bigbou, real *big, fint *ft));

Cextern void timer_ ANSI((fint*));

 void
#ifdef KR_headers
timer_(L) fint *L;
#else
timer_(fint *L)
#endif
{
	*L = (fint)(100.*xectim_());
	}

struct {
    real objnor, rhsnor, scdiff;
    fint scpass, scalmet;
} ascal_ = { 1e2, 0, 1., 5, 2 };

struct {
    real climit, ccorr;
} compl_ = { 1., 1e-5 };

struct {
    real tfixvar, tfixslack, slklim;
} drop_ = { 1e-16, 1e-16, 1e-16 };

struct {
    real tpiv1, tpiv2, tabs, trabs, lam, tfind, order, supdens;
} factor_ = { 1e-3, 1e-8, 1e-12, 1e-15, 1.1e-5, 25., 2.0, 300. };

struct {
    real prmin, upmax, dumin;
    fint stamet, safmet, premet, regul;
} initv_ = { 100., 50000., 100., 2, -3, 511, 0 };

struct {
    real tresx, tresy;
    fint maxref;
} itref_ = { 1e-9, 1e-9, 5 };

struct {
    real varadd, slkadd, scfree;
} mscal_ = { 1e-12, 1e16, 1e-6 };

struct {
    real tplus, tzer;
} numer_ = { 1e-10, 1e-35 };

struct {
    real palpha, dalpha;
} param_ = { 0.999, 0.999 };

struct {
    real target, tsmall, tlarge, center, corstp;
    fint mincc, maxcc;
} predc_ = { 9e-2, 0.2, 20., 5., 1.01, 0, 9 };


struct {
    real ccstop, barset, bargrw, barmin;
    fint mincor, maxcor, inibar;
} predp_ = { 1.01, 0.25, 100., 1e-10, 1, 1, 0 };

struct {
    real maxdense, densgap;
    fint setlam, denslen;
} setden_ = { 0.15, 3.0, 0, 10 };

struct {
    fint psupn, ssupn, maxsnz;
} sprnod_;

struct {
    real tsdir, topt1, topt2, tfeas1, tfeas2, feas1, feas2, pinfs,
	    dinfs, inftol;
    fint maxiter;
} toler_ = {1e-16, 1e-8, 1e-16, 1e-7, 1e-7, 1e-2, 1e-2, 1e-6, 1e-6, 1e4, 99 };

 static char *
#ifdef KR_headers
L1_val(oi, kw, value) Option_Info *oi; keyword *kw; char *value;
#else
L1_val(Option_Info *oi, keyword *kw, char *value)
#endif
{
	char *rv;
	Long L, *Lp = (Long*)kw->info;
	Not_Used(oi);
	L = (Long)strtol(value, &rv, 10) & 0xff;
	*Lp &= ~0xff;
	*Lp |= L;
	return rv;
	}

 static char *
#ifdef KR_headers
L2_val(oi, kw, value) Option_Info *oi; keyword *kw; char *value;
#else
L2_val(Option_Info *oi, keyword *kw, char *value)
#endif
{
	char *rv;
	Long L, *Lp = (Long*)kw->info;
	Not_Used(oi);
	L = (Long)strtol(value, &rv, 10) & 0x7f;
	*Lp &= 0xff;
	*Lp |= L << 8;
	return rv;
	}

 static char *
#ifdef KR_headers
F_val(oi, kw, v) Option_Info *oi; keyword *kw; char *v;
#else
F_val(Option_Info *oi, keyword *kw, char *v)
#endif
{
	char buf[40], *rv;
	fint ignore;

	Not_Used(oi);
	Not_Used(kw);

	for(rv = v; *rv > ' '; rv++);
	if (rv > v)
		readpar_(v, &outlev, &minmax, &ignore, buf, buf, buf, buf, buf,
			&bigb, &ignore, (ftnlen)(rv - v), (ftnlen)40,
			(ftnlen)8, (ftnlen)8, (ftnlen)8, (ftnlen)8);
	return rv;
	}

static keyword keywds[] = {
 KW("bargrow", D_val, &predp_.bargrw, "barrier growth bound"),
 KW("barmin", D_val, &predp_.barmin, "min. barrier threshold"),
 KW("barset", D_val, &predp_.barset, "barrier setup limit"),
 KW("bigbound", D_val, &bigb, "limit on finite bounds and ranges"),
 KW("center", D_val, &predc_.center, "centrality force"),
 KW("complimit", D_val, &compl_.climit, "improve if compl. < complim"),
 KW("comppar", D_val, &compl_.ccorr, "improve if compl. < comppar * avg. compl."),
 KW("corstop", D_val, &predc_.corstp, "correction stop factor"),
 KW("densgap", D_val, &setden_.densgap, "density gap param"),
 KW("denslen", L_val, &setden_.denslen, "nonzeros in \"dense\" columns"),
 KW("dinfs", D_val, &toler_.dinfs, "min. dual infeas. reduction"),
 KW("dslacklim", D_val, &drop_.slklim, "dual slack variable limit"),
 KW("dudare", D_val, &param_.dalpha, "max. dual steplength"),
 KW("dumin", D_val, &initv_.dumin, "min. initial slack value"),
 KW("feas1", D_val, &toler_.feas1, "abs. primal feas. tol."),
 KW("feas2", D_val, &toler_.feas2, "abs. dual feas. tol."),
 KW("inftol", D_val, &toler_.inftol, 0),
 KW("inibarr", L_val, &predp_.inibar, "use initial barrier param."),
 KW("lam", D_val, &factor_.lam, "minimum value of lambda"),
 KW("maxccorr", L_val, &predc_.maxcc, "number of maximum corrections"),
 KW("maxcorr", L_val, &predp_.maxcor, "number of max. corrections"),
 KW("maxdense", D_val, &setden_.maxdense, "max. dense column rate"),
 KW("maxiter", L_val, &toler_.maxiter, "iteration limit"),
 KW("maxref", L_val, &itref_.maxref, "max. refinements"),
 KW("maxsnz", L_val, &sprnod_.maxsnz, "max. nonzeros in one supernode"),
 KW("memadd", I_val, &ma, "increment (default 0) for scratch array len"),
 KW("memfac", D_val, &mf, "factor (default 2) for scratch array len"),
 KW("minccorr", L_val, &predc_.mincc, "number of minimum corrections"),
 KW("mincorr", L_val, &predp_.mincor, "number of min. corrections"),
 KW("mindiff", D_val, &ascal_.scdiff, "min. norm diff. to exit scaling"),
 KW("minmax", L_val, &minmax, "1 = minimize, -1 = maximize"),
 KW("objno", I_val, &nobj, "objective number: 0 = none, 1 (default) = first"),
 KW("objnorm", D_val, &ascal_.objnor, "scale objective maxnorm to this value"),
 KW("ordering", D_val, &factor_.order, "1 = min.deg., 2 (dflt) = min. fill, 3 = min. fill. in aug. sys."),
 KW("outlev", L_val, &outlev, "1 = print iteration summary, etc."),
 KW("parfile", F_val, 0, "name of BPMPD parameter file to read now"),
 KW("pinfs", D_val, &toler_.pinfs, "min. primal infeas. reduction"),
 KW("prdare", D_val, &param_.palpha, "max. primal steplength"),
 KW("presolv", L_val, &initv_.premet, "presolve methods: default = 511 (all but extended dual test)"),
 KW("prmin", D_val, &initv_.prmin, "min. initial variable value"),
 KW("psupnode", L_val, &sprnod_.psupn, "primary supernode length"),
 KW("regularize", L_val, &initv_.regul, "1 = introduce dummy ranges"),
 KW("rhsnorm", D_val, &ascal_.rhsnor, "scale RHS maxnorm to this value"),
 KW("safemet", L_val, &initv_.safmet, "safe method (1, 2, or 3); default -3"),
 KW("scfree", D_val, &mscal_.scfree, 0),
 KW("setlam", L_val, &setden_.setlam, "0 = dflt; 1, -1 = alt. values"),
 KW("slackadd", D_val, &mscal_.slkadd, 0),
 KW("smethod", L_val, &initv_.stamet, "starting method (1 or 2)"),
 KW("smethod1", L1_val, &ascal_.scalmet, "scaling method before aggregator"),
 KW("smethod2", L2_val, &ascal_.scalmet, "scaling method after aggregator"),
 KW("spasses1", L1_val, &ascal_.scpass, "max. passes (< 128) before aggregator"),
 KW("spasses2", L2_val, &ascal_.scpass, "max. passes (< 128) after aggregator"),
 KW("ssupnode", L_val, &sprnod_.ssupn, "secondary supernode length"),
 KW("stopcor", D_val, &predp_.ccstop, "correction stop par."),
 KW("supdens", D_val, &factor_.supdens, "super-dense column length"),
 KW("tabs", D_val, &factor_.tabs, "abs. pivot tol. in first factorization"),
 KW("target", D_val, &predc_.target, "trial steplength improvement"),
 KW("tfeas1", D_val, &toler_.tfeas1, "rel. primal feas. tol."),
 KW("tfeas2", D_val, &toler_.tfeas2, "rel. dual feas. tol."),
 KW("tfind", D_val, &factor_.tfind, "pivot search loop count"),
 KW("tfixslack", D_val, &drop_.tfixslack, "slack reset param."),
 KW("tfixvar", D_val, &drop_.tfixvar, "variable reset param."),
 KW("tlarge", D_val, &predc_.tlarge, "large complementarity bound"),
 KW("topt1", D_val, &toler_.topt1, "rel. duality gap tol."),
 KW("topt2", D_val, &toler_.topt2, "avg. complementarity gap tol."),
 KW("tpiv1", D_val, &factor_.tpiv1, "first threshold pivot tol."),
 KW("tpiv2", D_val, &factor_.tpiv2, "second threshold pivot tol."),
 KW("tplus", D_val, &numer_.tplus, "rel. addition tol."),
 KW("trabs", D_val, &factor_.trabs, "abs. pivot tol. during the algorithm"),
 KW("tresx", D_val, &itref_.tresx, "primal residual tol."),
 KW("tresy", D_val, &itref_.tresy, "dual residual tol."),
 KW("tsdir", D_val, &toler_.tsdir, "search dir. max. norm tol."),
 KW("tsmall", D_val, &predc_.tsmall, "small complementarity bound"),
 KW("tzer", D_val, &numer_.tzer, "rel. zero tol."),
 KW("upmax", D_val, &initv_.upmax, "max. initial variable value"),
 KW("varadd", D_val, &mscal_.varadd, 0),
 KW("wantsol", WS_val, 0, WS_desc_ASL+5)
 };


static Option_Info Oinfo = { "bpmpd", "BPMPD 2.11", "bpmpd_options",
				keywds, nkeywds };

 extern struct {
    fint n, n1, m, mn, nz, cfree, pivotn, denwin, rfree;
} dims_;

 extern struct {
    fint loglog, lfile;
} logprt_;


 extern char **xargv;

 void
#ifdef KR_headers
MAIN__(VOID) VOID;
#else
MAIN__(VOID)
#endif
{
	ASL *asl;
	char buf[256], *stub;
	FILE *nl;
	int cfree, i, j, m, mn, n, n1, nneg, nz, rfree;
	int ip[15], *irn, *neg, rp[25];
	fint *IM, *fcs, *frn;
	fint code, correct, dropn, fixn, fnzmax, fnzmin, ft, iter;
	real *RM, addobj, big, *lb, mbig, objsign, opt, *rhs, *ub, *x, *y;
	ograd *og;
	typedef struct { char *msg; int code; } Sol_info;
	static Sol_info solinfo[] = {
		{ "Execution stopped", 510 },
		{ "Optimal solution found", 0 },
		{ "Dual infeasible (or badly scaled)", 300 },
		{ "Primal infeasible (or badly scaled)", 200 }
		};

	asl = ASL_alloc(ASL_read_f);
	stub = getstops(xargv, &Oinfo);
	nl = jac0dim(stub, (fint)strlen(stub));

	dims_.n = n = n_var;
	dims_.n1 = n1 = n + 1;
	dims_.m = m = n_con;
	dims_.nz = nz = nzc;
	dims_.mn = mn = m + n;
	ip[0] = 0;
	ip[1] = n;
	ip[2] = mn;
	ip[3] = ip[2] + n1;
	for(i = 4; i < 12; i++)
		ip[i] = ip[i-1] + mn;
	ip[12] = ip[11] + nz;
	rp[0] = 0;
	rp[1] = n;
	rp[2] = mn;
	for(i = 3; i < 17; i++)
		rp[i] = rp[i-1] + mn;
	for(; i < 21; i++)
		rp[i] = rp[i-1] + m;
	for(; i < 24; i++)
		rp[i] = rp[i-1] + mn;

	if (mf < 2.)
		mf = 2.;
	if (ma < 0)
		ma = 0;
	dims_.cfree = cfree = mf*(rp[23] + nz) + ma;
	dims_.rfree = rfree = ((cfree + rp[23])>>1) + 11*m + 8*n;

	rp[24] = rp[23] + cfree;	/* ultimate length of real data */
	ip[13] = ip[12] + cfree;
	ip[14] = ip[13] + rfree;	/* ultimate length of integer data */

	IM = (fint*)Malloc(ip[14]*sizeof(fint));
	fcs = IM + ip[2];
	frn = IM + ip[12];

	if (sizeof(int) == sizeof(fint)) {
		A_colstarts = (int*)fcs;
		A_rownos = (int*)frn;
		}

	RM = (real*)Malloc(rp[24]*sizeof(real));
	A_vals = RM + rp[23];
	LUv = RM + rp[2];
	LUrhs = LUv + n;
	Uvx = RM + rp[3];
	Urhsx = Uvx + n;
	Fortran = 1;

	f_read(nl,0);

	if (sizeof(int) != sizeof(fint)) {
		irn = A_colstarts;
		for(i = 0; i < n1; i++)
			fcs[i] = irn[i];
		irn = A_rownos;
		for(i = 0; i < nz; i++)
			frn[i] = irn[i];
		}

	/* cost vector comes first in RM */

	objsign = addobj = 0;
	if (--nobj >= 0 && nobj < n_obj) {
		if (minmax)
			objsign = minmax > 0 ? 1. : -1.;
		else
			objsign = objtype[nobj] ? -1. : 1.;
		addobj = objsign * objconst(nobj);
		for(og = Ograd[nobj]; og; og = og->next)
			RM[og->varno] = objsign*og->coef;
		}

	/* compute rhs from lb and ub */

	rhs = RM + rp[1];
	lb =  RM + rp[2] + n;
	ub =  RM + rp[3] + n;
	for(i = nneg = 0; i < m; i++) {
		if (lb[i] > negInfinity) {
			ub[i] -= (rhs[i] = lb[i]);
			lb[i] = 0;
			}
		else if (ub[i] < Infinity) {
			if (!nneg)
				neg = (int*)Malloc((m-i)*sizeof(int));
			neg[nneg++] = i;
			rhs[i] = ub[i];
			ub[i] = 0;
			}
		}

	logprt_.loglog = outlev > 0;
	big = 1e30;
	mbig = -big;

	lb -= n;
	ub -= n;
	for(i = 0; i < mn; i++) {
		if (lb[i] < mbig)
			lb[i] = mbig;
		if (ub[i] > big)
			ub[i] = big;
		}

	solver_(RM+rp[0],  RM+rp[1],  RM+rp[2],  RM+rp[3],  RM+rp[4],
		RM+rp[5],  RM+rp[6],  RM+rp[7],  RM+rp[8],  RM+rp[9],
		RM+rp[10], RM+rp[11], RM+rp[12], RM+rp[13], RM+rp[14],
		RM+rp[15], RM+rp[16], RM+rp[17], RM+rp[18], RM+rp[19],
		RM+rp[20], RM+rp[21], RM+rp[22], RM+rp[23], IM+ip[0],
		IM+ip[1],  IM+ip[2],  IM+ip[3],  IM+ip[4],  IM+ip[5],
		IM+ip[6],  IM+ip[7],  IM+ip[8],  IM+ip[9],  IM+ip[10],
		IM+ip[11], IM+ip[12], IM+ip[13], &code,	    &opt,
		&iter,	   &correct,  &fixn,	 &dropn,    &fnzmax,
		&fnzmin,   &addobj,   &bigb,	 &big,	    &ft);

	x = y = 0;
	i = Sprintf(buf, "%s: ", Oinfo.bsname);
	if (code == -2) {
		i += sprintf(buf+i,
			"Ran out of memory: current mf = %g, ma = %d",
			mf, ma);
		solve_result_num = 520;
		}
	else if (--code < 0 || code >= 4) {
		i += sprintf(buf+i, "Unexpected return code %ld", code+1);
		solve_result_num = 500;
		}
	else {
		i += sprintf(buf+i, solinfo[code].msg);
		solve_result_num = solinfo[code].code;
		if (code) {
			x = RM + rp[6];
			y = RM + rp[16];
			switch(code) {
			  case 1: i += sprintf(buf+i, ", objective %.*g",
					obj_prec(), objsign*opt);
				  break;
			  case 2: y = 0; break;
			  case 3: x = 0;
			  }
			if (y) {
				while(nneg-- > 0)
					y[*neg++] *= -1.;
				if (objsign < 0)
					for(j = 0; j < m; j++)
						y[j] *= -1.;
				if (dropn) {
					i += sprintf(buf+i, "\n%ld%s\n%s",
						(long)dropn,
		" dual value(s) may be wrong; to get correct",
		"dual values, add \" presolv=0\" to $bpmpd_options.");
					if (!solve_result_num)
						solve_result_num = 1;
					}
				}
			}
		i += sprintf(buf+i, "\n%ld iterations, %ld corrections",
			(long)iter, (long)correct);
		if (iter > 1)
			sprintf(buf+i, " (%.2f per iter.)",
				(double)correct / iter);
		}
	write_sol(buf, x, y, &Oinfo);
	}
