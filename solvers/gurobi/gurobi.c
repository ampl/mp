/****************************************************************
Copyright (C) 2009, 2012 AMPL Optimization LLC
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

AMPL Optimization LLC DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS
SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS.  IN NO EVENT SHALL AMPL Optimization LLC OR ANY OF ITS
ENTITIES BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES
OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,
ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS
SOFTWARE.
****************************************************************/
#include "gurobi_c.h"
#include "getstub.h"
#include "signal.h"
#include "time.h"

#ifndef Sig_ret_type
#define Sig_ret_type void
#define SigRet /*nothing*/
#endif

#ifndef SigRet
#define SigRet /*nothing*/
#endif

#ifndef Sig_func_type
typedef void sig_func_type(int);
#endif

#undef ALLOW_GUROBI_SERVER
#if GRB_VERSION_MAJOR > 5 || (GRB_VERSION_MAJOR == 5 && GRB_VERSION_MINOR >= 5)
#ifndef DISALLOW_GUROBI_SERVER
#define ALLOW_GUROBI_SERVER
#endif
#endif

 typedef struct
Dims {
	double	*c;
	double	*x;
	double	*y;
	double	*y0;
	int	*cstat;
	int	*rstat;
	SufDesc	*csd;
	SufDesc *rsd;
	char	*mb, *mbend;
	int	kiv;
	int	missing;
	int	nc0, nv0;
	int	objsense;
	} Dims;

 typedef struct
Filename {
	struct Filename *next;
	char *name;
	} Filename;

 typedef struct
Ext_info {
	char *ext;
	int aftersol;
	} Ext_info;

 typedef struct
mint_values {
	int L;
	int U;
	int val;
	} mint_values;

 enum { /* sf_mint f values */
	set_iis		= 0,
	set_relax	= 1,
	set_mipstval	= 2,
	set_objno	= 3,
	set_sos		= 4,
	set_sos2	= 5,
	set_timing	= 6,
	set_basis	= 7,
	set_intstart	= 8,
	set_outlev	= 9,
	set_bestbound	= 10,
	set_solnsens	= 11,
	set_retmipgap	= 12,
	set_rays	= 13,
	set_priorities	= 14,
	set_feasrelax	= 15,
	set_warmstart	= 16
	};

 static mint_values
mint_val[17] = {
	/* set_iis */		{0, 1, 0},
	/* set_relax */		{0, 1, 0},
	/* set_mipstval */	{0, 1, 1},
	/* set_objno */		{0, 0/*n_obj*/,	1},
	/* set_sos */		{0, 1, 1},
	/* set_sos2 */		{0, 1, 1},
	/* set_timing */	{0, 3, 0},
	/* set_basis */		{0, 3, 3},
	/* set_intstart */	{0, 1, 1},
	/* set_outlev */	{0, 1, 0},
	/* set_bestbound */	{0, 1, 0},
	/* set_solnsens */	{0, 1, 0},
	/* set_retmipgap */	{0, 7, 0},
	/* set_rays */		{0, 3, 3},
	/* set_priorities */	{0, 1, 1},
	/* set_feasrelax */	{0, 6, 0},
	/* set_warmstart */	{0, 2, 1}
	};

#define want_iis	mint_val[0].val
#define relax		mint_val[1].val
#define mipstval	mint_val[2].val
#define nobjno		mint_val[3].U
#define objno		mint_val[3].val
#define sos		mint_val[4].val
#define sos2		mint_val[5].val
#define time_flag	mint_val[6].val
#define basis		mint_val[7].val
#define intstart	mint_val[8].val
#define outlev		mint_val[9].val
#define bestbound	mint_val[10].val
#define solnsens	mint_val[11].val
#define retmipgap	mint_val[12].val
#define rays		mint_val[13].val
#define priorities	mint_val[14].val
#define feasrelax	mint_val[15].val
#define warmstart	mint_val[16].val

 static Filename *Wflist[3];
 static GRBmodel *grbmodel;
 static char *logfile, verbuf[64];
 static double Times[5];
 static int breaking, wantlog;
 static jmp_buf Jb;

#ifdef ALLOW_GUROBI_SERVER /*{*/
 static char *server, *server_passwd, *serverlic;
 static int server_port = DEFAULT_CS_PORT, server_priority = DEFAULT_CS_PRIORITY;
 static double server_timeout = -1.;
 static char
	server_desc[] = "Comma-separated list of Gurobi compute servers, specified\n\
		either by name or by IP address.  Default: run Gurobi locally\n\
		(i.e., do not use a remote Gurobi server).",

	server_passwd_desc[] = "Password (if needed) for specified Guruobi compute server(s).",

	server_port_desc[] = "IP port to use for Gurobi compute server(s);\n\
		-1 ==> use default.",

	server_priority_desc[] = "Priority for Gurobi compute server(s).  Default = 1.\n\
		Highest priority = 100.",

	server_timeout_desc[] = "Report job as rejected by Gurobi compute server if the\n\
		job is not started within server_timeout seconds.\n\
		Default = -1 (no limit).",

	serverlic_desc[] = "Name of file containing \"server = ...\" and possibly\n\
		values for server_password, server_port, and server_timeout";

 static int
badserverlic(ASL *asl, const char *what, int rn)
{
	char *s;
	size_t L;

	L = strlen(serverlic) + strlen(what) + 32;
	asl->i.uinfo = s = M1alloc(L);
	snprintf(s, L, "%s serverlic file \"%s\".", what, serverlic);
	solve_result_num = rn;
	return 1;
	}

 static int
server_licread(ASL *asl)
{
	FILE *f;
	char buf[4096], *s;

	static keyword lrkeywds[] = {
	{ "server", C_val, &server, server_desc },
	{ "server_password", C_val, &server_passwd, server_passwd_desc },
	{ "server_port", I_val, &server_port, server_port_desc },
	{ "server_priority", I_val, &server_priority, server_priority_desc },
	{ "server_timeout", D_val, &server_timeout, server_timeout_desc } };
	static Option_Info lrinfo = { 0,0,0, lrkeywds,
		(int)(sizeof(lrkeywds)/sizeof(keyword)), ASL_OI_never_echo };

	if (!(f = fopen(serverlic, "r")))
		return badserverlic(asl, "Cannot open", 530);
	lrinfo.asl = asl;
 nextline:
	while(lrinfo.n_badopts == 0 && fgets(s = buf, sizeof(buf), f)) {
		while(*s && lrinfo.n_badopts == 0) {
			while(*s <= ' ') {
				if (!*s++)
					goto nextline;
				}
			if (*s == '#')
				goto nextline;
			s = get_opt_ASL(&lrinfo, s);
			}
		}
	fclose(f);
	if (lrinfo.n_badopts)
		return badserverlic(asl, "Bad assignment in", 531);
	return 0;
	}
#endif /*}*/

#if GRB_VERSION_MAJOR >= 3
 static double ams_eps, ams_epsabs;
 static int ams_limit;
 static char *ams_stub;
#endif

#if GRB_VERSION_MAJOR >= 5 /*{*/
 static real lbpen = 1., rhspen = 1., ubpen = 1.;
#ifdef GRB_INT_PAR_TUNEOUTPUT /*{*/
  static char *tunebase;

 static int
tunerun(ASL *asl, GRBenv *env, GRBmodel *mdl, char **tunemsg)
{
	char *b, *s, *tbuf, tbuf0[4096];
	const char *fmt, *fmt2, *em, *what;
	int i, j, k, m, n;
	size_t L;
	static char prm[] = ".prm";

	what = "GRBtunemodel";
	fmt2 = "Return %d from %s: %s().";
	L = 64;
	em = 0;
	if ((i = GRBtunemodel(mdl))) {
 trouble:
		fmt = "Return %d from %s().";
		if ((em = GRBgeterrormsg(env))) {
			L += strlen(em);
			fmt = fmt2;
			}
 trouble2:
		*tunemsg = s = (char*)M1alloc(L);
		snprintf(s, L, fmt, i, what, em);
 badret:
		solve_result_num = 532;
		return 1;
		}
	n = -1;
	fmt = "No tuning results available: return %d from %s().";
	fmt2 = "No tuning results available: return %d from %s():\n%s.";
	L = 80;
	if ((i = GRBgetintattr(mdl, GRB_INT_ATTR_TUNE_RESULTCOUNT, &n))) {
		what = "GRBgetintattr";
		goto trouble;
		}
	if (n <= 0) {
		L = 32;
		fmt = "No tuning results found.";
		goto trouble2;
		}
	for(s = tunebase, b = 0; *s; ++s) {
		if (*s == '.' && s[1] == 'p' && s[2] == 'r' && s[3] == 'm') {
			b = s;
			s += 3;
			}
		}
	if (b)
		j = b - tunebase;
	else {
		j = s - tunebase;
		b = prm;
		}
	tbuf = tbuf0;
	if ((L = s - tunebase + 16) > sizeof(tbuf0))
		tbuf = (char*)M1alloc(L);
	if (j > 0)
		memcpy(tbuf, tunebase, j);
	for(k = m = 0; n > 0; ++k) {
		if ((i = GRBgettuneresult(mdl, --n))) {
			what = "GRBgettuneresult";
			if (!k)
				goto trouble;
			L = L + 80;
			*tunemsg = s = (char*)M1alloc(L);
			snprintf(s, L, "Surprise return %d from %s() after writing"
				" %d %.*s*%s files.", i, what, k, j, tunebase, b);
			goto badret;
			}
		m = snprintf(tbuf+j, L-j, "%d%s", n+1, b);
		what = "GRBwriteparams";
		if ((i = GRBwriteparams(env, tbuf)))
			goto trouble2;
		}
	*tunemsg = s = (char*)M1alloc(L = 2*m + 64);
	switch(k) {
	 case 1:
		snprintf(s, L, "Wrote tuning parameter file \"%s\".", tbuf);
		break;
	 case 2:
		snprintf(s, L,
			"Wrote tuning parameter files \"%s\" and \"%.*s2%s\".",
			tbuf, j, tunebase, b);
		break;
	 default:
		snprintf(s, L,
			"Wrote %d tuning parameter files \"%s\" ... \"%.*s%d%s\".",
			k, tbuf, j, tunebase, k, b);
	 }
	return 0;
	}
#endif /*}*/
#endif /*}*/

 static void
badretfmt(int rc, char *fmt, ...)
{
	ASL *asl = cur_ASL;
	va_list ap;
	char buf[8192], *s;
	int k;

	va_start(ap, fmt);
	k = Vsnprintf(buf, sizeof(buf)-1, fmt, ap) + 1;
	if (rc) {
		solve_result_num = rc;
		memcpy(s = (char*)M1alloc(k), buf, k);
		asl->i.uinfo = s;
		}
	else
		fprintf(Stderr, "%s\n", buf);
	va_end(ap);
	}

 static void
failed(GRBenv *env, const char *what)
{
	const char *s = GRBgeterrormsg(env);
	if (s)
		badretfmt(501, "%s failed:\n\t%s.\n", what, s);
	else
		badretfmt(501, "%s failed.\n", what);
	longjmp(Jb,1);
	}

 static void enamefailed(GRBenv *env, const char *what, const char *name);

 static void
namefailed(const char *what, const char *name)
{
	badretfmt(506, "%s(\"%s\") failed.", what, name);
	longjmp(Jb,1);
	}

 static void
badival(Option_Info *oi, keyword *kw, int t, int L, int U)
{
	printf("rejecting %s %d; must be between %d and %d\n",
		kw->name, t, L, U);
	badopt_ASL(oi);
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
	return rv;
	}

 static char *
sf_dpar(Option_Info *oi, keyword *kw, char *v)
{
	GRBenv *env;
	double p[4], t;
	char *parname, *rv;

	env = (GRBenv*)oi->uinfo;
	parname = (char*)kw->info;

	if (*v == '?' && v[1] <= ' ') {
		if (GRBgetdblparam(env, parname, &t))
			namefailed("GRBgetdblparam", parname);
		printf("%s=%.g\n", kw->name, t);
		oi->option_echo &= ~ASL_OI_echothis;
		return v + 1;
		}
	t = strtod(v, &rv);
	if (rv == v) {
		printf("Expected a numeric value for %s, not \"%s\"\n",
			kw->name, v);
		badopt_ASL(oi);
		return v;
		}
	if (GRBsetdblparam(env, parname, t)) {
		if (GRBgetdblparaminfo(env, parname, p, p+1, p+2, p+3))
			namefailed("GRBsetdblparam", parname);
		badretfmt(506, "%s must be >= %.g and <= %.g.", kw->name, p[1], p[2]);
		badopt_ASL(oi);
		}
	return rv;
	}

 static void
int_rangerr(Option_Info *oi, keyword *kw)
{
	GRBenv *env;
	char *fmt, *parname;
	int p[4];

	env = (GRBenv*)oi->uinfo;
	parname = (char*)kw->info;
	if (GRBgetintparaminfo(env, parname, p, p+1, p+2, p+3))
		namefailed("GRBsetintparam", parname);
	fmt = p[2] == p[1] + 1
		? "%s must be %d or %d."
		: "%s must be >= %d and <= %d.";
	badretfmt(506, fmt, kw->name, p[1], p[2]);
	badopt_ASL(oi);
	}

 static char *
sf_ipar(Option_Info *oi, keyword *kw, char *v)
{
	GRBenv *env;
	int t;
	char *parname, *rv;

	env = (GRBenv*)oi->uinfo;
	parname = (char*)kw->info;

	if (*v == '?' && v[1] <= ' ') {
		if (GRBgetintparam(env, parname, &t))
			namefailed("GRBgetintparam", parname);
		printf("%s=%d\n", kw->name, t);
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
	if (GRBsetintparam(env, parname, t))
		int_rangerr(oi, kw);
	return rv;
	}

 static char*
sf_iparlog(Option_Info *oi, keyword *kw, char *v)
{
	GRBenv *env;
	int t;
	char *parname, *rv;

	env = (GRBenv*)oi->uinfo;
	parname = (char*)kw->info;

	if (*v == '?' && v[1] <= ' ') {
		if (GRBgetintparam(env, parname, &t))
			namefailed("GRBgetintparam", parname);
		printf("%s=%d\n", kw->name, t);
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
	if (GRBsetintparam(env, parname, t))
		int_rangerr(oi, kw);
	else if (t)
		++wantlog;
	return rv;
	}

 static Filename *
fn_value(char **pv, const char *what, keyword *kw)
{
	ASL *asl = cur_ASL;
	Filename *f;
	char *s, *t, *v;
	int c, q;
	size_t L;

	v = *pv;
	q = *v;
	if (q == '"' || q == '\'') {
		s = ++v;
		for(;;) {
			if (!(c = *s)) {
				printf("Bad %s \"%s\" for %s\n", what, v, kw->name);
				*pv = v;
				return 0;
				}
			++s;
			if (c == q && *s != q)
				break;
			}
		}
	else {
		q = 0;
		for(s = v; *s > ' '; ++s);
		if (s == v) {
			printf("Missing %s for %s\n", what, kw->name);
			return 0;
			}
		}
	L = s - v;
	f = M1alloc(sizeof(Filename) + L + 1);
	f->name = t = (char*)(f + 1);
	if (q) {
		for(s = v;; ++s, ++t) {
			if ((*t = *s) == q) {
				if (*++s == q)
					continue;
				break;
				}
			}
		}
	else {
		memcpy(t, v, L);
		t += L;
		}
	*t = 0;
	*pv = s;
	return f;
	}

#if GRB_VERSION_MAJOR > 1 /*{*/
#define GRB_MAJ2(x) x

 static char *
sf_spar(Option_Info *oi, keyword *kw, char *v)
{
	Filename *f;
	GRBenv *env;
	char *parname, tbuf[GRB_MAX_STRLEN + 8];

	env = (GRBenv*)oi->uinfo;
	parname = (char*)kw->info;

	if (*v == '?' && v[1] <= ' ') {
		memset(tbuf, 0, sizeof(tbuf));
		if (GRBgetstrparam(env, parname, tbuf))
			namefailed("GRBgetstrparam", parname);
		printf("%s=\"%s\"\n", kw->name, tbuf);
		oi->option_echo &= ~ASL_OI_echothis;
		return v + 1;
		}
	if ((f = fn_value(&v, "value", kw))
	 && GRBsetstrparam(env, parname, f->name))
		enamefailed(env, "GRBsetstrparam", parname);
	return v;
	}
#else
#define GRB_MAJ2(x) /*nothing*/
#endif /*}*/

 static char *
sf_wfile(Option_Info *oi, keyword *kw, char *v)
{
	Ext_info *e;
	Filename *f, **pf;
	char *dot, *t;
	int i, q;
	static Ext_info W_ext[] = {
		{"bas",1},
		{"fix_lp",2},
		{"fix_mps",2},
		{"lp",0},
		{"mps",0},
		{"prm",0},
		{"sol",1},
		{0,0}};

	q = *v;
	if (q == '?' && v[1] <= ' ') {
		for(i = 0; i < 3; ++i)
			for(f = Wflist[i]; f; f = f->next)
				printf("%s=\"%s\"\n", kw->name, f->name);
		oi->option_echo &= ~ASL_OI_echothis;
		return v + 1;
		}
	if (!(f = fn_value(&v, "file name", kw))) {
		printf("Bad file name \"%s\" for %s\n", v, kw->name);
		badopt_ASL(oi);
		return v;
		}
	dot = 0;
	for(t = f->name; *t; )
		if (*t++ == '.')
			dot = t;
	if (dot)
		for(e = W_ext; e->ext; ++e)
			if (!strcmp(e->ext, dot))
				goto good_ext;
	printf("File name for %s must end in one of\n", kw->name);
	for(e = W_ext; e->ext; ++e)
		printf("\t.%s\n", e->ext);
	badopt_ASL(oi);
	goto ret;
 good_ext:
	pf = &Wflist[e->aftersol];
	f->next = *pf;
	*pf = f;
	if (e->aftersol == 2)	/* replace _ by . */
		dot[3] = '.';
 ret:
	return v;
	}

#if GRB_VERSION_MAJOR >= 3 /*{*/

 static char *
sf_pf(Option_Info *oi, keyword *kw, char *v)
{
	FILE *f;
	Filename *fn;
	GRBenv *env;
	char buf[512], *fname, *s, *s1, *se;
	int lineno;
	static char extra[]  = "Line %d of paramfile \"%s\":\n\t"
		"expected a name and value, but got \"%s\".";
	static char failed[] = "Line %d of paramfile \"%s\":\n\t"
		"GRBsetstrparam(\"Dummy\", \"%s\") failed:\n\t%s.";
	static char missing[] = "Missing value in line %d of paramfile \"%s\".";

	env = (GRBenv*)oi->uinfo;

	if ((fn = fn_value(&v, "value", kw))) {
		fname = fn->name;
		if (!(f = fopen(fname, "r"))) {
			badretfmt(511, "Cannot open paramfile \"%s\".", fname);
			longjmp(Jb,1);
			}
		lineno = 0;
 nextline:
		while(fgets(buf, sizeof(buf), f)) {
			++lineno;
			for(s = buf; *s <= ' '; ++s)
				if (!*s)
					goto nextline;
			if (*s == '#')
				goto nextline;
			for(s1 = s; *++s1 > ' '; );
			while(*s1 <= ' ')
				if (!*s1++) {
					badretfmt(512, missing, lineno, fname);
					longjmp(Jb,1);
					}
			for(se = s1; *++se; );
			while(--se > s1 && *se <= ' ');
			se[1] = 0;
			while(*++s1 > ' ');
			if (*s1) {
				for(se = s1; *++se; ) {
					if (*se > ' ') {
						badretfmt(513, extra, lineno, fname, s);
						longjmp(Jb,1);
						}
					}
				*s1 = 0;
				}
			if (GRBsetstrparam(env, "Dummy", s)) {
				badretfmt(514, failed, lineno, fname, s, GRBgeterrormsg(env));
				longjmp(Jb,1);
				}
			}
		fclose(f);
		}
	return v;
	}

 static char aggfill_desc[] = "amount of fill allowed during aggregation during\n\
			gurobi's presolve "
#if GRB_VERSION_MAJOR > 5 || (GRB_VERSION_MAJOR == 5 && GRB_VERSION_MINOR >= 1)
	"(default -1)";
#else
	"(default 10)";
#endif
#endif /*}*/

 static char aggregate_desc[] = "whether to use aggregation during Gurobi presolve:\n\
			0 = no (sometimes reduces numerical errors)\n\
			1 = yes (default)";

#if GRB_VERSION_MAJOR >= 3 /*{*/
 static char ams_eps_desc[] = "relative tolerance for reporting alternate MIP solutions\n\
			(default = no limit)";

 static char ams_epsabs_desc[] = "absolute tolerance for reporting alternate MIP solutions\n\
			(default = no limit)";

 static char ams_limit_desc[] = "limit on number of alternate MIP solutions written\n\
			(default = number of available alternate solutions)";

 static char ams_stub_desc[] = "stub for alternate MIP solutions.  The number of\n\
			alternative MIP solution files written is determined\n\
			by three keywords:\n\
			  ams_limit gives the maximum number of files written;\n\
			  ams_eps gives a relative tolerance on the objective\n\
				values of alternative solutions; and\n\
			  ams_epsabs gives an absolute tolerance on how much\n\
				worse the objectives can be.";

 static char barconvtol_desc[] = "tolerance on the relative difference between the\n\
			primal and dual objectives for stopping the barrier\n\
			algorithm (default 1e-8)";

 static char barcorrectors_desc[] = "Limit on the number of central corrections done in\n\
			each barrier iteration (default -1 = automatic choice)";

#ifdef GRB_INT_PAR_BARHOMOGENEOUS /*{*/
 static char barhomogeneous_desc[] = "Whether to use the homogeneous barrier algorithm\n\
		(e.g., when method=2 is specified):\n\
			-1 = only when solving a MIP node relaxation (default)\n\
			 0 = never\n\
			 1 = always.\n\
		The homogeneous barrier algorithm can detect infeasibility or\n\
		unboundedness directly, without crossover, but is a bit slower\n\
		than the nonhomogeneous barrier algorithm.";
#endif /*}*/

 static char bariterlimit_desc[] = "Limit on the number of barrier iterations (default none)";

 static char barorder_desc[] = "Ordering used to reduce fill in sparse-matrix factorizations\n\
				during the barrier algorithm:\n\
		       -1 = automatic choice\n\
			0 = approximate minimum degree\n\
			1 = nested dissection";
#endif /*}*/

#ifdef GRB_DBL_PAR_BARQCPCONVTOL
 static char barqcptol_desc[] = "convergence tolerance on the relative difference between\n\
		primal and dual objective values for barrier algorithms when\n\
		solving problems with quadratic constraints (default 1e-6)";
#endif

 static char basis_desc[] = "whether to use or return a basis:\n\
			0 = no\n\
			1 = use incoming basis (if provided)\n\
			2 = return final basis\n\
			3 = both (1 + 2 = default)"
#if GRB_VERSION_MAJOR >= 5
			"\n\
		For problems with integer variables and quadratic constraints,\n\
		basis = 0 is assumed quietly unless qcpdual=1 is specified."
#endif
;

 static char bestbound_desc[] = "whether to return suffix .bestbound for the\n\
		best known bound on the objective value:\n\
			0 = no (default)\n\
			1 = yes";

#ifdef GRB_INT_PAR_CONCURRENTMIP
 static char concurrentmip_desc[] =
		"how many independent MIP solves to allow at once when multiple\n\
		threads are available.  The available threads are divided as\n\
		evenly as possible among the concurrent solves.  Default = 1.";
#endif

#if GRB_VERSION_MAJOR >= 3 /*{*/
 static char crossover_desc[] = "how to transform a barrier solution to a basic one:\n\
		       -1 = automatic choice (default)\n\
			0 = none: return an interior solution\n\
			1 = push dual vars first, finish with primal simplex\n\
			2 = push dual vars first, finish with dual simplex\n\
			3 = push primal vars first, finish with primal simplex\n\
			4 = push primal vars first, finish with dual simplex";
 static char crossoverbasis_desc[] = "strategy for initial basis construction during crossover:\n\
			0 = favor speed (default)\n\
			1 = favor numerical stability";
#endif /*}*/

#if (GRB_VERSION_MAJOR == 4 && GRB_VERSION_MINOR >= 5) || GRB_VERSION_MAJOR >= 5 /*{*/
 static char branchdir_desc[] = "Which child node to explore first when branching:\n\
			-1 = explore \"down\" branch first\n\
			 0 = explore \"most promising\" branch first (default)\n\
			 1 = explore \"up\" branch first";
#endif /*}*/
 static char cutagg_desc[] = "maximum number of constraint aggregation passes\n\
		during cut generation (-1 = default = no limit);\n\
		overrides \"cuts\"";

 static char cutoff_desc[] = "If the optimal objective value is no better than cutoff,\n\
		report \"objective cutoff\" and do not return a solution.\n\
		Default: -Infinity for minimizing, +Infinity for maximizing.";

#if (GRB_VERSION_MAJOR == 4 && GRB_VERSION_MINOR >= 5) || GRB_VERSION_MAJOR >= 5 /*{*/
 static char cutpasses_desc[] = "maximum number of cutting-plane passes to do\n\
		during root-cut generation; default = -1 ==> automatic choice";
#endif /*}*/

#ifdef GRB_INT_PAR_DUALREDUCTIONS  /*{*/
 static char dualreductions_desc[] = "whether Gurobi's presolve should use dual reductions, which\n\
		may be useful on a well-posed problem but can prevent\n\
		distinguishing whether a problem is infeasible or unbounded:\n\
			0 = no\n\
			1 = yes (default)";
#endif /*}*/

 static char cuts_desc[] = "global cut generation control, valid unless overridden\n\
		by individual cut-type controls:\n\
		       -1 = automatic choice (default)\n\
			0 = no cuts\n\
			1 = conservative cut generation\n\
			2 = aggressive cut generation"
			GRB_MAJ2("\n\t\t\t3 = very aggressive cut generation")
		;

#if GRB_VERSION_MAJOR >= 5 /*{*/
 static char feasrelax_desc[] = "Whether to modify the problem into a feasibility\n\
		relaxation problem:\n\
			0 = no (default)\n\
			1 = yes, minimizing the weighted sum of violations\n\
			2 = yes, minimizing the weighted count of violations\n\
			3 = yes, minimizing the sum of squared violations\n\
			4-6 = same objective as 1-3, but also optimize the\n\
				original objective, subject to the violation\n\
				objective being minimized\n\
		Weights are given by suffixes .lbpen and .ubpen on variables\n\
		and .rhspen on constraints (when positive), else by keywords\n\
		lbpen, ubpen, and rhspen, respectively (default values = 1).\n\
		Weights <= 0 are treated as Infinity, allowing no violation.";
#ifdef GRB_DBL_PAR_FEASRELAXBIGM
 static char feasrelaxbigm_desc[] =
		"Value of \"big-M\" sometimes used with constraints when doing\n\
		a feasibility relaxation.  Default = 1e6.";
#endif
#endif /*}*/

 static char feastol_desc[] = "primal feasibility tolerance (default 1e-6)";

 static char gomory_desc[] = "maximum number of Gomory cut passes during cut generation\n\
		(-1 = default = no limit); overrides \"cuts\"";

 static char heurfrac_desc[] = "fraction of time to spend in MIP heuristics (default 0.05)";

 static char iisfind_desc[] = "whether to return an IIS (via suffix .iis) when\n\
		the problem is infeasible:\n\
			0 = no (default)\n\
			1 ==> yes";

#if GRB_VERSION_MAJOR > 1 /*{*/
 static char iismethod_desc[] = "which method to use when finding an IIS (irreducible\n\
		infeasible set of constraints, including variable bounds):\n\
		       -1 = automatic choice (default)\n\
			0 = often faster than method 1\n\
			1 = can find a smaller IIS than method 0";
#endif /*}*/

#ifdef GRB_DBL_PAR_IMPROVESTARTGAP
 static char isg_desc[] = "optimality gap below which the MIP solver switches from\n\
		trying to improve the best bound to trying to find better\n\
		feasible solutions (default 0)";
#endif

#ifdef GRB_DBL_PAR_IMPROVESTARTNODES
 static char isn_desc[] = "number of MIP nodes after which the solution strategy\n\
		will change from improving the best bound to finding better\n\
		feasible solutions (default 0)";
#endif

#ifdef GRB_DBL_PAR_IMPROVESTARTTIME
 static char ist_desc[] = "execution seconds after which the MIP solver switches from\n\
		trying to improve the best bound to trying to find better\n\
		feasible solutions (default Infinity)";
#endif

 static char intfeastol_desc[] = "feasibility tolerance for integer variables (default 1e-05)";

 static char intstart_desc[] = "when there are integer variables, whether to use\n\
		an initial guess (if available):\n\
			0 = no\n\
			1 = yes (default)";

 static char logfile_desc[] = "name of file to receive log lines (default: none)"
		GRB_MAJ2(";\n\t\t\timplies outlev = 1")
		;

 static char logfreq_desc[] = "interval in seconds between log lines (default 5)";

 static char maxmipsub_desc[] = "maximum number of nodes for RIMS heuristic to explore\n\
		on MIP problems (default 500)";

#if (GRB_VERSION_MAJOR == 4 && GRB_VERSION_MINOR >= 5) || GRB_VERSION_MAJOR >= 5 /*{*/
 static char minrelnodes_desc[] = "number of nodes for the Minimum Relaxation heuristic\n\
		to explore at the MIP root node when a feasible solution has\n\
		not been found by any other heuristic (default 0)";
#endif /*}*/
#if GRB_VERSION_MAJOR >= 3 /*{*/
 static char mipfocus_desc[]  = "MIP solution strategy:\n\
			0 = balance finding good feasible solutions and\n\
			    proving optimality (default)\n\
			1 = favor finding feasible solutions\n\
			2 = favor proving optimality\n\
			3 = focus on improving the best objective bound";
#endif /*}*/

 static char mipstart_desc[] = "whether to use initial guesses in problems with\n\
		integer variables:\n\
			0 = no\n\
			1 = yes (default)";

#if GRB_VERSION_MAJOR >= 4 /*{*/
 static char nodemethod_desc[] = "algorithm used to solve relaxed MIP node problems:\n\
			0 = primal simplex\n\
			1 = dual simplex (default)\n\
			2 = barrier";
#endif /*}*/

#ifdef GRB_INT_PAR_NUMERICFOCUS
 static char numericfocus_desc[] = "how much to try detecting and managing numerical issues:\n\
			0 = automatic choice (default)\n\
			1-3 = increasing focus on more stable computations";
#endif


 static char objno_desc[] = "objective to optimize:\n\
			0 = none\n\
			1 = first (default, if available),\n\
			2 = second (if available), etc.";

#if GRB_VERSION_MAJOR > 1 /*{*/
 static char objscale_desc[] = "how to scale the objective:\n\
			0 ==> automatic choice (default)\n\
			negative >= -1 ==> divide by max abs. coefficient\n\
					   raised to this power\n\
			positive ==> divide by this value";
#endif /*}*/

 static char opttol_desc[] = "optimality tolerance on reduced costs (default 1e-6)";

 static char outlev_desc[] = "whether to write Gurobi log lines (chatter) to stdout:\n\
			0 = no (default)\n\
			1 = yes (see logfreq)";

#define Overrides_cuts "overrides \"cuts\"; choices as for \"cuts\""
 static char overrides_cuts[] = Overrides_cuts;

 static char perturb_desc[] = "magnitude of simplex perturbation (when needed; default 2e-4)";

#if GRB_VERSION_MAJOR >= 3 /*{*/
 static char param_desc[] = "general way to specify values of both documented and\n\
		undocumented Gurobi parameters; value should be a quoted string\n\
		(delimited by ' or \") containing a parameter name, a space, and\n\
		the value to be assigned to the parameter.  Can appear more\n\
		than once.  Cannot be used to query current parameter values.";

 static char paramfile_desc[] = "name of file (surrounded by 'single' or \"double\" quotes if the\n\
		name contains blanks) of parameter names and values for them.\n\
		Lines that start with # are ignored.  Otherwise, each nonempty\n\
		line should contain a name and a value, separated by a space.";

 static char predeprow_desc[] = "whether Gurobi's presolve should remove linearly\n\
		dependent constraint-matrix rows:\n\
		       -1 = only for continuous models\n\
			0 = never\n\
			1 = for all models";

 static char predual_desc[] = "whether gurobi's presolve should form the dual of a\n\
				continuous model:\n\
		       -1 = automatic choice (default)\n\
			0 = no\n\
			1 = yes\n\
			2 = form both primal and dual and use two threads to\n\
			     choose heuristically between them";
#endif /*}*/

 static char presolve_desc[] = "whether to use Gurobi's presolve:\n\
		       -1 (default) = automatic choice\n\
			0 = no\n\
			1 = conservative presolve\n\
			2 = aggressive presolve";

 static char pricing_desc[] = "pricing strategy:\n\
		       -1 = automatic choice (default)\n\
			0 = partial pricing\n\
			1 = steepest edge\n\
			2 = Devex\n\
			3 = quick-start steepest edge";

#ifdef GRB_INT_ATTR_BRANCHPRIORITY
 static char priorities_desc[] =
		"Whether to use the variable.priority suffix with MIP problems.\n\
		When several branching candidates are available, a variable\n\
		with the highest .priority is chosen for the next branch.\n\
		Priorities are nonnegative integers (default 0).\n\
		Possible values for \"priorities\":\n\
			0 = ignore .priority; assume priority 0 for all vars\n\
			1 = use .priority if present (default).";
#endif

#if GRB_VERSION_MAJOR >= 4 /*{*/
 static char psdtol_desc[] = "maximum diagonal perturbation to correct indefiniteness\n\
		in quadratic objectives (default 1e-6)";
#endif /*}*/

#if GRB_VERSION_MAJOR >= 3 /*{*/
 static char pumppasses_desc[] = "number of feasibility-pump passes to do after the\n\
			MIP root when no other root heuristoc found a\n\
			feasible solution (default 0)";
#endif /*}*/

#ifdef GRB_INT_PAR_QCPDUAL
 static char qcpdual_desc[] = "Whether to compute dual variables when the problem\n\
		has quadratic constraints (which can be expensive):\n\
			0 = no (default)\n\
			1 = yes";
#endif

#if (GRB_VERSION_MAJOR == 4 && GRB_VERSION_MINOR >= 5) || GRB_VERSION_MAJOR >= 5 /*{*/
 static char rays_desc[] = "Whether to return suffix .unbdd if the objective is unbounded\n\
		or suffix .dunbdd if the constraints are infeasible:\n\
			0 = neither\n\
			1 = just .unbdd\n\
			2 = just .dunbdd\n\
			3 = both (default)";
#endif /*}*/

 static char relax_desc[] = "whether to relax integrality:\n\
			0 = no (default)\n\
			1 = yes: treat integer and binary variables\n\
				as continuous";

 static char return_mipgap_desc[] =
		"Whether to return mipgap suffixes or include mipgap values\n\
		(|objectve - best_bound|) in the solve_message:  sum of\n\
			1 = return relmipgap suffix (relative to |obj|);\n\
			2 = return absmipgap suffix (absolute mipgap);\n\
			4 = suppress mipgap values in solve_message.\n\
		Default = 0.  The suffixes are on the objective and problem.\n\
		Returned suffix values are +Infinity if no integer-feasible\n\
		solution has been found, in which case no mipgap values are\n\
		reported in the solve_message.";

#ifdef GRB_INT_PAR_SEED
 static char seed_desc[] =
		"random number seed (default 0), affecting perturbations that\n\
		may influence the solution path.";
#endif

#ifdef GRB_INT_PAR_SIFTING /*{ new in 4.6 */
 static char sifting_desc[] =
		"whether to use sifting within the dual simplex algorithm,\n\
		which can be useful when there are many more variables than\n\
		constraints:\n\
			-1 = automatic choice (default)\n\
			 0 = no\n\
			 1 = yes, moderate sifting\n\
			 2 = yes, aggressive sifting.";

  static char siftmethod_desc[] =
		"Algorithm to use for sifting with the dual simplex method:\n\
			-1 = automatic choice (default)\n\
			 0 = primal simplex\n\
			 1 = dual simplex\n\
			 2 = barrier.";
#endif /*}*/


 static char scale_desc[] = "whether to scale the problem:\n\
			0 = no\n\
			1 = yes (default)";

 static char simplex_desc[] =
#if GRB_VERSION_MAJOR < 4
			"which algorithm to use:"
#else
			"which algorithm to use for non-MIP problems or for the root\n\
		node of MIP problems:"
#endif
#if (GRB_VERSION_MAJOR == 4 && GRB_VERSION_MINOR >= 5) || GRB_VERSION_MAJOR >= 5 /*{*/
			"\n\
			-1 automatic (default): 3 for LP, 2 for QP, 1 for MIP\n\
				root node"
#endif /*}*/
			"\n\
			0 = primal simplex\n\
			1 = dual simplex (default)"
#if GRB_VERSION_MAJOR >= 3 /*{*/
			"\n\
			2 = barrier"
#if (GRB_VERSION_MAJOR == 4 && GRB_VERSION_MINOR >= 5) || GRB_VERSION_MAJOR >= 5 /*{*/
			"\n\
			3 = nondeterministic concurrent (several solves in\n\
				parallel)\n\
			4 = deterministic concurrent"
#endif
#endif /*}*/
			;

 static char solnsens_desc[] = "whether to return suffixes for solution sensitivities, i.e.,\n\
		ranges of values for which the optimal basis remains optimal:\n\
			0 = no (default)\n\
			1 = yes:  suffixes return on variables are\n\
				.sensobjlo = smallest objective coefficient\n\
				.sensobjhi = greatest objective coefficient\n\
				.senslblo = smallest variable lower bound\n\
				.senslbhi = greatest variable lower bound\n\
				.sensublo = smallest variable upper bound\n\
				.sensubhi = greatest variable upper bound\n\
			suffixes for constraints are\n\
				.sensrhslo = smallest right-hand side value\n\
				.sensrhshi = greatest right-hand side value"
#if GRB_VERSION_MAJOR >= 5
			"\n\
		For problems with integer variables and quadratic constraints,\n\
		solnsens = 0 is assumed quietly."
#endif
;

 static char sos_desc[] = "whether to honor declared suffixes .sosno and .ref describing\n\
		SOS sets:\n\
			0 = no\n\
			1 = yes (default):  each distinct nonzero .sosno\n\
				value designates an SOS set, of type 1 for\n\
				positive .sosno values and of type 2 for\n\
				negative values.  The .ref suffix contains\n\
				corresponding reference values.";

 static char sos2_desc[] = "whether to tell Gurobi about SOS2 constraints for nonconvex\n\
		piecewise-linear terms:\n\
			0 = no\n\
			1 = yes (default), using suffixes .sos and .sosref\n\
				provided by AMPL.";

 static char threads_desc[] = "maximum threads to use on MIP problems\n\
		(default 0 ==> max possible)";

 static char timing_desc[] = "whether to report timing:\n\
			0 (default) = no\n\
			1 = report times on stdout\n\
			2 = report times on stderr";

#ifdef GRB_INT_PAR_TUNEOUTPUT /*{*/

 static char tunebase_desc[] = "base name for results of running Gurobi's search for better\n\
		parameter settings.  The search is run only when tuneparbase is\n\
		specified.  Results are written to files with names derived\n\
		from tunebase by appending \".prm\" if \".prm\" does not occur in\n\
		tuneparbase and inserting 1, 2, ... (for the first, second, ...\n\
		set of parameter settings) before the right-most \".prm\".\n\
		The file with \"1\" inserted is the best set and the solve\n\
		results returned are for this set.  In a subsequent \"solve;\",\n\
		you can use paramfile=... to apply the settings in results\n\
		file ... .";

 static char tuneoutput_desc[] = "amount of tuning output when tunebase is specified:\n\
			0 = none\n\
			1 = summarize each new best parameter set\n\
			2 = summarize each set tried (default)\n\
			3 = summary plus detailed solver output for each trial";

 static char tuneresults_desc[] = "limit on the number of tuning result files to write\n\
		when tunerbase is specified.  The default (-1) is to write\n\
		results for all parameter sets on the efficient frontier.";

 static char tunetimelim_desc[] = "time limit (in seconds) on tuning when tunebase\n\
		is specified.  Default -1 ==> automatic choice of time limit.";

 static char tunetrials_desc[] = "number of trials for each parameter set when tunebase\n\
		is specified, each with a different random seed value.\n\
		Default = 2.";

#endif /*}*/

 static char varbranch_desc[] = "MIP branch variable selection strategy:\n\
		       -1 = automatic choice (default)\n\
			0 = pseudo reduced-cost branching\n\
			1 = pseudo shadow-price branching\n\
			2 = maximum infeasibility branching\n\
			3 = strong branching";

 static char version_desc[] = "Report version details before solving the problem.  This is a\n\
		single-word \"phrase\" that does not accept a value assignment.";

 static char writeprob_desc[] = "name of a GUROBI-format file to be written (for debugging);\n\
		must end in one of \".bas\", \".lp\", \".mps\", \".prm\", \".sol\", or\n\
		for the \"fixed\" model used to recover a basis or dual values\n\
		for problems with integer variables or quadratic constraints,\n\
		\".fix_lp\" or \".fix_mps\"; the '_' will be replaced by '.' in\n\
		the name of the file written for \".fix_lp\" or \".fix_mps\".\n\
		Can appear more than once (with different filenames).";

 /* WS_desc_ASL = modified solvers/ws_desc.c: extra initial tab; must not be static. */
 char WS_desc_ASL[] = "=... solution report without -AMPL: sum of\n\
			1 ==> write .sol file\n\
			2 ==> print primal variable values\n\
			4 ==> print dual variable values\n\
			8 ==> do not print solution message";

#if GRB_VERSION_MAJOR > 1 /*{*/

 static char multprice_norm_desc[] = "choice of norm used in multiple pricing:\n\
		       -1 = automatic choice (default)\n\
			0, 1, 2, 3 = specific choices:  hard to describe,\n\
				but sometimes a specific choice will perform\n\
				much better than the automatic choice.";

 static char nodefiledir_desc[] = "directory where MIP tree nodes are written after memory\n\
		for them exceeds nodefilestart; default \".\"";

 static char nodefilestart_desc[] = "gigabytes of memory to use for MIP tree nodes;\n\
		default = Infinity (no limit, i.e., no node files written)";

#ifdef GRB_INT_PAR_PREMIQPMETHOD /*{*/
 static char premiqpmethod_desc[] =
#ifdef GRB_INT_PAR_PREQLINEARIZE
			"Deprecated; replaced by preqlinearize.\n\
		Same possible values as preqlinearize.";
#else
			"how Gurobi's presolve should treat MIQP problems:\n\
			-1 = automatic choice (default)\n\
			 0 = leave the problem as an MIQP\n\
			 1 = try to transform an MIQP to an MILP";
#endif
#endif /*}*/

#if GRB_VERSION_MAJOR >= 3 /*{*/
 static char prepasses_desc[] = "limit on the number of Gurobi presolve passes:\n\
		       -1 = automatic choice (automatic)\n\
			n >= 0: at most n passes";
#endif /*}*/

#ifdef GRB_INT_PAR_PREQLINEARIZE /*{*/
 static char preqlinearize_desc[] = "How Gurobi's presolve should treat quadratic problems:\n\
			-1 = automatic choice (default)\n\
			 0 = do not modify the quadratic part(s)\n\
			 1 = try to linearize quadratic parts";
#endif /*}*/

#ifdef GRB_INT_PAR_PRESPARSIFY
 static char presparsify_desc[] =
		"Whether Gurobi's presolve should use its \"sparsify reduction\",\n\
		which sometimes results in significant problem-size reductions:\n\
			0 = no (default)\n\
			1 = yes.";
#endif

 static char quad_desc[] = "whether simplex should use quad-precision:\n\
		       -1 = automatic choice (default)\n\
			0 = no\n\
			1 = yes";

#if GRB_VERSION_MAJOR >= 3 /*{*/
 static char resultfile_desc[] = "name of a file of extra information written after\n\
				completion of optimization.  The name's suffix\n\
				determines what is written:\n\
			.sol	solution vector\n\
			.bas	simplex basis\n\
			.mst	integer variable solution vector";

 static char rins_desc[] = "how often to apply the RINS heuristic for MIP problems:\n\
		       -1 = automatic choice (default)\n\
			0 = never\n\
			n > 0: every n-th node";
#endif /*}*/

#if GRB_VERSION_MAJOR < 4 /*{*/
 static char rootmethod_desc[] = "algorithm for MIP root relaxation:\n\
			0 = primal simplex\n\
			1 = dual simplex (default)"
#if GRB_VERSION_MAJOR >= 3
			"\n\
			2 = barrier"
#endif
			;
#endif /*}*/
#endif /*}*/

#if GRB_VERSION_MAJOR >= 3 /*{*/
 static char symmetry_desc[] = "MIP symmetry detection:\n\
		       -1 = automatic choice (default)\n\
			0 = none\n\
			1 = conservative\n\
			2 = agressive";
#endif /*}*/

#if GRB_VERSION_MAJOR >= 3 /*{*/
 static char warmstart_desc[] = "Whether to use incoming primal and dual variable values\n\
		(if both are available) in a simplex warm start:\n\
			0 = no\n\
			1 = yes if there is no incoming basis (default)\n\
		Note that specifying basis=0 or basis=2 causes there to be\n\
		no incoming basis.";
#endif /*}*/

#ifdef GRB_INT_PAR_ZEROOBJNODES
 static char zeroobjnodes_desc[] =
		"Number of nodes to explore at the root MIP node if no other\n\
		heuristic has found a feasible solution.  Default = 0.";
#endif

#define VP (void*)

#if GRB_VERSION_MAJOR >= 4
#define Method "Method"
#else
#define Method "LPMethod"
#endif

 static keyword
keywds[] = {	/* must be in alphabetical order */

#if GRB_VERSION_MAJOR >= 3
	{ "aggfill", sf_ipar, "AggFill", aggfill_desc },
#endif
	{ "aggregate", sf_ipar, "Aggregate", aggregate_desc },
#if GRB_VERSION_MAJOR >= 3 /*{*/
	{ "ams_eps", D_val, &ams_eps, ams_eps_desc },
	{ "ams_epsabs", D_val, &ams_epsabs, ams_epsabs_desc },
	{ "ams_limit", I_val, &ams_limit, ams_limit_desc },
	{ "ams_stub", C_val, &ams_stub, ams_stub_desc },
	{ "barconvtol", sf_dpar, "BarConvTol", barconvtol_desc },
	{ "barcorrectors", sf_ipar, "BarCorrectors", barcorrectors_desc },
#ifdef GRB_INT_PAR_BARHOMOGENEOUS
	{ "barhomogeneous", sf_ipar, "BarHomogeneous", barhomogeneous_desc },
#endif
	{ "bariterlim",  sf_ipar, "BarIterLimit", bariterlimit_desc },
	{ "barorder", sf_ipar, "BarOrder", barorder_desc },
#endif /*}*/
#ifdef GRB_DBL_PAR_BARQCPCONVTOL
	{ "barqcptol", sf_dpar, "BarQCPConvTol", barqcptol_desc },
#endif
	{ "basis", sf_mint, VP set_basis, basis_desc },
	{ "bestbound", sf_mint, VP set_bestbound, bestbound_desc },
#if (GRB_VERSION_MAJOR == 4 && GRB_VERSION_MINOR >= 5) || GRB_VERSION_MAJOR >= 5 /*{*/
	{ "branchdir", sf_ipar, "BranchDir", branchdir_desc },
#endif /*}*/
	{ "cliquecuts", sf_ipar, "CliqueCuts", overrides_cuts },
#ifdef GRB_INT_PAR_CONCURRENTMIP
	{ "concurrentmip", sf_ipar, GRB_INT_PAR_CONCURRENTMIP, concurrentmip_desc },
#endif
	{ "covercuts", sf_ipar, "CoverCuts", overrides_cuts },
#if GRB_VERSION_MAJOR >= 3 /*{*/
	{ "crossover", sf_ipar, "Crossover", crossover_desc },
	{ "crossoverbasis", sf_ipar, "CrossoverBasis", crossoverbasis_desc },
#endif /*}*/
	{ "cutagg", sf_ipar, "CutAggPasses", cutagg_desc },
	{ "cutoff", sf_dpar, "Cutoff", cutoff_desc },
#if (GRB_VERSION_MAJOR == 4 && GRB_VERSION_MINOR >= 5) || GRB_VERSION_MAJOR >= 5 /*{*/
	{ "cutpasses", sf_ipar, "CutPasses", cutpasses_desc },
#endif /*}*/
	{ "cuts", sf_ipar, "Cuts", cuts_desc },
#ifdef GRB_INT_PAR_DUALREDUCTIONS
	{ "dualreductions", sf_ipar, "DualReductions", dualreductions_desc },
#endif
#if GRB_VERSION_MAJOR >= 5 /*{*/
	{ "feasrelax", sf_mint, VP set_feasrelax, feasrelax_desc },
#ifdef GRB_DBL_PAR_FEASRELAXBIGM
	{ "feasrelaxbigm", sf_dpar, GRB_DBL_PAR_FEASRELAXBIGM, feasrelaxbigm_desc },
#endif
#endif /*}*/
	{ "feastol", sf_dpar, "FeasibilityTol", feastol_desc },
	{ "flowcover", sf_ipar, "FlowCoverCuts", "flowcover cuts:  " Overrides_cuts },
	{ "flowpath", sf_ipar, "FlowPathCuts", "flowpath cuts:  " Overrides_cuts },
	{ "gomory", sf_ipar, "GomoryPasses", gomory_desc },
	{ "gubcover", sf_ipar, "GUBCoverCuts", "gubcover cuts:  " Overrides_cuts },
	{ "heurfrac", sf_dpar, "Heuristics", heurfrac_desc },
	{ "iisfind", sf_mint, VP set_iis, iisfind_desc },
#if GRB_VERSION_MAJOR > 1 /*{*/
	{ "iismethod", sf_ipar, "IISMethod", iismethod_desc },
#endif /*}*/
	{ "implied", sf_ipar, "ImpliedCuts", "implied cuts:  " Overrides_cuts },
#ifdef GRB_DBL_PAR_IMPROVESTARTGAP
	{ "improvegap", sf_dpar, GRB_DBL_PAR_IMPROVESTARTGAP, isg_desc },
#endif
#ifdef GRB_DBL_PAR_IMPROVESTARTTIME
	{ "improvetime", sf_dpar, GRB_DBL_PAR_IMPROVESTARTTIME, ist_desc },
#endif
#ifdef GRB_DBL_PAR_IMPROVESTARTNODES
	{ "impstartnodes", sf_dpar, GRB_DBL_PAR_IMPROVESTARTNODES, isn_desc },
#endif
	{ "intfeastol", sf_dpar, "IntFeasTol", intfeastol_desc },
	{ "intstart", sf_mint, VP set_intstart, intstart_desc },
	{ "iterlim", sf_dpar, "IterationLimit", "iteration limit (default: no limit)" },
#if GRB_VERSION_MAJOR >= 5
	{ "lbpen", D_val, &lbpen, "See feasrelax." },
#endif
	{ "logfile", C_val, &logfile, logfile_desc },
	{ "logfreq", sf_iparlog, "DisplayInterval", logfreq_desc },
	{ "lpmethod", sf_ipar, Method, "synonym for \"method\"" },
	{ "maxmipsub", sf_ipar, "SubMIPNodes", maxmipsub_desc },
	{ "method", sf_ipar, Method, simplex_desc },
#if (GRB_VERSION_MAJOR == 4 && GRB_VERSION_MINOR >= 5) || GRB_VERSION_MAJOR >= 5 /*{*/
	{ "minrelnodes", sf_ipar, "MinRelNodes", minrelnodes_desc },
#endif /*}*/
#if GRB_VERSION_MAJOR >= 3
	{ "mipfocus", sf_ipar, "MIPFocus", mipfocus_desc },
#endif
	{ "mipgap", sf_dpar, "MipGap", "max. relative MIP optimality gap (default 1e-4)" },
#if GRB_VERSION_MAJOR >= 3
	{ "mipgapabs", sf_dpar, "MipGapAbs", "absolute MIP optimality gap (default 1e-10)" },
#endif
	{ "mipsep", sf_ipar, "MIPSepCuts", "MIPsep cuts:  " Overrides_cuts },
	{ "mipstart", sf_mint, VP set_mipstval, mipstart_desc },
	{ "mircuts", sf_ipar, "MIRCuts", "MIR cuts:  " Overrides_cuts },
#if GRB_VERSION_MAJOR >= 4
	{ "modkcuts", sf_ipar, "ModKCuts", "mod-k cuts:  " Overrides_cuts },
#endif
#if GRB_VERSION_MAJOR > 1 /*{*/
	{"multprice_norm", sf_ipar, "NormAdjust", multprice_norm_desc},
#if GRB_VERSION_MAJOR >= 3
	{ "networkcuts", sf_ipar, "NetworkCuts", "Network cuts:  " Overrides_cuts },
#endif
	{"nodefiledir", sf_spar, "NodefileDir", nodefiledir_desc},
	{"nodefilestart", sf_dpar, "NodefileStart", nodefilestart_desc},
#endif /*}*/
	{ "nodelim", sf_dpar, "NodeLimit", "maximum MIP nodes to explore (default: no limit)" },
#if GRB_VERSION_MAJOR >= 4 /*{*/
	{ "nodemethod", sf_ipar, "NodeMethod", nodemethod_desc },
#endif /*}*/
#if GRB_VERSION_MAJOR > 1 /*{*/
	{ "normadjust", sf_ipar, "NormAdjust", "synonym for multprice_norm" },
#endif /*}*/
#ifdef GRB_INT_PAR_NUMERICFOCUS
	{ "numericfocus", sf_ipar, GRB_INT_PAR_NUMERICFOCUS, numericfocus_desc },
#endif
	{ "objno", sf_mint, VP set_objno, objno_desc },
#if GRB_VERSION_MAJOR > 1 /*{*/
	{ "objscale", sf_dpar, "ObjScale", objscale_desc },
#endif /*}*/
	{ "opttol", sf_dpar, "OptimalityTol", opttol_desc },
	{ "outlev", sf_mint, VP set_outlev, outlev_desc },
#if GRB_VERSION_MAJOR >= 3
	{ "param", sf_spar, "Dummy", param_desc },
	{ "paramfile", sf_pf, 0, paramfile_desc },
#endif
	{ "perturb", sf_dpar, "PerturbValue", perturb_desc },
	{ "pivtol", sf_dpar, "MarkowitzTol", "Markowitz pivot tolerance (default 7.8125e-3)" },
	/*GRB_MAJ2(({"precrush", sf_ipar, "PreCrush", precrush_desc},))*/
#if GRB_VERSION_MAJOR >= 3
	{ "predeprow", sf_ipar, "PreDepRow", predeprow_desc },
	{ "predual", sf_ipar, "PreDual", predual_desc },
#ifdef GRB_INT_PAR_PREMIQPMETHOD
	{ "premiqpmethod", sf_ipar, "PreMIQPMethod", premiqpmethod_desc },
#endif
	{ "prepases", sf_ipar, "PrePasses", "deprecated synonym for prepasses" },
	{ "prepasses", sf_ipar, "PrePasses", prepasses_desc },
#endif
#ifdef GRB_INT_PAR_PREQLINEARIZE
	{ "preqlinearize", sf_ipar, GRB_INT_PAR_PREQLINEARIZE, preqlinearize_desc },
#endif
	{ "presolve", sf_ipar, "Presolve", presolve_desc },
#ifdef GRB_INT_PAR_PRESPARSIFY /* new in 4.6 */
	{ "presparsify", sf_ipar, GRB_INT_PAR_PRESPARSIFY, presparsify_desc },
#endif
	{ "pricing", sf_ipar, "SimplexPricing", pricing_desc },
#ifdef GRB_INT_ATTR_BRANCHPRIORITY
	{ "priorities", sf_mint, VP set_priorities, priorities_desc },
#endif
#if GRB_VERSION_MAJOR >= 4
	{ "psdtol", sf_dpar, "PSDTol", psdtol_desc },
#endif
#if GRB_VERSION_MAJOR >= 3
	{ "pumppasses", sf_ipar, "PumpPasses", pumppasses_desc },
#endif
#ifdef GRB_INT_PAR_QCPDUAL
	{ "qcpdual", sf_ipar, "QCPDual", qcpdual_desc },
#endif
#if GRB_VERSION_MAJOR > 1 /*{*/
	{ "quad", sf_ipar, "Quad", quad_desc},
#endif /*}*/
#if (GRB_VERSION_MAJOR == 4 && GRB_VERSION_MINOR >= 5) || GRB_VERSION_MAJOR >= 5 /*{*/
	{ "rays", sf_mint, VP set_rays, rays_desc },
#endif /*}*/
	{ "relax", sf_mint, VP set_relax, relax_desc },
#if GRB_VERSION_MAJOR >= 3
	{ "resultfile", sf_spar, "ResultFile", resultfile_desc },
#endif
	{ "return_mipgap", sf_mint, VP set_retmipgap, return_mipgap_desc },
#if GRB_VERSION_MAJOR >= 5
	{ "rhspen", D_val, &rhspen, "See feasrelax." },
#endif
#if GRB_VERSION_MAJOR >= 3
	{ "rins", sf_ipar, "RINS", rins_desc },
#endif
#if GRB_VERSION_MAJOR > 1 && GRB_VERSION_MAJOR < 4 /*{*/
	{"rootmethod", sf_ipar, "RootMethod", rootmethod_desc},
#endif /*}*/
	{ "scale", sf_ipar, "ScaleFlag", scale_desc },
#ifdef GRB_INT_PAR_SEED
	{ "seed", sf_ipar, GRB_INT_PAR_SEED, seed_desc },
#endif
#ifdef ALLOW_GUROBI_SERVER
	{ "server", C_val, &server, server_desc },
	{ "server_password", C_val, &server_passwd, server_passwd_desc },
	{ "server_port", I_val, &server_port, server_port_desc },
	{ "server_priority", I_val, &server_priority, server_priority_desc },
	{ "server_timeout", D_val, &server_timeout, server_timeout_desc },
	{ "serverlic", C_val, &serverlic, serverlic_desc },
#endif
#ifdef GRB_INT_PAR_SIFTING /* new in 4.6 */
	{ "sifting", sf_ipar, GRB_INT_PAR_SIFTING, sifting_desc },
	{ "siftmethod", sf_ipar, GRB_INT_PAR_SIFTMETHOD, siftmethod_desc },
#endif
	{ "simplex", sf_ipar, Method, "synonym for lpmethod" },
	{ "solnlimit", sf_ipar, "SolutionLimit", "maximum MIP solutions to find (default 2e9)" },
	{ "solnsens", sf_mint, VP set_solnsens, solnsens_desc },
	{ "sos", sf_mint, VP set_sos, sos_desc },
	{ "sos2", sf_mint, VP set_sos2, sos2_desc },
#if GRB_VERSION_MAJOR >= 3
	{ "submipcuts", sf_ipar, "SubMIPCuts", "sub-MIP cuts:  " Overrides_cuts },
	{ "symmetry", sf_ipar, "Symmetry", symmetry_desc },
#endif
	{ "threads", sf_ipar, "Threads", threads_desc },
	{ "timelim", sf_dpar, "TimeLimit", "limit on solve time (in seconds; default: no limit)" },
	{ "timing", sf_mint, VP set_timing, timing_desc },
#ifdef GRB_INT_PAR_TUNEOUTPUT
	{ "tunebase", C_val, &tunebase, tunebase_desc },
	{ "tuneoutput", sf_ipar, GRB_INT_PAR_TUNEOUTPUT, tuneoutput_desc },
	{ "tuneresults", sf_ipar, GRB_INT_PAR_TUNERESULTS, tuneresults_desc },
	{ "tunetimelimit", sf_dpar, GRB_DBL_PAR_TUNETIMELIMIT, tunetimelim_desc },
	{ "tunetrials", sf_ipar, GRB_INT_PAR_TUNETRIALS, tunetrials_desc },
#endif
#if GRB_VERSION_MAJOR >= 5
	{ "ubpen", D_val, &ubpen, "See feasrelax." },
#endif
	{ "varbranch", sf_ipar, "VarBranch", varbranch_desc },
	{ "version", Ver_val, 0, version_desc },
	{ "wantsol", WS_val, 0, WS_desc_ASL+5 },
#if GRB_VERSION_MAJOR >= 5
	{ "warmstart", sf_mint, VP set_warmstart, warmstart_desc },
#endif
	{ "writeprob", sf_wfile, 0, writeprob_desc }
#if GRB_VERSION_MAJOR > 1 /*{*/
	,{"zerohalfcuts", sf_ipar, "ZeroHalfCuts", "zero-half cuts:  " Overrides_cuts }
#endif /*}*/
#ifdef GRB_INT_PAR_ZEROOBJNODES /* new in 4.6 */
	,{"zeroobjnodes", sf_ipar, GRB_INT_PAR_ZEROOBJNODES, zeroobjnodes_desc }
#endif
	};

 static Option_Info
Oinfo = { "gurobi", verbuf, "gurobi_options", keywds, nkeywds, 0, verbuf,
	   0,0,0,0,0, 20130419 };

 static void
enamefailed(GRBenv *env, const char *what, const char *name)
{
	fprintf(Stderr, "%s(\"%s\") failed:\n\t%s.\n", what, name, GRBgeterrormsg(env));
	++Oinfo.n_badopts;
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
	{ "absmipgap", 0, ASL_Sufkind_obj   | ASL_Sufkind_outonly },
	{ "absmipgap", 0, ASL_Sufkind_prob  | ASL_Sufkind_outonly },
	{ "bestbound", 0, ASL_Sufkind_obj   | ASL_Sufkind_outonly },
	{ "bestbound", 0, ASL_Sufkind_prob  | ASL_Sufkind_outonly },
#if (GRB_VERSION_MAJOR == 4 && GRB_VERSION_MINOR >= 5) || GRB_VERSION_MAJOR >= 5 /*{*/
	{ "dunbdd", 0, ASL_Sufkind_con | ASL_Sufkind_outonly },
#endif /*}*/
	{ "iis", iis_table, ASL_Sufkind_var | ASL_Sufkind_outonly },
	{ "iis", 0, ASL_Sufkind_con | ASL_Sufkind_outonly },
#if GRB_VERSION_MAJOR >= 5
	{ "lbpen", 0, ASL_Sufkind_var | ASL_Sufkind_real },
#endif
#ifdef GRB_INT_ATTR_BRANCHPRIORITY
	{  "priority", 0, ASL_Sufkind_var },
#endif
	{ "ref", 0, ASL_Sufkind_var | ASL_Sufkind_real },
	{ "relmipgap", 0, ASL_Sufkind_obj   | ASL_Sufkind_outonly },
	{ "relmipgap", 0, ASL_Sufkind_prob  | ASL_Sufkind_outonly },
#if GRB_VERSION_MAJOR >= 5
	{ "rhspen", 0, ASL_Sufkind_con | ASL_Sufkind_real },
#endif
	{ "senslbhi",  0, ASL_Sufkind_var | ASL_Sufkind_real | ASL_Sufkind_outonly },
	{ "senslblo",  0, ASL_Sufkind_var | ASL_Sufkind_real | ASL_Sufkind_outonly },
	{ "sensobjhi", 0, ASL_Sufkind_var | ASL_Sufkind_real | ASL_Sufkind_outonly },
	{ "sensobjlo", 0, ASL_Sufkind_var | ASL_Sufkind_real | ASL_Sufkind_outonly },
	{ "sensrhshi", 0, ASL_Sufkind_con | ASL_Sufkind_real | ASL_Sufkind_outonly },
	{ "sensrhslo", 0, ASL_Sufkind_con | ASL_Sufkind_real | ASL_Sufkind_outonly },
	{ "sensubhi",  0, ASL_Sufkind_var | ASL_Sufkind_real | ASL_Sufkind_outonly },
	{ "sensublo",  0, ASL_Sufkind_var | ASL_Sufkind_real | ASL_Sufkind_outonly },
	{ "sos", 0, ASL_Sufkind_var },
	{ "sos", 0, ASL_Sufkind_con },
	{ "sosno", 0, ASL_Sufkind_var | ASL_Sufkind_real },
	{ "sosref", 0, ASL_Sufkind_var | ASL_Sufkind_real },
	{ "sstatus", 0, ASL_Sufkind_var, 1 },
	{ "sstatus", 0, ASL_Sufkind_con, 1 }
#if GRB_VERSION_MAJOR >= 5
	,{ "ubpen", 0, ASL_Sufkind_var | ASL_Sufkind_real }
#endif
#if (GRB_VERSION_MAJOR == 4 && GRB_VERSION_MINOR >= 5) || GRB_VERSION_MAJOR >= 5 /*{*/
	,{ "unbdd", 0, ASL_Sufkind_var | ASL_Sufkind_outonly }
#endif /*}*/
	};

#ifdef GRB_INT_ATTR_BRANCHPRIORITY
 static void
add_priorities(ASL *asl, GRBenv *env, GRBmodel *mdl)
{
	SufDesc *dp;
	int *p;

	if ((dp = suf_get("priority", ASL_Sufkind_var))
	 && (p = dp->u.i)
	 && GRBsetintattrarray(mdl, GRB_INT_ATTR_BRANCHPRIORITY, 0, n_var, p))
		failed(env, "GRBsetintattrarray(\"BranchPriority\")");
	}
#endif

 static void
show_times(void)
{
	FILE *f;
	int i;

	Times[3] = xectim_();
	Times[4] = time(0) - Times[4];
	for(i = 1; i <= 2; i++)
	    if (time_flag & i) {
		f = i == 1 ? stdout : Stderr;
		fprintf(f, "\nTimes (seconds):\nInput =  %g"
			"\nSolve =  %g (summed over threads)"
			"\nOutput = %g\nElapsed ",
			Times[1] - Times[0], Times[2] - Times[1],
			Times[3] - Times[2]);
		fprintf(f, Times[4] < 1. ? "< 1\n" : "= %g\n", Times[4]);
		}
	}

 Sig_ret_type
intcatch(int n)
{
	printf("\n<BREAK> (gurobi)\n", n);
	fflush(stdout);
	if (++breaking > 3)
		longjmp(Jb, 2);
	signal(SIGINT, intcatch);
	if (grbmodel)
		GRBterminate(grbmodel);
	SigRet;
	}

 static const char*
retiis(ASL *asl, GRBenv *env, GRBmodel *mdl, Dims *d, const char *what, int *srp)
{
	char buf[128], *rv;
	int *c, i, j, k, kv, m, n, nr, *s, *v;

	m = n_con;
	n = n_var;
	nr = n + nranges;
	if (GRBgetintattrarray(mdl, "IISConstr", 0, m, s = d->rstat))
		failed(env, "GRBgetintattrarray(\"IISConstr\")");
	c = v = 0;
	for(i = k = kv = 0; i < m; ++i) {
		if (s[i]) {
			c = (int*)M1zapalloc(m*sizeof(int));
			for(; i < m; ++i)
				if (s[i]) {
					c[i] = 4;
					++k;
					}
			break;
			}
		}
	if (GRBgetintattrarray(mdl, "IISLB", 0, nr, s = d->cstat))
		failed(env, "GRBgetintattrarray(\"IISLB\")");
	for(i = 0; i < n; ++i) {
		if (s[i]) {
			v = (int*)M1zapalloc(n*sizeof(int));
			for(; i < n; ++i)
				if (s[i]) {
					v[i] = 1;
					++kv;
					}
			break;
			}
		}
	for(i = n; i < nr; ++i) {
		if (s[i]) {
			if (!c)
				c = (int*)M1zapalloc(m*sizeof(int));
			for(; i < nr; ++i)
				if (s[i] && !c[i]) {
					c[i] = 4;
					++k;
					}
			break;
			}
		}
	if (GRBgetintattrarray(mdl, "IISUB", 0, nr, s = d->cstat))
		failed(env, "GRBgetintattrarray(\"IISUB\")");
	for(i = 0; i < n; ++i) {
		if (s[i]) {
			if (!v)
				v = (int*)M1zapalloc(n*sizeof(int));
			for(; i < n; ++i)
				if (s[i] && !v[i]) {
					v[i] = 3;
					++kv;
					}
			break;
			}
		}
	for(i = n; i < nr; ++i) {
		if (s[i]) {
			if (!c)
				c = (int*)M1zapalloc(m*sizeof(int));
			for(; i < nr; ++i)
				if (s[i] && !c[i]) {
					c[i] = 4;
					++k;
					}
			break;
			}
		}
	if (c)
		suf_iput("iis", ASL_Sufkind_con, c);
	if (v)
		suf_iput("iis", ASL_Sufkind_var, v);
	*srp = 201;
	if (k) {
		j = Snprintf(buf, sizeof(buf), "%s\nReturning an IIS of %d constraints",
			what, k);
		j += Snprintf(buf+j, sizeof(buf)-j, kv ? " and %d variables." : ".", kv);
		}
	else if (kv)
		j = Snprintf(buf, sizeof(buf), "%s\nReturning an IIS of %d variables.",
			what, kv);
	else {
		j = Snprintf(buf, sizeof(buf), "%s; empty IIS!");
		*srp = 202;
		}
	rv = (char*)M1alloc(++j);
	memcpy(rv, buf, j);
	return rv;
	}

 static void
dpf(Dims  *d, const char *fmt, ...)
{
	size_t L;
	va_list ap;

	if ((L = d->mbend - d->mb) > 0) {
		va_start(ap, fmt);
		d->mb += Vsnprintf(d->mb, L, fmt, ap);
		va_end(ap);
		}
	}

 static void
missing_msg(Dims *d)
{
	static const char *missing[3] = {
		"primal",
		"dual",
		"primal or dual"
		};
	dpf(d, "\nNo %s variables returned.", missing[d->missing - 1]);
	}

#if (GRB_VERSION_MAJOR == 4 && GRB_VERSION_MINOR >= 5) || GRB_VERSION_MAJOR >= 5 /*{*/
#ifdef RAYDEBUG
#define Debug(x) x
#else
#define Debug(x) /*nothing*/
#endif

 static int
send_ray(ASL *asl, Dims *d, GRBenv *env, GRBmodel *mdl)
{
	int i, n, n0;
	real *c, t, *y;

	n = n_var;
	n0 = d->nv0;
	y = (real *)M1zapalloc(n0 * 2*sizeof(real));
	c = y + n;
	if (n > n0)
		n = n0;
	if (GRBgetdblattrarray(mdl, GRB_DBL_ATTR_UNBDRAY, 0, n, y)) {
		Debug(printf("Get UnbdRay failed: %s\n", GRBgeterrormsg(env)));
		return 0;
		}
	if (!GRBgetdblattrarray(mdl, GRB_DBL_ATTR_OBJ, 0, n, c)) {
		t = 0.;
		for(i = 0; i < n; ++i)
			t += c[i]*y[i];
		if (d->objsense * t > 0.)
			for(i = 0; i < n; ++i)
				y[i] = -y[i];
		}
	suf_rput("unbdd", ASL_Sufkind_var, y);
	return 1;
	}

 static int
send_dray(ASL *asl, Dims *d, GRBenv *env, GRBmodel *mdl)
{
	char *sense;
	int i, n, n0;
	real *rhs, t, t1, *y;

	n = n_con;
	n0 = d->nc0;
	y = (real *)M1zapalloc(n*sizeof(real) + n0*(sizeof(real) + 1));
	rhs = y + n;
	sense = (char*)(rhs + n0);
	if (n > n0)
		n = n0;
	if (GRBgetdblattrarray(mdl, GRB_DBL_ATTR_FARKASDUAL, 0, n, y)) {
		Debug(printf("Get FarkasDual failed: %s\n", GRBgeterrormsg(env)));
		return 0;
		}
	if (GRBgetdblattrarray(mdl, GRB_DBL_ATTR_RHS, 0, n, rhs))
		Debug(printf("Get RHS failed: %s\n", GRBgeterrormsg(env)));
	else if (GRBgetcharattrarray(mdl, GRB_CHAR_ATTR_SENSE, 0, n, sense))
		Debug(printf("Get SENSE failed: %s\n", GRBgeterrormsg(env)));
	else {
		t = 0.;
		for(i = 0; i < n; ++i) {
			t1 = y[i]*rhs[i];
			if (sense[i] == '>')
				t1 = -t1;
			t += t1;
			}
		if (d->objsense * t < 0.)
			for (i = 0; i < n; ++i)
				y[i] = -y[i];
		}
	suf_rput("dunbdd", ASL_Sufkind_con, y);
	return 1;
	}

#undef Debug
#endif /*}*/

#if GRB_VERSION_MAJOR >= 5 /*{*/
 static void
pen_set(real *x, int n, SufDesc *d, real pen)
{
	real *s, t, *xe;

	if (pen <= 0.)
		pen = GRB_INFINITY;
	xe = x + n;
	if (d && (s = d->u.r)) {
		while(x < xe)
			*x++ = (t = *s++) > 0. ? t : pen;
		return;
		}
	while(x < xe)
		*x++ = pen;
	}

 static real *
do_feasrelax(ASL *asl, GRBenv *env, GRBmodel *mdl, const char **objqual, real *fto)
{
	SufDesc *dlb, *drhs, *dub;
	int mr, nc, nlb, nr, nrhs, nub, nv, nvr, rt;
	real pen, t;
	real *lbp, *lbr, *lu, *lue, *p, *rhp, *rv, *ubp, *ubr, *x;
	size_t L;

	nlb = lbpen > 0.;
	nub = ubpen > 0.;
	nrhs = rhspen > 0.;
	if ((dlb = suf_get("lbpen", ASL_Sufkind_var | ASL_Sufkind_input)))
		++nlb;
	if ((drhs = suf_get("rhspen", ASL_Sufkind_con | ASL_Sufkind_input)))
		++nrhs;
	if ((dub = suf_get("ubpen", ASL_Sufkind_var | ASL_Sufkind_input)))
		++nub;
	if (!(nlb + nub + nrhs))
		return 0;
	nc = n_con;
	nv = n_var;
	if ((nr = nranges) && nrhs) {
		++nlb;
		++nub;
		}
	nvr = nv + nr;
	L = 0;
	if (nlb)
		L += nvr;
	if (nub)
		L += nvr;
	if (nrhs)
		L += nc;
	lbp = rhp = ubp = 0;
	x = (real*)Malloc(L*sizeof(real));
	if (nlb) {
		pen_set(lbp = x, nv, dlb, lbpen);
		x += nvr;
		}
	if (nub) {
		pen_set(ubp = x, nv, dub, ubpen);
		x += nvr;
		}
	if (nrhs) {
		pen_set(rhp = x, nc, drhs, rhspen);
		if (nr) {
			lbr = lbp + nv;
			ubr = ubp + nv;
			if ((pen = rhspen) <= 0.)
				pen = Infinity;
			lu = LUrhs;
			lue = lu + 2*nc;
			if (drhs && (p = drhs->u.r)) do {
				t = *p++;
				if (lu[0] > negInfinity
				 && lu[1] < Infinity
				 && lu[0] < lu[1]) {
					if (t <= 0.)
						t = pen;
					*lbr++ = *ubr++ = t;
					}
				} while((lu += 2) < lue);
			else do
				*lbr++ = *ubr++ = pen;
				while(--nr > 0);
			}
		}
	if ((rt = feasrelax - 1) >= 3) {
		rt -= 3;
		mr = 1;
		rv = fto;
		}
	else {
		mr = 0;
		rv = 0;
		*objqual = "feastol ";
		}
	if (GRBfeasrelax(mdl, rt, mr, lbp, ubp, rhp, fto))
		failed(env, "GRBfeasrelax");
	return rv;
	}
#endif /*}*/

 static const char*
statmsg(ASL *asl, GRBenv *env, GRBmodel *mdl, int i, Dims *d, int *wantobj)
{
	char buf[64], *rv1;
	const char *rv;
	int m, n, nc, nv, nvr, objwant, sr, srd, srp;
	real *x, *y;
	size_t L;

	*wantobj = 0;
	if (i) {
		solve_result_num = 502;
		d->x = d->y = 0;
		switch(i) {
		  case GRB_ERROR_OUT_OF_MEMORY:
			rv = "ran out of memory";
			break;
		  case GRB_ERROR_NO_LICENSE:
			rv = "invalid license";
			break;
		  case GRB_ERROR_SIZE_LIMIT_EXCEEDED:
			rv = "problem size limit exceeded";
			break;
		  case GRB_ERROR_IIS_NOT_INFEASIBLE:
			rv = "bug: IIS problem is infeasible";
			break;
#if GRB_VERSION_MAJOR >= 4 /*{*/
		  case GRB_ERROR_Q_NOT_PSD:
			rv = "quadratic objective is not positive definite";
#if GRB_VERSION_MAJOR >= 5
			if (nlc)
				rv = objno > 0 && objno <= nlo
					? "quadratic objective or constraint is not positive definite"
					: "quadratic constraint is not positive definite";
#endif
			solve_result_num = 524;
			break;
#endif /*}*/
		  default:
			Snprintf(buf, sizeof(buf), "surprise return %d from GRBoptimize", i);
			rv = rv1 = M1alloc(strlen(buf)+1);
			strcpy(rv1, buf);
		  }
		return rv;
		}
	if (GRBgetintattr(mdl, GRB_INT_ATTR_STATUS, &i))
		failed(env, "GRBgetintattr(STATUS)");
	nc = n_con;
	nv = n_var;
	nvr = nv + nranges;
	objwant = 1;
	sr = 0;
	switch(i) {
	  case GRB_OPTIMAL:
		rv = "optimal solution";
		break;
	  case GRB_INFEASIBLE:
		objwant = 0;
		nc = srd = 0;
#if (GRB_VERSION_MAJOR == 4 && GRB_VERSION_MINOR >= 5) || GRB_VERSION_MAJOR >= 5 /*{*/
		if (rays & 2)
			srd = send_dray(asl, d, env, mdl);
#endif /*}*/
		if (!want_iis) {
			sr = 200;
			rv = "infeasible";
			}
		else if (GRBcomputeIIS(mdl)) {
			sr = 202;
			rv = "infeasible; no IIS found";
			}
		else
			rv = retiis(asl, env, mdl, d, "infeasible", &sr);
#if (GRB_VERSION_MAJOR == 4 && GRB_VERSION_MINOR >= 5) || GRB_VERSION_MAJOR >= 5 /*{*/
		if (srd) {
 have_srd:
			sr += 3;
			rv1 = (char*)M1alloc(L = strlen(rv) + 40);
			snprintf(rv1, L, "%s; constraint.dunbdd returned.", rv);
			rv = rv1;
			}
#endif /*}*/
		break;
	  case GRB_INF_OR_UNBD:
		objwant = 0;
		nc = nv = 0;
#if GRB_VERSION_MAJOR >= 4 && GRB_VERSION_MINOR >= 5 /*{*/
		srd = srp = 0;
		if (rays & 1)
			srp = send_ray(asl, d, env, mdl);
		if (rays & 2 && srp != 1)
			srd = send_dray(asl, d, env, mdl);
#endif /*}*/
		if (!want_iis) {
			sr = 300;
			rv = "infeasible or unbounded";
			}
		else if (GRBcomputeIIS(mdl)) {
			sr = 301;
			rv = "infeasible or unbounded; no IIS";
			srp = srd = 0;
			}
		else
			rv = retiis(asl, env, mdl, d, "infeasible or unbounded", &sr);
#if (GRB_VERSION_MAJOR == 4 && GRB_VERSION_MINOR >= 5) || GRB_VERSION_MAJOR >= 5 /*{*/
		if (srd)
			goto have_srd;
		if (srp)
			goto have_srp;
#endif /*}*/
		break;
	  case GRB_UNBOUNDED:
		rv = "unbounded";
		objwant = 0;
		nc = nv = srp = 0;
		sr = 300;
#if (GRB_VERSION_MAJOR == 4 && GRB_VERSION_MINOR >= 5) || GRB_VERSION_MAJOR >= 5 /*{*/
		if (rays & 1) {
			srp = send_ray(asl, d, env, mdl);
 have_srp:
			sr += 2;
			rv1 = (char*)M1alloc(L = strlen(rv) + 30);
			snprintf(rv1, L, "%s; variable.unbdd returned.", rv);
			rv = rv1;
			}
#endif /*}*/
		break;
	  case GRB_CUTOFF:
		rv = "objective cutoff";
		sr = 400;
		break;
	  case GRB_ITERATION_LIMIT:
		rv = "iteration limit";
		sr = 401;
		break;
	  case GRB_NODE_LIMIT:
		rv = "node limit";
		sr = 402;
		break;
	  case GRB_TIME_LIMIT:
		rv = "time limit";
		sr = 403;
		break;
	  case GRB_SOLUTION_LIMIT:
		rv = "solution limit";
		sr = 404;
		break;
	  case GRB_INTERRUPTED:
		rv = "interrupted";
		sr = 600;
		break;
	  case GRB_NUMERIC:
		rv = "numeric error";
		sr = 520;
		break;
#ifdef GRB_SUBOPTIMAL
	  case GRB_SUBOPTIMAL:
		rv = "suboptimal";
		sr = 100;
		break;
#endif
	  default:
		Snprintf(buf, sizeof(buf), "surprise status %d after GRBoptimize", i);
		rv = rv1 = (char*)M1alloc(strlen(buf)+1);
		sr = 530;
		strcpy(rv1, buf);
	  }
	solve_result_num = sr;
	x = y = 0;
	if ((n = nvr + nc) > 0) {
		x = (real*)M1alloc(n*sizeof(real));
		d->y0 = y = x + nvr;
		}
	m = 0;
	if (!nv)
		x = 0;
	else if (GRBgetdblattrarray(mdl, GRB_DBL_ATTR_X, 0, nvr, x)) {
		x = 0;
		m = 1;
		}
	if (!nc)
		y = 0;
#if GRB_VERSION_MAJOR >= 5
	else if ((nlc && GRBgetdblattrarray(mdl, GRB_DBL_ATTR_QCPI, 0, nlc, y))
		|| (nc > nlc && GRBgetdblattrarray(mdl, GRB_DBL_ATTR_PI, 0, nc-nlc, y+nlc)))
#else
	else if (GRBgetdblattrarray(mdl, GRB_DBL_ATTR_PI, 0, nc, y))
#endif
		{
		y = 0;
		m += 2;
		}
	d->missing = m;
	d->x = x;
	d->y = y;
	if (!x)
		objwant = 0;
	*wantobj = objwant;
	return rv;
	}

 static void
stat_map(int *stat, int n, int *map, int mx, char *what)
{
	int bad, i, i1, j, j1;
	static char badfmt[] = "gurobi driver: %s[%d] = %d\n";

	bad = i1 = j1 = 0;
	for(i = 0; i < n; i++) {
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
		"gurobi driver: %d messages about bad %s values suppressed.\n",
				bad-1, what);
		}
	}

 static int
get_input_statuses(ASL *asl, GRBenv *env, GRBmodel *mdl, Dims *d)
{
	SufDesc *sd;
	int i, m, n, nvr, *rs, *rsta;
	real *lu;
	static int vmap[] = {-3, 0, -3, -1, -2, -3, -3};
	static int cmap[] = {-1, 0, -1, -1, -1, -1, -1};

	sd = d->csd;
	n = n_var;
	if (!(sd->kind & ASL_Sufkind_input))
		return 0;
	sd = d->rsd;
	m = n_con;
	if (!(sd->kind & ASL_Sufkind_input))
		return 0;
	stat_map(d->cstat, n, vmap, 6, "incoming cstat");
	stat_map(d->rstat, m, cmap, 6, "incoming rstat");
	nvr = n + nranges;
	if (nvr > n) {
		rs = d->rstat;
		rsta = d->cstat + n;
		lu = LUrhs;
		for(i = 0; i < m; ++i, lu += 2) {
			if (lu[0] > negInfinity && lu[0] < lu[1] && lu[1] < Infinity) {
				if (rs[i] == 0)
					*rsta = -1;
				++rsta;
				}
			}
		}
	if (GRBsetintattrarray(mdl, GRB_INT_ATTR_VBASIS, 0, nvr, d->cstat))
		failed(env, "GRBsetintattrarray(\"VBasis\")");
	if (GRBsetintattrarray(mdl, GRB_INT_ATTR_CBASIS, 0, m, d->rstat))
		failed(env, "GRBsetintattrarray(\"CBasis\")");
	return 1;
	}

 static void
intbasis_fail(Dims *d, const char *call)
{
	dpf(d, "\nintbasis trouble: GRB%s failed.", call);
	}

 static GRBmodel*
fixed_model(ASL *asl, GRBmodel *mdl0, Dims *d)
{
	Filename *fn;
	GRBenv *env;
	GRBmodel *mdl;
	double f, *y;
	int i, k;
#if GRB_VERSION_MAJOR >= 5
	int m, m1, nqc;
#endif
	static char *statusname[] = {
		"infeasible",
		"infeasible or unbounded",
		"unbounded",
		"cutoff",
		"iteration limit",
		"node limit",
		"time limit",
		"solution limit",
		"interrupted",
		"numeric difficulty",
		"suboptimal"
		};

	if (!(mdl = GRBfixedmodel(mdl0)))
		return 0;
	if (!(env = GRBgetenv(mdl))) {
		dpf(d, "\nGRBgetenv failed in fixed_model().");
 badret:
		GRBfreemodel(mdl);
		return 0;
		}
	if (GRBsetintparam(env, "Presolve", 0)) {
		intbasis_fail(d, "setintparam(\"Presolve\")");
		goto badret;
		}
	if (!GRBgetintparam(env, Method, &k) && (k >= 2 || k < 0))
		GRBsetintparam(env, Method, 1);
	if ((fn = Wflist[2])) {
		GRBupdatemodel(mdl);
		do {
			if (GRBwrite(mdl, fn->name))
				enamefailed(env, "GRBwrite", fn->name);
			} while((fn = fn->next));
		}
	if (GRBoptimize(mdl)) {
		intbasis_fail(d, "optimize()");
		goto badret;
		}
	if (GRBgetintattr(mdl, GRB_INT_ATTR_STATUS, &i)) {
		intbasis_fail(d, "getintattr()");
		goto badret;
		}
	if (i != GRB_OPTIMAL) {
		if (i >= GRB_INFEASIBLE && i <= GRB_SUBOPTIMAL)
			dpf(d, "\nGRBoptimize of fixed model: %s.",
				statusname[i-GRB_INFEASIBLE]);
		else
			dpf(d, "\nSurprise status %d after GRBoptimize of fixed model.",
				i);
		goto badret;
		}
	if (d->missing & 2 && (y = d->y0)) {
#if GRB_VERSION_MAJOR < 5
		if (!GRBgetdblattrarray(mdl, GRB_DBL_ATTR_PI, 0, n_con, y))
#else
		m = n_con;
		nqc = nlc;
		k = 0;
		if (nqc > 0)
			k = GRBgetdblattrarray(mdl, GRB_DBL_ATTR_QCPI, 0, nqc, y);
		if ((m1 = m - nqc) > 0 && !k)
			k = GRBgetdblattrarray(mdl, GRB_DBL_ATTR_PI, 0, m1, y+nqc);
		if (!k)
#endif
			{
			d->y = y;
			d->missing &= ~2;
			}
		}
	if (!GRBgetdblattr(mdl, GRB_DBL_ATTR_ITERCOUNT, &f)) {
		if (f > 0.)
			dpf(d, "\nplus %.0f simplex iteration%s for intbasis",
				f, "s" + (f == 1.));
		}
	return mdl;
	}
#undef Method

 static void
get_output_statuses(ASL *asl, GRBmodel *mdl, Dims *d)
{
	int i, j, m, n, nr, rv, *s;
	static int vmap[4] = {2, 4, 3, 1};

	m = n_con;
	n = n_var;
	nr = n + nranges;
	rv = 1;
	if (GRBgetintattrarray(mdl, GRB_INT_ATTR_VBASIS, 0, nr, d->cstat)
	 || GRBgetintattrarray(mdl, GRB_INT_ATTR_CBASIS, 0, m,  d->rstat)) {
		/*failed(env, "GRBgetintattrarray(\"VBasis\")");*/
		d->csd->kind &= ~ASL_Sufkind_output;
		d->rsd->kind &= ~ASL_Sufkind_output;
		goto ret;
		}
	rv = 0;
	s = d->cstat;
	for(i = 0; i < n; ++i) {
		if ((j = s[i] + 3) >= 0 && j <= 3)
			s[i] = vmap[j];
		else {
			badretfmt(504, "Surprise VBasis[%d] = %d.", i, j-3);
			goto ret;
			}
		}
	s = d->rstat;
	for(i = 0; i < m; ++i) {
		j = ++s[i];
		if (j < 0 || j > 1) {
			badretfmt(505, "Surprise CBasis[%d] = %d.", i, j-1);
			goto ret;
			}
		}
 ret:
	if (rv)
		dpf(d, "\nNo basis.");
	}

 static void
nl_iv_adj(ASL *asl, int j, int k, char *vtype, real *x)
{
	/* This will be needed once gurobi can handle nonlinear discrete variables, */
	/* e.g., in QPs. */

	int i, i0;
	real *L, *U;

	L = LUv;
	U = Uvx;
	i0 = k - j;
	if (vtype)
		for(i = i0; i < k; ++i)
			vtype[i] = L[i] == 0. && U[i] == 1. ? 'B' : 'I';
	if (x)
		for(i = i0; i < k; ++i) {
			if (x[i] < L[i])
				x[i] = L[i];
			else if (x[i] > U[i])
				x[i] = U[i];
			else
				x[i] = floor(x[i] + .5);
			}
	}

 typedef struct
Sensname {
	char *aname, *gname;
	int iscon;
	} Sensname;

 static void
put_sens(ASL *asl, GRBmodel *mdl)
{
	Sensname *sn, *sne;
	int len[2], nc, nv;
	real *a;
	static Sensname Snames[] = {
		{ "senslbhi",  "SALBUp", 0 },
		{ "senslblo",  "SALBLow", 0 },
		{ "sensobjhi", "SAObjUp", 0 },
		{ "sensobjlo", "SAObjLow", 0 },
		{ "sensrhshi", "SARHSUp", 1 },
		{ "sensrhslo", "SARHSLow", 1 },
		{ "sensubhi",  "SAUBUp", 0 },
		{ "sensublo",  "SAUBLow", 0 }};
	static int ak[2] = {ASL_Sufkind_var, ASL_Sufkind_con};

	len[0] = nv = n_var;
	len[1] = nc = n_con;
	a = (real*)M1alloc((6*nv + 2*nc)*sizeof(real));
	for(sn = Snames, sne = sn + sizeof(Snames)/sizeof(Sensname); sn < sne; ++sn) {
		if (GRBgetdblattrarray(mdl, sn->gname, 0, len[sn->iscon], a))
  	  		namefailed("GRBgetdblattrarray", sn->gname);
		suf_rput(sn->aname, ak[sn->iscon], a);
		a += len[sn->iscon];
		}
	}

#if GRB_VERSION_MAJOR >= 3 /*{*/
 static int
ams_write(ASL *asl, GRBenv *env, GRBmodel *mdl, Dims *dims, int nsols,
	  real bestobj, int objprec)
{
	Option_Info oi;
	char *fname, *fname_end, msg[64];
	enum {fname_endlen = 32};
	int havetol, i, j, nvr;
	real ba, *c, obj, t, ta, *x;
	size_t L;

	memset(&oi, 0, sizeof(oi));
	oi.wantsol = 9;
	if (--nsols > ams_limit && ams_limit > 0)
		nsols = ams_limit;
	nvr = n_var + nranges;
	L = strlen(ams_stub);
	x = (real*)Malloc(nvr*sizeof(real) + L + fname_endlen);
	fname = (char*)(x + nvr);
	fname_end = fname + L;
	memcpy(fname, ams_stub, L);
	havetol = 0;
	ba = 0.; /* silence erroreous warning */
	if (ams_eps > 0. || ams_epsabs > 0.) {
		havetol = 1;
		if ((ba = bestobj) < 0.)
			ba = -ba;
		}
	c = dims->c;
	for(i = 1; i <= nsols; ++i) {
		if (GRBsetintparam(env, "SolutionNumber", i))
			namefailed("GRBsetintparam", "SolutionNumber");
		if (GRBgetdblattrarray(mdl, GRB_DBL_ATTR_Xn, 0, nvr, x))
			namefailed("GRBgetdblattrarray", GRB_DBL_ATTR_Xn);
		if (!GRBgetdblattr(mdl, GRB_DBL_ATTR_OBJVAL, &obj)) {
			obj = 0.;
			for(j = 0; j < nvr; ++j)
				obj += c[j]*x[j];
			}
		if (havetol) {
			t = dims->objsense*(obj - bestobj);
			if (ams_epsabs > 0. && t > ams_epsabs)
				break;
			if (ams_eps > 0. && t > 0.) {
				if ((ta = obj) < 0.)
					ta = -ta;
				if (ta < ba)
					ta = ba;
				if (t/ta > ams_eps)
					break;
				}
			}
		Snprintf(msg, sizeof(msg), "Alternative MIP solution %d, objective = %.*g",
			i, objprec, obj);
		Snprintf(fname_end, fname_endlen, "%d.sol", i);
		if (write_solf_ASL(asl, msg, x, 0, &oi, fname))
			break;
		dpf(dims, "\n%s", msg);
		}
	free(x);
	if (GRBsetintparam(env, "SolutionNumber", 0)) /* restore */
		namefailed("GRBsetintparam", "SolutionNumber");
	return i - 1;
	}
#endif /*}*/

#define MBL 8192

 static GRBenv *env0;

 FILE *GRB_nl_ASL;	/* Foil gcc optimization bug: with -O, the "nl = 0" */
			/* assignment below does not happen, and after longjmp */
			/* we get an erroneous attempt to fclose(nl). */

 int
main(int argc, char **argv)
{
	ASL *asl;
	Dims dims;
	/*FILE *nl;*/ #define nl GRB_nl_ASL
	Filename *fn;
	GRBenv *env;
	GRBmodel *mdl, *mdl1;
	char mbuf[MBL];
	char *hx0, *sense, *smsg, *sostype, *stub, *vtype;
	const char *solmsg;
	int i, j, k, lvi, nc, nfree, nlvi, nqc, nsosnz, nrange, nsos;
	int nv, nvr, nz, nzcr, objprec, rc, wantobj;
	int *cs, *csr, *rnr, *rsta, *sosbeg, *sosind, *sostypes;
	int vinfo[3], *vlen, *vlenr;
	ograd *og;
	real absmipgap, f, obj, oc, relmipgap, t;
	real *A, *Ar, *bb, *lu, *lxr, *rhs, *sosref, *uxr, *x, *y;
	sig_func_type *oic;
	size_t L;
	static int sos_types[2] = { GRB_SOS_TYPE1, GRB_SOS_TYPE2 };
#if GRB_VERSION_MAJOR >= 4 /*{*/
	fint *colqf, nelqf, *rowqf;
	int *colq, i1, j1, nelq, *rowq;
	real *qmat;
#if GRB_VERSION_MAJOR >= 5 /*{*/
	cgrad *cg, **cgp;
	char qsense;
	const char *objqual = "";
	fint *colqc, nelqc, *rowqc;
	int *ia, *rn, *zc, *zr;
	int k1, n, nlnz;
	real *a1, fto, *pfto, *qmatc, qrhs, *resid;
#ifdef GRB_INT_PAR_TUNEOUTPUT
	char *tunemsg;
#endif
#endif /*}*/
#endif /*}*/

	nelqf = 0;
	Times[0] = xectim_();
	Times[4] = time(0);
#ifdef LICENSE_FILE
	if (!(stub = getenv("GRB_LICENSE_FILE")) || !*stub)
		putenv("GRB_LICENSE_FILE=" LICENSE_FILE);
#endif

	oic = 0;
	env0 = 0;
	mdl = 0;
	memset(&dims, 0, sizeof(dims));
	vinfo[0] = vinfo[1] = vinfo[2] = 0;
	GRBversion(&vinfo[0], &vinfo[1], &vinfo[2]);
	Snprintf(verbuf, sizeof(verbuf), "Gurobi %d.%d.%d", vinfo[0], vinfo[1], vinfo[2]);
	Lic_info_add_ASL = "Portions Copyright Gurobi Optimization, Inc., 2008.";

#ifdef main
	if (!(asl = asl1))
#endif
	asl = ASL_alloc(ASL_read_fg);
	if (!(stub = getstub(&argv, &Oinfo)))
		usage_ASL(&Oinfo, 1);
	nl = jac0dim(stub, 0);
	nqc = nlc;
#if GRB_VERSION_MAJOR < 5
	if (nqc > 0) {
		asl->i.uinfo = "Gurobi can't handle nonlinear constraints.";
		solve_result_num = 522;
		rc = 1;
		goto bailout;
		}
#endif
	if (n_cc > 0) {
		asl->i.uinfo = "Gurobi can't handle complementarity constraints.";
		solve_result_num = 567;
		rc = 1;
		goto bailout;
		}
	if (!(nobjno = n_obj))
		objno = 0;
	rc = 1;
	/* Passing 0 file logfile; would make logfile an option if we could	*/
	/* specify it after processing $gurobi_options, but we cannot do so,	*/
	/* and we need to have env before reading $gurobi_options, as values	*/
	/* in env may be changed by $gurobi_options. */
	if ((!env0 && GRBloadenv(&env0,0)) || !env0) {
		solve_result_num = 500;
		solmsg = "Could not create the gurobi environment.";
		goto ws_now;
		}
	Oinfo.uinfo = (char*)env0;
	nlvi = nlvbi + nlvci + nlvoi;
	lvi = nbv + niv;
	dims.kiv = nlvi + lvi;
	GRBsetintparam(env0, "OutputFlag", 0); /* possibly changed by getopts... */
	rc = setjmp(Jb);
	if (rc) {
 bailout:
		if (nl)
			fclose(nl);
		--rc;
		if (solve_result_num > 0 && asl->i.uinfo){
			if (amplflag | (Oinfo.wantsol & 1))
				rc = 0;
			L = strlen(Oinfo.bsname) + strlen(asl->i.uinfo);
			solmsg = smsg = (char*)M1alloc(L+3);
			sprintf(smsg, "%s: %s", Oinfo.bsname, asl->i.uinfo);
 ws_now:
			write_sol(solmsg, 0, 0, &Oinfo);
			}
		goto done;
		}
	if (getopts(argv, &Oinfo)) {
		if ((solmsg = asl->i.uinfo))
			goto ws_now;
		solve_result_num = 503;
		solmsg = "Bad $gurobi_options.";
		goto ws_now;
		}
#ifdef ALLOW_GUROBI_SERVER
	if (serverlic && server_licread(asl)) {
		solmsg = asl->i.uinfo;
		goto ws_now;
		}
	if (server) {
		if (env0) {
			GRBfreeenv(env0);
			env0 = 0;
			}
		if ((i = GRBloadclientenv(&env0, logfile, server, server_port,
				server_passwd, server_priority, server_timeout))) {
			switch(i) {
			 case GRB_ERROR_NETWORK:
				solmsg = "Could not talk to Gurobi Compute Server.";
				solve_result_num = 601;
				break;
			 case GRB_ERROR_JOB_REJECTED:
				solmsg = "Job rejected by Gurobi Compute Server.";
				solve_result_num = 602;
				break;
			 case GRB_ERROR_NO_LICENSE:
				solmsg = "No license for specified Gurobi Compute Server.";
				solve_result_num = 603;
				break;
			 default:
				solmsg = mbuf;
				snprintf(mbuf, sizeof(mbuf),
					"Surprise return %d from GRBloadclientenv().", i);
				solve_result_num = 604;
			 }
			goto ws_now;
			}
		Oinfo.uinfo = (char*)env0;
		logfile = 0;
		}
#endif
	if (relax)
		dims.kiv = 0;
	breaking = 3;
	oic = signal(SIGINT, intcatch);
	nrange = nranges;
	dims.nv0 = nv = n_var;
	nvr = nv + nrange;
	dims.nc0 = nc = n_con;
	nz = nzc;
	nzcr = nz + nrange;
	L = (2*(nvr+nc)+nzcr+nrange) * sizeof(real)
		+ (2*nvr + nzcr + 1)*sizeof(int) + nc + nvr + nv;
	A_vals = A = (real*)Malloc(L);
	Ar = A + nz;
	LUv = lxr = A + nzcr;
	Uvx = uxr = lxr + nvr;
	y = uxr + nvr;
	rhs = y + nc;
	vlen = vlenr = (int*)(rhs + nc);
	A_colstarts = cs = vlen + nvr;
	A_rownos = rnr = cs + nvr + 1;
	sense = (char*)(rnr + nzcr);
	vtype = sense + nc;
	havex0 = hx0 = vtype + nvr;
	suf_declare(suftab, sizeof(suftab)/sizeof(SufDecl));
	dims.cstat = (int*)M1zapalloc((nvr+nc+2)*sizeof(int));
	dims.rstat = dims.cstat + nvr + 1;
	rsta = dims.cstat + nv;
	dims.csd = suf_iput("sstatus", ASL_Sufkind_var, dims.cstat);
	dims.rsd = suf_iput("sstatus", ASL_Sufkind_con, dims.rstat);

#if GRB_VERSION_MAJOR >= 5
	resid = 0;
	if (dims.kiv) {
		want_xpi0 = 1;
		warmstart = 0;
		}
	else if (warmstart) {
		want_xpi0 = 3;
		if (nrange) {
			pi0 = 0;
			resid = y;
			asl->i.nsufext[ASL_Sufkind_var] += nrange;
			/* include space for range slacks in X0 */
			}
		}
#else
	if (dims.kiv)
		want_xpi0 = 1;
#endif
#if GRB_VERSION_MAJOR >= 4
	qp_read(nl,0);
#else
	fg_read(nl,0);
#endif
	nl = 0;	/* was closed by qp_read */
	nfree = nsos = 0;
	dims.c = 0;
	dims.objsense = 1;
	oc = 0.;
	if (objno > 0) {
		i = objno - 1;
#if GRB_VERSION_MAJOR < 4
		if (i >= n_obj - nlo) {
			asl->i.uinfo = "Gurobi cannot handle nonlinear objectives.";
			solve_result_num = ASL_readerr_nonlin;
			rc = 1;
			goto bailout;
			}
#else
		/* must call mqpcheck() after qp_read so objconst() will work right */
		nelqf = mqpcheck(i, &rowqf, &colqf, &qmat);
		if (nelqf < 0) {
			if (nelqf == -2) {
				solve_result_num = 523;
				asl->i.uinfo = "Cannot handle a quadratic objective involving division by 0";
				}
			else {
				solve_result_num = 521;
				asl->i.uinfo = "Gurobi cannot handle general nonlinear objectives.";
				}
			rc = 1;
			goto bailout;
			}
#endif
		if (objtype[i])
			dims.objsense = -1;
		oc = objconst(i);
		dims.c = M1zapalloc(nvr * sizeof(real));
		for(og = Ograd[i]; og; og = og->next)
			dims.c[og->varno] = og->coef;
		}
	if (dims.kiv || (sos
	    && suf_get("sosno", ASL_Sufkind_var | ASL_Sufkind_input)
	    && suf_get("ref", ASL_Sufkind_var | ASL_Sufkind_input))) {
		i = ASL_suf_sos_explict_free;
		if (!sos)
			i |= ASL_suf_sos_ignore_sosno;
		if (!sos2)
			i |= ASL_suf_sos_ignore_amplsos;
		if ((nsos = suf_sos(i, &nsosnz, &sostype, 0, 0,
				&sosbeg, &sosind, &sosref))) {
			nv = n_var;
			nvr = nv + nrange;
			nc = n_con;
			nz = nzc;
			nzcr = nz + nrange;
			lvi = nbv + niv;
			}
		}
	Ar = A + nz;
	rnr = A_rownos + nz;
	lxr = LUv + nv;
	uxr = Uvx + nv;
	rsta = dims.cstat + nv;
#if GRB_VERSION_MAJOR >= 5 /*{*/
	if (!pi0) {
		if (nc) {
			warmstart = 0;
			memset(y, 0, nc*sizeof(real));
			}
		}
	else if (warmstart && nrange && (x = X0)) {
		memset(resid, 0, nc*sizeof(real));
		rn = A_rownos;
		for(i = j = 0; i < nv; ) {
			t = x[i];
			for(k = cs[++i]; j < k; ++j)
				resid[rn[j]] += t*A[j];
			}
		}
#endif /*}*/
	for(i = 0; i < nv; ++i)
		*vlenr++ = cs[i+1] - cs[i];
	csr = cs + nv;
	for(lu = LUrhs, i = 0; i < nc; ++i, lu += 2) {
		rhs[i] = lu[0];
		if (lu[0] <= negInfinity) {
			if (lu[1] >= Infinity)
				++nfree;
			else {
				sense[i] = '<';
				rhs[i] = lu[1];
				}
			}
		else if (lu[1] >= Infinity)
			sense[i] = '>';
		else {
			sense[i] = '=';
			if (lu[1] > lu[0]) {
				*rnr++ = i;
				*csr++ = Ar - A;
				*Ar++ = -1.;
				*vlenr++ = 1;
				*lxr++ = 0.;
				*uxr++ = lu[1] - lu[0];
				*rsta++ = y[i] > 0. ? -1 : y[i] < 0. ? -2 : 0;
				}
			}
		}
	if (nfree) {
		fprintf(stderr, "Botch: gurobi cannot handle %d free rows.\n", nfree);
		return 1;
		}
	if (nvr)
		*csr = Ar - A;

	memset(vtype, 'C', nvr);
	if (dims.kiv) {
		if (nlvi) {
			k = nlvb;
			if ((j = nlvbi))
				nl_iv_adj(asl, j, k, vtype, X0);
			k = nlvc;
			if ((j = nlvci))
				nl_iv_adj(asl, j, k, vtype, X0);
			k += nlvo - nlvc;
			if ((j = nlvoi))
				nl_iv_adj(asl, j, k, vtype, X0);
			}
		k = nv - lvi;
		if ((j = nbv)) {
			memset(vtype+k, 'B', j);
			k += j;
			if (X0)
				nl_iv_adj(asl, j, k, 0, X0);
			}
		if ((j = niv)) {
			memset(vtype+k, 'I', j);
			k += j;
			if (X0)
				nl_iv_adj(asl, j, k, 0, X0);
			}
		}
#if GRB_VERSION_MAJOR >= 5
	if (nqc) {
		if (dims.kiv && (GRBgetintparam(env0, "QCPDual", &i) || i == 0))
			basis = solnsens = 0;
		ia = A_rownos;
		j = nzc;
		for(i = nlnz = 0; i < j; i++)
			if (ia[i] < nqc)
				++nlnz;
		if ((n = nlvc) < nlvo)
			n = nlvo;
		for(i = k = 0; i < nqc; ++i) {
			nelqc = mqpcheck(-(i+1), 0, 0, 0);
			if (k < nelqc)
				k = nelqc;
			}
		a1 = (real*)Malloc(nqc*sizeof(cgrad*) + nlnz*sizeof(cgrad)
				+ nv*(sizeof(int) + sizeof(real)) + k*(2*sizeof(int)));
		cg = (cgrad*)(a1 + nv);
		Cgrad = cgp = (cgrad**)(cg + nlnz);
		zc = (int*)(cgp + nqc);
		zr = zc + k;
		rn = zr + k;
		memset(cgp, 0, nqc*sizeof(cgrad*));
		k = cs[i = nvr];
		while(--i >= 0) {
			j = cs[i];
			j1 = k1 = k;
			while(k > j) {
				if ((i1 = ia[--k]) < nqc) {
					cg->varno = i;
					cg->coef = A[k];
					cg->next = cgp[i1];
					cgp[i1] = cg++;
					}
				else  {
					A[--j1] = A[k];
					ia[j1] = i1 - nqc;
					}
				}
			cs[i] = j1;
			vlen[i] = k1 - j1;
			}
		if (GRBloadmodel(env0, &mdl, "foo",
				nvr, nc-nqc, dims.objsense, oc, dims.c,
				sense+nqc, rhs+nqc, cs, vlen, ia, A,
				LUv, Uvx, vtype,  0,0) || !mdl)
			failed(env0, "GRBloadmodel");
		lu = LUrhs;
		for(i = 0; i < nqc; ++i, lu += 2) {
			nelqc = mqpcheck(-(i+1), &rowqc, &colqc, &qmatc);
			if (nelqc <= 0) {
				switch(nelqc) {
				  case 0:
					solve_result_num = 502;
					asl->i.uinfo =
					 "no quadratic terms in a \"quadratic\" constraint.";
					break;
				  case -2:
					solve_result_num = 525;
					asl->i.uinfo =
					 "a quadratic constraint involving division by 0.";
					break;
				  default:
					solve_result_num = 522;
					asl->i.uinfo =
					 "Gurobi can't handle nonquadratic nonlinear constraints.";
					}
				rc = 1;
				goto bailout;
				}
			if (lu[0] <= negInfinity) {
				qsense = '<';
				qrhs = lu[1];
				}
			else if (lu[1] >= Infinity) {
				qsense = '>';
				qrhs = lu[0];
				}
			else {
				if (lu[0] == lu[1]) {
					solve_result_num = 526;
					asl->i.uinfo =
					 "Gurobi cannot handle quadratic equality constraints.";
					}
				else {
					solve_result_num = 527;
					asl->i.uinfo =
					 "Gurobi cannot handle quadratic range constraints.";
					}
				rc = 1;
				goto bailout;
				}
			for(i1 = j = k = 0; j < n; ++j) {
				for(j1 = colqc[j+1]; i1 < j1; ++i1) {
					k1 = rowqc[i1];
					if (k1 <= j) {
						zc[k] = j;
						zr[k] = k1;
						if (j == k1)
							qmatc[i1] *= 0.5;
						qmatc[k++] = qmatc[i1];
						}
					}
				}
			for(j = 0, cg = cgp[i]; cg; cg = cg->next)
				if (cg->coef != 0.) {
					rn[j] = cg->varno;
					a1[j++] = cg->coef;
					}
			if (GRBaddqconstr(mdl, j, rn, a1,
						k, zr, zc, qmatc, qsense, qrhs, 0))
				failed(env0, "GRBaddqconstr");
			free(colqc);
			free(rowqc);
			free(qmatc);
			}
		free(a1);
		}
	else
#endif
	if (GRBloadmodel(env0, &mdl, "foo",
			nvr, nc, dims.objsense, oc, dims.c,
			sense+nqc, rhs+nqc, cs, vlen, A_rownos, A_vals,
			LUv, Uvx, vtype,  0,0) || !mdl)
		failed(env0, "GRBloadmodel");
	x = 0;
	if (X0 && (!dims.kiv || mipstval)) {
		x = X0;
		for(i = 0; i < nv; ++i)
			if (!hx0[i]) {
				x[i] = GRB_UNDEFINED;
				warmstart = 0;
				}
		}

	if (!(env = GRBgetenv(mdl)))
		failed(env0, "GRBgetenv");
	Oinfo.uinfo = (char*)env;

#if GRB_VERSION_MAJOR >= 4
	if (nelqf) {
		for(i = j = nelq = 0; i < nv; ++i) {
			for(k = colqf[i+1]; j < k; ++j)
				if (rowqf[j] <= i)
					++nelq;
			}
		colq = (int*)Malloc(nelq*sizeof(int));
		for(i = j = j1 = 0; i < nv; ++i) {
			for(k = colqf[i+1]; j < k; ++j) {
				if ((i1 = rowqf[j]) <= i) {
					colq[j1] = i;
					rowqf[j1] = i1;
					qmat[j1] = qmat[j];
					if (i1 == i)
						qmat[j1] *= 0.5;
					++j1;
					}
				}
			}
		if (sizeof(fint) == sizeof(int))
			rowq = (int*)rowqf;
		else {
			rowq = (int*)Malloc(nelq*sizeof(int));
			for(i = 0; i < nelq; ++i)
				rowq[i] = rowqf[i];
			free(rowqf);
			}
		free(colqf);
		if (GRBaddqpterms(mdl, nelq, rowq, colq, qmat))
			failed(env, "GRBaddqpterms");
		free(rowq);
		free(qmat);
		}
#endif
	if (nsos) {
		sostypes = (int*)Malloc(nsos*sizeof(int));
		for(i = 0; i < nsos; ++i)
			sostypes[i] = sos_types[sostype[i] - '1'];
		if (GRBaddsos(mdl, nsos, nsosnz, sostypes, sosbeg, sosind, sosref))
			failed(env, "GRBaddsos");
		free(sostypes);
  	  	}
#if GRB_VERSION_MAJOR >= 5 /*{{*/
	if (!(basis & 1 && get_input_statuses(asl, env, mdl, &dims))
	 && warmstart && x && (nc == 0 || pi0 != 0)) {
		if (nrange) {
			lu = LUrhs;
			j = nv;
			for(i = 0; i < nc; ++i, lu += 2) {
				if (lu[0] > negInfinity
				 && lu[1] < Infinity
				 && lu[0] < lu[1])
					x[j++] = resid[i] - rhs[i];
				}
			}
		if (GRBsetdblattrarray(mdl, GRB_DBL_ATTR_PSTART, 0, nvr, x))
			failed(env, "GRBsetdblattrarray(PStart)");
		if (pi0
		 && nc > nqc
		 && GRBsetdblattrarray(mdl, GRB_DBL_ATTR_DSTART, 0, nc-nqc, pi0))
			failed(env, "GRBsetdblattrarray(DStart)");
		}
#else /*}{*/
	if (basis & 1)
		get_input_statuses(asl, env, mdl, &dims);
#endif /*}}*/
	free(A);
	if (x && GRBsetdblattrarray(mdl, GRB_DBL_ATTR_START, 0, nv, x))
		failed(env, "GRBsetdblattrarray(START)");
	if (logfile) {
		if (GRBsetintparam(env, "OutputFlag", 1))
			namefailed("GRBsetintparam", "OutputFlag");
#if GRB_VERSION_MAJOR == 1
		if (GRBsetlogfilename(env, logfile))
			failed(env, "GRBsetlogfilename");
#elif GRB_VERSION_MAJOR == 2
		if (GRBsetstrparam(env, "LogfileName", logfile))
			namefailed("GRBsetstrparam", "LogfileName");
#else
		if (GRBsetstrparam(env, GRB_STR_PAR_LOGFILE, logfile))
			namefailed("GRBsetstrparam", GRB_STR_PAR_LOGFILE);
#endif
		}
	else if (outlev && GRBsetintparam(env, "OutputFlag", 1))
		namefailed("GRBsetintparam", "OutputFlag");
	if ((fn = Wflist[0])) {
		GRBupdatemodel(mdl);
		do {
			if (GRBwrite(mdl, fn->name))
				enamefailed(env, "GRBwrite", fn->name);
			} while((fn = fn->next));
		}
#if (GRB_VERSION_MAJOR == 4 && GRB_VERSION_MINOR >= 5) || GRB_VERSION_MAJOR >= 5 /*{*/
	if (rays)
		GRBsetintparam(env, "InfUnbdInfo", 1);
#endif /*}*/
#ifdef GRB_INT_ATTR_BRANCHPRIORITY
	if (dims.kiv && priorities)
		add_priorities(asl, env, mdl);
#endif
#if GRB_VERSION_MAJOR >= 5
	pfto = 0;
	if (feasrelax)
		pfto = do_feasrelax(asl, env, mdl, &objqual, &fto);
#endif
	breaking = 1;
	Times[1] = xectim_();
#ifdef GRB_INT_PAR_TUNEOUTPUT
	tunemsg = 0;
	if (tunebase && tunerun(asl, env, mdl, &tunemsg)) {
		write_sol(tunemsg, 0, 0, &Oinfo);
		goto done;
		}
#endif
	grbmodel = mdl;
	i = GRBoptimize(mdl);
	grbmodel = 0;
	Times[2] = xectim_();
	solmsg = statmsg(asl, env, mdl, i, &dims, &wantobj);
	dims.mb = mbuf;
	dims.mbend = mbuf + sizeof(mbuf);
	dpf(&dims, "%s: %s", verbuf, solmsg);
	absmipgap = relmipgap = Infinity;
	objprec = 0;
	if (wantobj && !GRBgetdblattr(mdl, GRB_DBL_ATTR_OBJVAL, &obj)) {
#if GRB_VERSION_MAJOR >= 5
		dpf(&dims, "; %sobjective %.*g", objqual, objprec = obj_prec(), obj);
		if (pfto)
			dpf(&dims, "\nfeasrelax objective = %.*g", objprec, *pfto);
#else
		dpf(&dims, "; objective %.*g", objprec = obj_prec(), obj);
#endif
		if ((bestbound | (retmipgap^4)) && objno > 0) {
			if (GRBgetdblattr(mdl, GRB_DBL_ATTR_OBJBOUND, &f))
				f = Infinity * dims.objsense;
			else {
				if ((absmipgap = obj - f) < 0.)
					absmipgap = -absmipgap;
				if ((t = obj) < 0.)
					t = -t;
				relmipgap = absmipgap / (1e-10 + t);
				}
			if (retmipgap & 1) {
				bb = (real*)M1zapalloc(nobjno*sizeof(real));
				bb[objno - 1] = relmipgap;
				suf_rput("relmipgap", ASL_Sufkind_obj, bb);
				suf_rput("relmipgap", ASL_Sufkind_prob, bb);
				}
			if (retmipgap & 2) {
				bb = (real*)M1zapalloc(nobjno*sizeof(real));
				bb[objno - 1] = absmipgap;
				suf_rput("absmipgap", ASL_Sufkind_obj, bb);
				suf_rput("absmipgap", ASL_Sufkind_prob, bb);
				}
			if (bestbound) {
				bb = (real*)M1zapalloc(nobjno*sizeof(real));
				bb[objno - 1] = f;
				suf_rput("bestbound", ASL_Sufkind_obj, bb);
				suf_rput("bestbound", ASL_Sufkind_prob, bb);
				}
			}
		}
	else
		solnsens = 0;
#if GRB_VERSION_MAJOR >= 3
	if (!GRBgetintattr(mdl, GRB_INT_ATTR_BARITERCOUNT, &i) && i > 0)
		dpf(&dims, "\n%d barrier iterations", i);
#endif
	if (!GRBgetdblattr(mdl, GRB_DBL_ATTR_ITERCOUNT, &f) && f > 0.)
		dpf(&dims, "\n%.0f simplex iterations", f);
	if (dims.kiv && !GRBgetdblattr(mdl, GRB_DBL_ATTR_NODECOUNT, &f) && f > 0.)
		dpf(&dims, "\n%.0f branch-and-cut nodes", f);
	mdl1 = mdl;
	if (dims.kiv && ((basis & 2) | solnsens))
		mdl1 = fixed_model(asl, mdl, &dims);
	if (absmipgap > 0. && absmipgap < Infinity && !(retmipgap & 4))
		dpf(&dims, "\nabsmipgap = %.3g, relmipgap = %.3g", absmipgap, relmipgap);
	if (mdl1) {
		if (solnsens)
			put_sens(asl, mdl1);
		if (basis & 2)
			get_output_statuses(asl, mdl1, &dims);
		if (mdl1 != mdl)
			GRBfreemodel(mdl1);
		}
	if (dims.missing)
		missing_msg(&dims);
#if GRB_VERSION_MAJOR >= 3 /*{*/
	if (dims.kiv && ams_stub > 0 && wantobj
	 && !GRBgetintattr(mdl, GRB_INT_ATTR_SOLCOUNT, &i) && i > 1
	 && (k = ams_write(asl, env, mdl, &dims, i, obj, objprec))) {
		dpf(&dims, "\n%d alternative MIP solution%s written to %s1.sol",
			k, "s" + (k == 1), ams_stub);
		if (k >= 2)
			dpf(&dims, "\n%s %s%d.sol", k == 2 ? "and" : "...", ams_stub, k);
		if (k < --i)
			dpf(&dims, "\nIgnoring %d other inferior alternative MIP solutions.",
				i - k);
		}
#ifdef GRB_INT_PAR_TUNEOUTPUT
	if (tunemsg)
		dpf(&dims, "\n%s", tunemsg);
#endif
#endif /*}*/
	write_sol(mbuf, dims.x, dims.y, &Oinfo);
	for(fn = Wflist[1]; fn; fn = fn->next)
		if (GRBwrite(mdl, fn->name))
			enamefailed(env, "GRBwrite", fn->name);
 done:
	if (oic)
		signal(SIGINT, oic);
	if (mdl)
		GRBfreemodel(mdl);
	if (env0)
		GRBfreeenv(env0);
	if (!amplflag && solve_result_num >= 500)
		rc = 1;
	ASL_free(&asl);
	show_times();
	return rc;
	}
