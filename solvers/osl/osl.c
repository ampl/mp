/****************************************************************
Copyright (C) 1992-1999, 2001 Lucent Technologies
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

/* #define noDEBUG to omit the "fortran" keyword ("fortran=foo"
		makes osl write Fortran to file foo.  The Fortran
		currently omits the kludge needed for recovering
		the solution with "simplify" >= 0 when there are
		integer variables.) */

#ifndef MAX_PRIORITY
/* default = 2^31 - 2, to allow expressing MAX_PRIORITY  + 1*/
#define MAX_PRIORITY 2147483646
#endif

#ifdef KR_headers
#define remove unlink
#define Void /* */
#else
#define Void void
#endif
#include "limits.h"
#include "getstub.h"
#define asl cur_ASL

#ifdef Bad_DBL_MAX
/* Guard against bozo C compilers, e.g. Borland,
   for which DBL_MAX is not a constant expression.
*/
#undef DBL_MAX
#define DBL_MAX 1.7e38	/* safe even for VAX arithmetic */
#else
#include "float.h"	/* for DBL_MAX */
#endif

#ifndef OSL_V2
#undef  OSL_V3
#define OSL_V3
#endif/*OSL_V2*/

#include "osl.h"

#ifndef KR_headers
/* for getpid()... */
#ifdef WIN32
#include "process.h"
#else
#include "unistd.h"
#endif
#endif

#ifndef OSLLIBAPI
#define OSLLIBAPI extern
#endif

#ifndef OSLLINKAGE
#define OSLLINKAGE
#endif

#ifndef OSLCALLBACK
#define OSLCALLBACK OSLLINKAGE
#endif

#define OSLvoid OSLLIBAPI void OSLLINKAGE

typedef void OSLCALLBACK Cbfunc ANSI((double*, fint*, double*, fint*, fint*));

OSLvoid ekkbasi_ ANSI((fint*, double*, fint*));
OSLvoid ekkbaso_ ANSI((fint*, double*, fint*, fint*));
OSLvoid ekkbcdo_ ANSI((fint*, double*, fint*, fint*, fint*));
OSLvoid ekkbmpr_ ANSI((fint*, double*, fint*, fint*));
OSLvoid ekkbmps_ ANSI((fint*, double*, fint*));
OSLvoid ekkbslv_ ANSI((fint*, double*, fint*, fint*));
OSLvoid ekkclcb  ANSI((fint*));
OSLvoid ekkcrsh_ ANSI((fint*, double*, fint*));
OSLvoid ekkdsca_ ANSI((fint*, double*, fint*, fint*));
OSLvoid ekkdscm_ ANSI((fint*, double*, fint*, fint*));
OSLvoid ekkfcls_ ANSI((fint*, char*, fint*, fint*));
#if defined(OSLWCSTAT) || defined(OSLWCDLL)
/* idiocy */
OSLvoid ekkfopn_ ANSI((fint*, char*, char*, char*, char*, char*,
			fint*, fint*, fint*, fint*, fint*, fint*, fint*));
#define Ekkfcls(a1,a2,a3,a4) ekkfcls_(a1,a2,a4,a3)
#define Ekkfopn(a1,a2,a3,a4,a5,a6,a7,a8,a9,aa,ab,ac,ad) \
	ekkfopn_(a1,a2,a4,a6,a8,aa,ac,ad,a3,a5,a7,a9,ab)
#else
OSLvoid ekkfopn_ ANSI((fint*, char*, fint*, char*, fint*, char*, fint*,
			char*, fint*, char*, fint*, fint*, fint*));
#define Ekkfcls(a1,a2,a3,a4) ekkfcls_(a1,a2,a3,a4)
#define Ekkfopn(a1,a2,a3,a4,a5,a6,a7,a8,a9,aa,ab,ac,ad) \
	ekkfopn_(a1,a2,a3,a4,a5,a6,a7,a8,a9,aa,ab,ac,ad)
#endif
OSLvoid ekkiget_ ANSI((fint*, double*, fint*, fint*));
OSLvoid ekkimdl_ ANSI((fint*, double*, fint*, fint*, fint*, fint*, fint*,
			fint*, fint*, fint*, double *, double *));
OSLvoid ekkinit_ ANSI((fint*, double*));
OSLvoid ekkiset_ ANSI((fint*, double*, fint*, fint*));
OSLvoid ekklmdl_ ANSI((fint *rc, double *dspace, fint *type, fint *nr,
			fint *nc, fint *nel, double *c, double *clb,
			double *cub, double *lb, double *ub,
			fint *ia, fint *ja, double *a));
OSLvoid ekkmpre_ ANSI((fint*, double*, fint*));
OSLvoid ekkmset_ ANSI((fint *rc, double *dspace, fint *strtnum,
			fint *maxalw, fint *maxprt, fint *trace,
			fint *usrexit, fint *endnum, fint *nonum));
OSLvoid ekkmslv_ ANSI((fint*, double*, fint*, fint*, fint*));
OSLvoid ekkname_ ANSI((fint*, double*, fint*, char*, fint*, fint*, char*,
			fint*, fint*));
OSLvoid ekknfes_ ANSI((fint*, double*, fint*, fint*, fint*, fint*));
OSLvoid ekknget_ ANSI((fint*, double*, fint*, fint*));
OSLvoid ekknslv_ ANSI((fint*, double*, fint*, fint*));
OSLvoid ekknwmt_ ANSI((fint*, double*, fint*));
OSLvoid ekkprsl_ ANSI((fint*, double*, fint*, fint*));
OSLvoid ekkpssl_ ANSI((fint*, double*, fint*));
OSLvoid ekkqmdl_ ANSI((fint*, double*, fint*, fint*, fint*, fint*, double*));
OSLvoid ekkqslv_ ANSI((fint*, double*, fint*, fint*));
OSLvoid ekkrgcb  ANSI((fint*, Cbfunc*));
OSLvoid ekkrget_ ANSI((fint*, double*, double*, fint*));
OSLvoid ekkrset_ ANSI((fint*, double*, double*, fint*));
OSLvoid ekksbnd_ ANSI((fint*, double*, fint*, fint*));
OSLvoid ekkscal_ ANSI((fint*, double*));
OSLvoid ekksobj_ ANSI((fint*, double*));
OSLvoid ekksslv_ ANSI((fint*, double*, fint*, fint*));

extern char *getenv ANSI((const char *));

static fint conpri = 2, convex_qp = 1, mipststat = 1, objpri = 3, outlev = 1;
static double osl_infinity = 9.99999999e27, osl_neginfinity = -9.99999999e27;
static int havebasis = 1;
static char *netsign;
static int nshift, shift;
static int intrestore_bug, relaxed_infeas;

 static fint M, N, NO, NZ;
#ifndef noDEBUG
 static FILE *fort;
#endif
 static int	Nchatter, WantPrintf, WantPrintf1;
 static fint	iter_inc = 2000, mbfile = 1, pricing = 2;
 static fint	prestrat = 1, sos = 1, sos2 = 1;
 static fint	branch, bs_switch, dspinc, iisfind, method, mpsout;
 static fint	nelq, netinit, qpnodealg, scale, sensitivity, timing, trace;
 static fint binconv, bdisplay, endbasis, mipdisplay, simplexinit, startbasis;
 static fint mipfreq = 20;
 static double *dspace;
 static double mingap, rowinc;
 static fint dspqfac = 10.;
 static int need_minmax = 1;
 static int binig, binig1, mingapped, need_rm, negobj, negobj1, nircset;
 static int relax, solving;
 static fint barrier = -1, binsimp = 1, crash = -1;
 static fint netalg = 1, simplify = -1;
 typedef struct MBunit { fint unit[2]; fint recl; char *name; } MBunit;
 static MBunit	bunit = { 0, 3, 1024 },
		munit = { 0, 4, 4096 };
 static int dualused, intsols, nonconvex, primaldual = 1;
#define Icontrolsize 61
#define Rcontrolsize 45
 static double	rcontrol[Rcontrolsize], rc1[Rcontrolsize];
 static fint	icontrol[Icontrolsize], ic1[Icontrolsize];
 static char	icset[Icontrolsize+1], rcset[Rcontrolsize+1];
 extern char *progname;
 static fint passx0 = 1;
 static fint statuses = 3;
 static char versbuf[32];

 typedef unsigned Long ULong;

 typedef struct
Statinfo {
	int	*cstat;
	int	*rstat;
	SufDesc	*csd;
	SufDesc *rsd;
	real	*c;
	fint	nsets;
	int	objadj;
	} Statinfo;

 typedef struct
rtol {
	double	tmin, tmax;	/* min, max allowed values */
	double	*val;		/* rcontrol subscript */
	int	rsub;
	} rtol;

 typedef struct
itol {
	fint	tmin, tmax;
	fint	*val;
	int	isub;
	} itol;

 typedef struct
rknown {
	double	t;
	double	*val;
	int	rsub;
	} rknown;

static rtol
	Mingap		= { 0., DBL_MAX, &mingap },
	Rbcutoff	= { -1e20, DBL_MAX, &rc1[15], 16 },
	Rchangeweight	= { 1e-12, 1., &rc1[14], 15 },
	Rcholabstol	= { 1e-30, 1e-6, &rc1[9], 10 },
	Rcholtinytol	= { 1e-30, 1e-6, &rc1[10], 11 },
	Rdccutoff	= { -DBL_MAX, DBL_MAX, &rc1[34], 35 },
	Rdegscale	= { 0., DBL_MAX, &rc1[21], 22 },
	Rdensethr	= { -DBL_MAX, DBL_MAX, &rc1[31], 32 },
	Rdweight	= { 0., 1., &rc1[16], 17 },
	Rfixvar1	= { 0., 1e-3, &rc1[7], 8 },
	Rfixvar2	= { 0., 1e-3, &rc1[8], 9 },
	Rimprove	= { -DBL_MAX, DBL_MAX, &rc1[24], 25 },
	Riweight	= { 0., DBL_MAX, &rc1[23], 24 },
	Rmaxmin		= { -1., 1., &rc1[2], 3 },
	Rmufactor	= { 1e-6, 0.99999, &rc1[3], 4 },
	Rmuinit		= { 1e-20, 1e6, &rc1[30], 31 },
	Rmulimit	= { 1e-16, 1., &rc1[4], 5 },
	Rmulinfac	= { 0., 1., &rc1[11], 12 },
	Rnetsamp	= { 0., 1., &rc1[43], 44 },
	Robjweight	= { 0., 1e8, &rc1[32], 33 },
	Rowinc		= { 0., 1e6, &rowinc},
	Rpdgaptol	= { 1e-12, 0.1, &rc1[40], 41 },
	Rpdstepmult	= { 0.01, 0.999999, &rc1[41], 42 },
	Rpertdiag	= { 0., 1e-6, &rc1[42], 43 },
	Rprintcpu	= { 0, DBL_MAX, &rc1[44], 45 },
	Rprojtol	= { 0., 1., &rc1[12], 13 },
	Rpweight	= { 1e-12, 1e10, &rc1[13], 14 },
	Rrgfactor	= { 1e-6, 0.99999, &rc1[5], 6 },
	Rrglimit	= { 0., 1., &rc1[6], 7 },
	Rstepmult	= { 0.01, 0.99999, &rc1[29], 30 },
	Rtarget		= { -DBL_MAX, DBL_MAX, &rc1[25], 26 },
	Rthreshold	= { 0., DBL_MAX, &rc1[39], 40 },
	Rtoldinf	= { 1e-12, 0.1, &rc1[1], 2 },
	Rtolint		= { 1e-12, 0.1, &rc1[26], 27 },
	Rtolpinf	= { 1e-12, 0.1, &rc1[0], 1 }
	;

static itol
	Barrier		= { -1, 4, &barrier },
	Bdisplay	= { 0, LONG_MAX, &bdisplay },
	Binconv		= { 0, 2, &binconv },
	Binsimp		= { 0, 2, &binsimp },
	Branch		= { 0, 15, &branch },
	Bs_switch	= { 0, 4, &bs_switch },
	Conpri		= { 1, MAX_PRIORITY, &conpri },
	Convex_qp	= { 0, 1, &convex_qp },
	Dspinc		= { 0, LONG_MAX, &dspinc },
	Dspqfac		= { 0, 1000, &dspqfac },
#ifdef Use_basunits
	Endbasis	= { 9, 99, &endbasis },
#endif
	Iadjactype	= { 0, 1, &ic1[13], 14 },
	Icrashtype	= { 0, 6, &crash},
	Idensecol	= { 10, LONG_MAX, &ic1[15], 16 },
	Idevexmode	= { -3, 3, &ic1[16], 17 },
	Idroprowct	= { 1, 30, &ic1[18], 19 },
	Ifastits	= { LONG_MIN, LONG_MAX, &ic1[38], 39 },
	Iformntype	= { 0, 1, &ic1[14], 15 },
	Iheurpass	= { 0, LONG_MAX, &ic1[54], 55 },
	Iisfind		= { 0, 1, &iisfind },
	Ilinelen	= { 60, 150, &ic1[35], 36 },
	Ilogfreq	= { 0, LONG_MAX, &ic1[0], 1 },
	Imajorits	= { 0, 999, &ic1[47], 48 },
	Imaxfactor	= { 0, 999, &ic1[2], 3 },
	Imaxiter	= { 0, LONG_MAX, &ic1[4], 5 },
	Imaxiterb	= { 0, LONG_MAX, &ic1[12], 13 },
	Imaxnodes	= { 0, LONG_MAX, &ic1[39], 40 },
	Imaxprojns	= { 1, 10, &ic1[21], 22 },
	Imaxsols	= { 0, LONG_MAX, &ic1[40], 41 },
	Inullcheck	= { LONG_MIN, LONG_MAX, &ic1[17], 18 },
	Ipossbasis	= { 0, LONG_MAX, &ic1[20], 21 },
	Ipricetype	= { 0, 1, &ic1[60], 61 },
	Isupertol	= { 0, LONG_MAX, &ic1[55], 56 },
	Iter_inc	= { 0, LONG_MAX, &iter_inc },
	Ithreshold	= { 0, LONG_MAX, &ic1[53], 54 },
	Method		= { 0, 2, &method },
	MBfile		= { 0, 1, &mbfile },
#ifdef msgDebug
	Mipdisplay	= { 0, 4, &mipdisplay },
#else
	Mipdisplay	= { 0, 2, &mipdisplay },
#endif
	Mipfreq		= { 1, LONG_MAX, &mipfreq },
	Mpsout		= { 0, 2, &mpsout },
	Mipststat	= { 0, 1, &mipststat },
	Netalg		= { -1, 2, &netalg },
	Netinit		= { 0, 1, &netinit },
	Objno		= { 0, LONG_MAX, 0 },
	Objpri		= { 1, MAX_PRIORITY, &objpri },
	Outlev		= { 0, 4, &outlev },
	Passx0		= { 0, 1, &passx0 },
	Prestrat	= { 0, 3, &prestrat },
	Pricing		= { 1, 3, &pricing },
	Scale		= { 0, 1, &scale },
	Sensitivity	= { 0, 1, &sensitivity },
	Simplexinit	= { 0, 2, &simplexinit },
	Simplify	= { -1, 3, &simplify },
	Sos		= { 0, 1, &sos },
	Sos2		= { 0, 1, &sos2 },
#ifdef Use_basunits
	Startbasis	= { 9, 99, &startbasis },
#endif
	Statuses	= { 0, 3, &statuses },
	Trace		= { 0, 2, &trace },
	set_timing	= { 0, 3, &timing }
	;

static I_Known
	set_primal	= { 1, &primaldual },
	set_dual	= { 2, &primaldual },
	set_relax	= { 1, &relax };

static rknown
	set_min		= {  1.0, &rc1[2], 3 },
	set_max		= { -1.0, &rc1[2], 3 };

 static void
#ifdef KR_headers
ftncls_(n) fint *n;
#else
ftncls_(fint *n)
#endif
{
	static fint I4 = 4;
	fint ierr;
	Ekkfcls(n, "KEEP", &I4, &ierr);
	}

 static void
#ifdef KR_headers
rm(u) MBunit *u;
#else
rm(MBunit *u)
#endif
{
	static fint I6 = 6;
	fint ierr;
	Ekkfcls(&u->unit[1], "DELETE", &I6, &ierr);
	}

 static void
rmtemps(Void)
{
	if (need_rm) {
		need_rm = 0;
		if (mbfile) {
			rm(&munit);
			rm(&bunit);
			}
		}
	}

 static void
#ifdef KR_headers
fexit(rc) int rc;
#else
fexit(int rc)
#endif
{
	rmtemps();
	exit(rc);
	}

 static void
#ifdef KR_headers
cantopen(name, what) char *name, *what;
#else
cantopen(char *name, char *what)
#endif
{
	fprintf(Stderr, "Can't open %s (%s file)\n", name, what);
	fexit(2);
	}

 static void
opn(MBunit *u, char *td, int pid)
{
	fint ierr, nlen;
	char *s;
	static fint I11 = 11, I4 = 4, I6 = 6, I7 = 7;

	s = u->name;
	if (td) {
		strcpy(s, td);
		s += strlen(s);
		if (s[-1] != '/')
			*s++ = '/';
		}
#ifdef WIN32
#define fmt "%x.%ld"
#else
#define fmt "fort%u.%ld"
#endif
	s += Sprintf(s, fmt, pid, (long)u->unit[1]);
	nlen = s - u->name;
	ierr = 0;
	Ekkfopn(u->unit+1, u->name, &nlen, "UNKNOWN", &I7,
		"DIRECT", &I6, "UNFORMATTED", &I11,
		"NULL", &I4, &u->recl, &ierr);
	if (ierr)
		cantopen(u->name, "scratch");
	}

 static char*
#ifdef KR_headers
sf_rknown(oi, kw, v) Option_Info *oi; keyword *kw; char *v;
#else
sf_rknown(Option_Info *oi, keyword *kw, char *v)
#endif
{
	rknown *r = (rknown *)kw->info;
	Not_Used(oi);

	*r->val = r->t;
	need_minmax = 0;
	rcset[r->rsub] = 1;
	return v;
	}

 static void
#ifdef KR_headers
no_num(oi, kw) Option_Info *oi; keyword *kw;
#else
no_num(Option_Info *oi, keyword *kw)
#endif
{
	fprintf(Stderr, "Bad option \"%s\" (expected a number to follow)\n",
		kw->name);
	badopt_ASL(oi);
	}

 static void
#ifdef KR_headers
badnum(oi, kw, rv, v) Option_Info *oi; keyword *kw; char *rv; char *v;
#else
badnum(Option_Info *oi, keyword *kw, char *rv, char *v)
#endif
{
	fprintf(Stderr,
		"Bad option \"%s%s%.*s\" (not in ", /*)*/
		kw->name, oi->eqsign, rv-v, v);
	badopt_ASL(oi);
	}

 static char*
#ifdef KR_headers
sf_int(oi, kw, v) Option_Info *oi; keyword *kw; char *v;
#else
sf_int(Option_Info *oi, keyword *kw, char *v)
#endif
{
	fint t;
	char *rv;
	itol *it = (itol *)kw->info;

	t = strtol(v, &rv, 10);
	if (rv == v)
		no_num(oi, kw);
	else if (t < it->tmin || t > it->tmax) {
		badnum(oi, kw, rv, v);
		fprintf(Stderr, /*(*/ "[%ld, %ld]).\n", (long)it->tmin,
				(long)it->tmax);
		}
	else {
		*it->val = t;
		icset[it->isub] = 1;
		}
	return rv;
	}

 static char*
#ifdef KR_headers
sf_objno(oi, kw, v) Option_Info *oi; keyword *kw; char *v;
#else
sf_objno(Option_Info *oi, keyword *kw, char *v)
#endif
{
	fint t;
	char *rv;
	register itol *it = (itol *)kw->info;

	it->tmax = n_obj;
	it->val = &obj_no;
	t = obj_no;
	obj_no = -1;
	rv = sf_int(oi, kw, v);
	if (obj_no >= 0)
		--obj_no;
	else
		obj_no = t;
	return rv;
	}

 static char*
#ifdef KR_headers
sf_dbl(oi, kw, v) Option_Info *oi; keyword *kw; char *v;
#else
sf_dbl(Option_Info *oi, keyword *kw, char *v)
#endif
{
	double t;
	char *rv;
	register rtol *rt = (rtol *)kw->info;

	t = strtod(v, &rv);
	if (rv == v)
		no_num(oi,kw);
	else if (t < rt->tmin || t > rt->tmax) {
		badnum(oi, kw, rv, v);
		fprintf(Stderr, /*(*/ "[%g, %g]).\n", rt->tmin, rt->tmax);
		}
	else {
		*rt->val = t;
		rcset[rt->rsub] = 1;
		}
	return rv;
	}

 static char*
#ifdef KR_headers
sf_minmax(oi, kw, v) Option_Info *oi; keyword *kw; char *v;
#else
sf_minmax(Option_Info *oi, keyword *kw, char *v)
#endif
{
	need_minmax = 0;
	return sf_dbl(oi, kw, v);
	}

 struct Basinfo {
	fint create;
	fint unit;
	fint *up;
	};
 typedef struct Basinfo Basinfo;
 static Basinfo
	Startbasis = { 0, 9, &startbasis },
	Endbasis = { 1, 10, &endbasis };

#ifdef KR_headers
 static char*
sf_bas(oi, kw, v) Option_Info *oi; keyword *kw; char *v;
#else
 static char*
sf_bas(Option_Info *oi, keyword *kw, char *v)
#endif
{
	Basinfo *bi = (Basinfo *)kw->info;
	fint ierr, len;
	static fint I0 = 0, I10 = 10, I13 = 13, I3 = 3, I4 = 4, I9 = 9;
	char *b, *be, buf[256];

	if (*v <= ' ') {
		fprintf(Stderr, "\"%s%s\": missing file name\n",
			kw->name, oi->eqsign);
		return v;
		}
	be = buf + sizeof(buf) - 1;
	for(b = buf; *v > ' '; v++)
		if (b < be)
			*b++ = *v;
	*b = 0;
	len = b - buf;
	ierr = 0;
	Ekkfopn(&bi->unit, buf, &len, "OLD", &I3,
		"SEQUENTIAL", &I10, "FORMATTED", &I9,
		"NULL", &I4, &I0, &ierr);
	if (ierr && bi->create) {
		ierr = 0;
		Ekkfopn(&bi->unit, buf, &len, "NEW", &I3,
			"SEQUENTIAL", &I10, "FORMATTED", &I9,
			"NULL", &I4, &I0, &ierr);
		}
	if (ierr)
		cantopen(buf, kw->name);
	else
		*bi->up = bi->unit;
	return v;
	}

 static char*
#ifdef KR_headers
sf_file(oi, kw, v) Option_Info *oi; keyword *kw; char *v;
#else
sf_file(Option_Info *oi, keyword *kw, char *v)
#endif
{
	FILE *f;
	char *b, *be, buf[256];

	Not_Used(oi);

	if (*v <= ' ') {
		*(FILE **)kw->info = stdout;
		return v;
		}
	be = buf + sizeof(buf) - 1;
	for(b = buf; *v > ' '; v++)
		if (b < be)
			*b++ = *v;
	*b = 0;
	f = fopen(buf, "w");
	if (!f)
		cantopen(buf, kw->name);
	*(FILE **)kw->info = f;
	return v;
	}

#define VP (Char*)
 static keyword keywds[] = { /* must be sorted */
	{ "adjactype",	sf_int,		VP &Iadjactype },
	{ "barrier",	sf_int,		VP &Barrier },
	{ "bb_bfile",	sf_int,		VP &MBfile },
	{ "bb_file",	sf_int,		VP &MBfile },
	{ "bb_mfile",	sf_int,		VP &MBfile },
	{ "bbcutoff",	sf_dbl,		VP &Rbcutoff },
	{ "bbdispfreq",	sf_int,		VP &Mipfreq },
	{ "bbdisplay",	sf_int,		VP &Mipdisplay },
	{ "bbimprove",	sf_dbl,		VP &Rimprove },
	{ "binconv",	sf_int,		VP &Binconv },
	{ "branch",	sf_int,		VP &Branch },
	{ "bs_switch",	sf_int,		VP &Bs_switch },
	{ "changeweight",sf_dbl,	VP &Rchangeweight },
	{ "cholabstol",	sf_dbl,		VP &Rcholabstol },
	{ "choltinytol",sf_dbl,		VP &Rcholtinytol },
	{ "convex_qp",	sf_int,		VP &Convex_qp },
	{ "crash",	sf_int,		VP &Icrashtype },
	{ "cutmult",	sf_dbl,		VP &Rowinc },
	{ "dccutoff",	sf_dbl,		VP &Rdccutoff },
	{ "degscale",	sf_dbl,		VP &Rdegscale },
	{ "densecol",	sf_int,		VP &Idensecol },
	{ "densethr",	sf_dbl,		VP &Rdensethr },
	{ "devexmode",	sf_int,		VP &Idevexmode },
	{ "droprowct",	sf_int,		VP &Idroprowct },
	{ "dspace",	sf_int,		VP &Dspinc },
	{ "dspqfac",	sf_int,		VP &Dspqfac },
	{ "dual",	IK_val,		VP &set_dual },
	{ "dweight",	sf_dbl,		VP &Rdweight },
	{ "endbasis",	sf_bas,		VP &Endbasis },
	{ "fastits",	sf_int,		VP &Ifastits },
	{ "fixvar1",	sf_dbl,		VP &Rfixvar1 },
	{ "fixvar2",	sf_dbl,		VP &Rfixvar2 },
	{ "formntype",	sf_int,		VP &Iformntype },
#ifndef noDEBUG
	{ "fortran",	sf_file,	VP &fort },
#endif
	{ "fracweight",	sf_dbl,		VP &Riweight },
	{ "iadjactype",	sf_int,		VP &Iadjactype },
	{ "idensecol",	sf_int,		VP &Idensecol },
	{ "idevexmode",	sf_int,		VP &Idevexmode },
	{ "idroprowct",	sf_int,		VP &Idroprowct },
	{ "ifastits",	sf_int,		VP &Ifastits },
	{ "iformntype",	sf_int,		VP &Iformntype },
	{ "iheurpass",	sf_int,		VP &Iheurpass },
	{ "iisfind",	sf_int,		VP &Iisfind },
	{ "ilinelen",	sf_int,		VP &Ilinelen },
	{ "imajorits",	sf_int,		VP &Imajorits },
	{ "inullcheck",	sf_int,		VP &Inullcheck },
	{ "ipossbasis",	sf_int,		VP &Ipossbasis },
	{ "ipricetype",	sf_int,		VP &Ipricetype },
	{ "isupertol",	sf_int,		VP &Isupertol },
	{ "iter_inc",	sf_int,		VP &Iter_inc },
	{ "ithreshold",	sf_int,		VP &Ithreshold },
	{ "linelen",	sf_int,		VP &Ilinelen },
	{ "logfreq",	sf_int,		VP &Ilogfreq },
	{ "logfreqb",	sf_int,		VP &Bdisplay },
	{ "majorits",	sf_int,		VP &Imajorits },
	{ "maxfactor",	sf_int,		VP &Imaxfactor },
	{ "maximize",	sf_rknown,	VP &set_max },
	{ "maxiter",	sf_int,		VP &Imaxiter },
	{ "maxiterb",	sf_int,		VP &Imaxiterb },
	{ "maxmin",	sf_minmax,	VP &Rmaxmin },
	{ "maxnodes",	sf_int,		VP &Imaxnodes },
	{ "maxprojns",	sf_int,		VP &Imaxprojns },
	{ "maxsols",	sf_int,		VP &Imaxsols },
	{ "mingap",	sf_dbl,		VP &Mingap },
	{ "minimize",	sf_rknown,	VP &set_min },
	{ "mipstartstatus", sf_int,	VP &Mipststat },
	{ "mpsout",	sf_int,		VP &Mpsout },
	{ "mslv_bfile",	sf_int,		VP &MBfile },
	{ "mslv_file",	sf_int,		VP &MBfile },
	{ "mslv_mfile",	sf_int,		VP &MBfile },
	{ "mufactor",	sf_dbl,		VP &Rmufactor },
	{ "muinit",	sf_dbl,		VP &Rmuinit },
	{ "mulimit",	sf_dbl,		VP &Rmulimit },
	{ "mulinfac",	sf_dbl,		VP &Rmulinfac },
	{ "netalg",	sf_int,		VP &Netalg },
	{ "netinit",	sf_int,		VP &Netinit },
	{ "netprice",	sf_int,		VP &Ipricetype },
	{ "netsamp",	sf_dbl,		VP &Rnetsamp },
	{ "nullcheck",	sf_int,		VP &Inullcheck },
	{ "objno",	sf_objno,	VP &Objno },
	{ "objweight",	sf_dbl,		VP &Robjweight },
	{ "outlev",	sf_int,		VP &Outlev },
	{ "passx0",	sf_int,		VP &Passx0 },
	{ "pdgaptol",	sf_dbl,		VP &Rpdgaptol },
	{ "pdstepmult",	sf_dbl,		VP &Rpdstepmult },
	{ "pertdiag",	sf_dbl,		VP &Rpertdiag },
	{ "plconpri",	sf_int,		VP &Conpri },
	{ "plobjpri",	sf_int,		VP &Objpri },
	{ "possbasis",	sf_int,		VP &Ipossbasis },
	{ "prepassbranch",sf_int,	VP &Isupertol },
	{ "prepassfix",	sf_int,		VP &Ithreshold },
	{ "prepassimprove",sf_dbl,	VP &Rthreshold },
	{ "prepassmax",	sf_int,		VP &Iheurpass },
	{ "prestrat",	sf_int,		VP &Prestrat },
	{ "pretype",	sf_int,		VP &Binsimp },
	{ "pricing",	sf_int,		VP &Pricing },
	{ "primal",	IK_val,		VP &set_primal },
	{ "printcpu",	sf_dbl,		VP &Rprintcpu },
	{ "projtol",	sf_dbl,		VP &Rprojtol },
	{ "pweight",	sf_dbl,		VP &Rpweight },
	{ "rbcutoff",	sf_dbl,		VP &Rbcutoff },
	{ "rchangeweight",sf_dbl,	VP &Rchangeweight },
	{ "rcholabstol",sf_dbl,		VP &Rcholabstol },
	{ "rcholtinytol",sf_dbl,	VP &Rcholtinytol },
	{ "rdccutoff",	sf_dbl,		VP &Rdccutoff },
	{ "rdegscale",	sf_dbl,		VP &Rdegscale },
	{ "rdensethr",	sf_dbl,		VP &Rdensethr },
	{ "rdweight",	sf_dbl,		VP &Rdweight },
	{ "relax",	IK_val,		VP &set_relax },
	{ "rfixvar1",	sf_dbl,		VP &Rfixvar1 },
	{ "rfixvar2",	sf_dbl,		VP &Rfixvar2 },
	{ "rgfactor",	sf_dbl,		VP &Rrgfactor },
	{ "rglimit",	sf_dbl,		VP &Rrglimit },
	{ "rimprove",	sf_dbl,		VP &Rimprove },
	{ "riweight",	sf_dbl,		VP &Riweight },
	{ "rmaxmin",	sf_minmax,	VP &Rmaxmin },
	{ "rmufactor",	sf_dbl,		VP &Rmufactor },
	{ "rmuinit",	sf_dbl,		VP &Rmuinit },
	{ "rmulimit",	sf_dbl,		VP &Rmulimit },
	{ "rmulinfac",	sf_dbl,		VP &Rmulinfac },
	{ "rnetsamp",	sf_dbl,		VP &Rnetsamp },
	{ "robjweight",	sf_dbl,		VP &Robjweight },
	{ "rowinc",	sf_dbl,		VP &Rowinc },
	{ "rpdgaptol",	sf_dbl,		VP &Rpdgaptol },
	{ "rpdstepmult",sf_dbl,		VP &Rpdstepmult },
	{ "rpertdiag",	sf_dbl,		VP &Rpertdiag },
	{ "rprintcpu",	sf_dbl,		VP &Rprintcpu },
	{ "rprojtol",	sf_dbl,		VP &Rprojtol },
	{ "rpweight",	sf_dbl,		VP &Rpweight },
	{ "rrgfactor",	sf_dbl,		VP &Rrgfactor },
	{ "rrglimit",	sf_dbl,		VP &Rrglimit },
	{ "rstepmult",	sf_dbl,		VP &Rstepmult },
	{ "rtarget",	sf_dbl,		VP &Rtarget },
	{ "rthreshold",	sf_dbl,		VP &Rthreshold },
	{ "rtoldinf",	sf_dbl,		VP &Rtoldinf },
	{ "rtolint",	sf_dbl,		VP &Rtolint },
	{ "rtolpinf",	sf_dbl,		VP &Rtolpinf },
	{ "scale",	sf_int,		VP &Scale },
	{ "sensitivity",sf_int,		VP &Sensitivity },
	{ "simplex",	sf_int,		VP &Method },
	{ "simplexinit",sf_int,		VP &Simplexinit },
	{ "simplify",	sf_int,		VP &Simplify },
	{ "sos",	sf_int,		VP &Sos },
	{ "sos2",	sf_int,		VP &Sos2 },
	{ "startbasis",	sf_bas,		VP &Startbasis },
	{ "statuses",	sf_int,		VP &Statuses },
	{ "stderr",	sf_file,	VP &Stderr },
	{ "stepmult",	sf_dbl,		VP &Rstepmult },
	{ "target",	sf_dbl,		VP &Rtarget },
	{ "timing",	sf_int, 	VP &set_timing },
	{ "toldinf",	sf_dbl,		VP &Rtoldinf },
	{ "tolint",	sf_dbl,		VP &Rtolint },
	{ "tolpinf",	sf_dbl,		VP &Rtolpinf },
	{ "trace",	sf_int,		VP &Trace },
	{ "version",	Ver_val,	VP 0, "include version details" },
	{ "wantsol",	WS_val,		0, "write .sol file (without -AMPL)" }
	};

 static Option_Info Oinfo = {"osl", "OSL 3.x", "osl_options", keywds, nkeywds, 0, "",
				0,0,0,0,0, 20020507 };

 static double obj_adj;

 static void
#ifdef KR_headers
nonlin(n, what) int n; char *what;
#else
nonlin(int n, char *what)
#endif
{
	if (n) {
		fprintf(Stderr, "%s contains %s\n", filename, what);
		fexit(4);
		}
	}

 static void
#ifdef KR_headers
rc_check(rc, who) fint rc; char *who;
#else
rc_check(fint rc, char *who)
#endif
{
	if (rc) {
		fflush(stdout);
		fprintf(Stderr, "\nReturn code %ld from %s\n", (long)rc, who);
		fexit(1);
		}
	}

#define chk(x) rc_check(rc, x)
#define objsen rcontrol[2]

 static void
irget(Void)
{
	fint i, rc;

	i = Icontrolsize;
	ekkiget_(&rc, dspace, icontrol, &i);
	chk("ekkiget");

	i = Rcontrolsize;
	ekkrget_(&rc, dspace, rcontrol, &i);
	chk("ekkrget");
	}

 static void
irset(Void)
{
	fint i, rc;

	i = Icontrolsize;
	ekkiset_(&rc, dspace, icontrol, &i);
	chk("ekkiset");

	i = Rcontrolsize;
	ekkrset_(&rc, dspace, rcontrol, &i);
	chk("ekkrset");
	}

 static fint
#ifdef KR_headers
netcheck(a0, ia0, ja0) double *a0; fint *ia0, *ja0;
#else
netcheck(double *a0, fint *ia0, fint *ja0)
#endif
{
	fint *ia, *ia1, *iae, *it, *it0, *ite, *ja, *jae, *jta, *jta0, *jta01,
		*todo, *todo0;
	fint i, j, j1, m1, nz;
	double *a, *a1, t;

	if (!netalg)
		return 0;
	jae = ja0 + N;
	for(a = a0, ja = ja0; ja < jae; ja++)
		switch(ja[1] - ja[0]) {
			case 2:
				if ((t = *a++) != 1. && t != -1.)
					return 0;
				/* no break */
			case 1:
				if ((t = *a++) == 1. || t == -1.)
					break;
			default:
				return 0;
			}
	nz = ja0[N] - 1;
	netsign = (char *)Malloc(m1 = M + 1);
	jta0 = (fint *)Malloc((nz+2*m1)*sizeof(fint));
	todo = todo0 = jta0 + m1;
	it0 = todo + m1;
	memset((char *)jta0, 0, m1*sizeof(fint));
	memset(netsign, 0, m1);
	iae = ia0 + nz;
	jta01 = jta0 - 1;
	for(ia = ia0; ia < iae; ia++)
		jta01[*ia]++;
	for(i = j = 0; ++j <= m1; i = jta01[j] += i);
	--ia0;
	--a0;
	for(i = N; --i >= 0; )
		for(iae = ia0 + *--ja; ia > iae; )
			it0[--jta01[*--ia]] = i;
	for(i = 0; i < M; i++)
		if (!netsign[i]) {
			netsign[i] = 1;
			*todo++ = i;
			while(todo > todo0) {
				jta = jta0 + (j = *--todo);
				it = it0 + jta[0];
				ite = it0 + jta[1];
				while(it < ite) {
					ja  = ja0 + *it++;
					ia  = ia0 + ja[0];
					iae = ia0 + ja[1];
					a   =  a0 + ja[0];
					if ((j1 = *ia - 1) == j) {
						if ((ia1 = ia + 1) >= iae)
							continue;
						a1 = a + 1;
						j1 = *ia1 - 1;
						}
					else {
						a1 = a++;
						ia1 = ia++;
						}
					if (!netsign[j1])
						netsign[*todo++ = j1] =
							*a == *a1
							  ? 3 - netsign[j]
							  : netsign[j];
					else if (*a !=
						  (netsign[j] == netsign[j1]
							? -*a1 : *a1)) {
						free(jta0);
						free(netsign);
						netsign = 0;
						return 0;
						}
					}
				}
			}
	free(jta0);
	return netalg;
	}

 static void
#ifdef KR_headers
msgzap(a, b) fint a,b;
#else
msgzap(fint a, fint b)
#endif
{
	static fint I0 = 0, Im1 = -1;
	fint rc;

	ekkmset_(&rc,dspace,&a,&I0,&Im1,&I0,&I0,&b,&I0);
	chk("ekkmset");
	}

 static void
#ifdef KR_headers
msgzap1(a) fint a;
#else
msgzap1(fint a)
#endif
{ msgzap(a,a); }

 static void
#ifdef KR_headers
Free(xp) real **xp;
#else
Free(real **xp)
#endif
{
	real *x;
	if (x = *xp) {
		*xp = 0;
		free(x);
		}
	}

 static void
#ifdef KR_headers
msgcatch(L) fint L;
#else
msgcatch(fint L)
#endif
{
	static fint I0 = 0, I2 = 2, I256 = 256, Im1 = -1;
	fint rc;

	ekkmset_(&rc,dspace,&L,&I256,&Im1,&I0,&I2,&L,&I0);
	chk("ekkmset");
	}

 static void
#ifdef KR_headers
msgcatch1(L) fint L;
#else
msgcatch1(fint L)
#endif
{
	static fint I0 = 0, I2 = 2, I256 = 256;
	fint rc;

	ekkmset_(&rc,dspace,&L,&I256,&I0,&I0,&I2,&L,&I0);
	chk("ekkmset");
	}

 static void
#ifdef KR_headers
msgcatch2(a, b) fint a, b;
#else
msgcatch2(fint a, fint b)
#endif
{
	static fint I0 = 0, I2 = 2, I256 = 256;
	fint rc;

	ekkmset_(&rc, dspace, &a, &I256, &I256, &I0,
		mipdisplay > 3 ? &I2 : &I0, &b, &I2);
	chk("ekkmset");
	}

  static void OSLCALLBACK
#ifdef KR_headers
ekkmsgu(dsp, Idsp, msgnum, nrvec, rvec, nivec, ivec, ncvec, cvec)
	double *dsp, *rvec;
	fint *Idsp, *msgnum, *nrvec, *nivec, *ivec, *ncvec;
	char *cvec;
#else
ekkmsgu(double *dsp, fint *Idsp, fint *msgnum, fint *nrvec,
	 double *rvec, fint *nivec, fint *ivec, fint *ncvec, char *cvec)
#endif
{
	/* OSL 1.2 passed an extra parameter, cvec_len, giving the length */
	/* of each character string in the cvec array.  OSL 2 omits this */
	/* convenience; we must simply guess that 128 is the correct value; */
	/* in practice, that's what it was with OSL 1.2. */
#define cvec_len 128
	static fint msgcount, ncpu, needhead, nextout;
	static double best;
	static char buf[32];
	double t;
	int j;
	char bufb[16];
	fint mnum = *msgnum;
	static char abandoning[] = "Abandoning QP algorithm:\n";
	static char not_convex[] =
		"the quadratic objective function is not convex because\n";
	static char *msg3063[4] = {
	 "%s%sthe Q matrix has a negative diagonal element for variable %ld.\n",
	 "%s%sthe Q matrix has a zero diagonal element.\n",
		/* Tried adding for var. %ld but ivec[1] can be 0. */
	 "%s%sa principal minor of the Q matrix has determinant < 0.\n",
	 "%s%sa principal minor of the Q matrix has low rank.\n"
		};
	Not_Used(dsp);
	Not_Used(Idsp);
	if (mnum == 259) {	/* get version information */
		memcpy(versbuf, cvec, sizeof(versbuf));
		return;
		}
#ifdef msgDebug
    if (mipdisplay > 2) {
	fint i, n, n_cv, n_iv, n_rv;

	printf("ekkmsgu(%ld,%ld,%ld,%ld)\n", mnum, n_rv = *nrvec,
		n_iv = *nivec, n_cv = *ncvec);
	if (n_iv) {
		n = n_iv /* < 5 ? n_iv : 5 */;
		for(i = 0; i < n; i++)
			printf("\tivec[%ld] = %ld\n", i, ivec[i]);
		}
	if (n_rv) {
		n = n_rv /* < 5 ? n_rv : 5 */;
		for(i = 0; i < n; i++)
			printf("\trvec[%ld] = %g\n", i, rvec[i]);
		}
	if (n_cv) {
		n = n_cv /* < 5 ? n_cv : 5 */;
		for(i = 0; i < n; i++) {
			char *s, *se;
			s = cvec + i*cvec_len;
			se = s + cvec_len;
			while(--se >= s && *se == ' ');
			printf("\tcvec[%ld] = \"%.*s\"\n", i, se+1-s, s);
			}
		}
	return;
	}
#else
	Not_Used(ncvec);
	Not_Used(nivec);
	Not_Used(nrvec);
#endif
	switch(mnum) {
	   case 3:
		intrestore_bug = 1;
	   case 50:
		return;
	   case 85: dualused = 1;
		return;
	   case 3011:
		fprintf(Stderr,
			"OSL bug: \"%ld duplicate matrix elements\".\n\
	*** Try rerunning with a larger dspace. ***\n", ivec[0]);
		fexit(4);
	   case 197:
		if (ncpu++) {
			char *se = cvec + cvec_len;
			while(--se >= cvec && *se == ' ');
			printf(
		" Total CPU time = %g; time in %.*s = %g.\n",
				rvec[0], se - cvec, cvec, rvec[1]);
			}
		return;
	  case 7060:
		ivec[0] = ivec[1];
	  case 3051:
		fprintf(Stderr, "%sdiagonal element %ld is negative.\n",
			abandoning, ivec[0]);
	  case 7003:
	  case 7034:
		fexit(4);
	  case 3063:
		if (!convex_qp || ivec[0] == 3) {
			j = 1 << ivec[0];
			if (!(nonconvex & j)) {
				nonconvex |= j;
				fprintf(Stderr, msg3063[ivec[0]-1],
					"Warning: ", not_convex, (long)ivec[1]);
				if (ivec[0] == 3)
					fprintf(Stderr,
	"The previous Warning message may be erroneous -- an OSL bug.\n");
				}
			return;
			}
		fprintf(Stderr, msg3063[ivec[0]-1], abandoning, not_convex,
			ivec[1]);
		fexit(4);
	  case 7068:
		fprintf(Stderr,
	"Abandoning QP algorithm because of negative curvature (nonconvexity).\n");
		fexit(4);
		}
	if (!msgcount++) {
		needhead = 1;
		if (mipdisplay == 2)
			nextout = mipfreq;
		}
	switch(mnum) {
	  case 153:
		irget();
		g_fmtp(buf, rcontrol[17], 8);
		if (needhead) {
			printf("Iter  Objective       Change\n%4ld  %s\n",
				icontrol[59], buf);
			needhead = 0;
			}
		else {
			if (icontrol[59] < nextout)
				return;
			t = best - rcontrol[17];
			if (negobj1)
				t = -t;
			printf("%4ld  %-14s % .1e\n", icontrol[59], buf, t);
			nextout += bdisplay;
			}
		best = rcontrol[17];
		return;
	  case 204:
		irget();
		if (needhead) {
			needhead = 0;
			printf("Iter  %-14s  %-14s  Gap\n", "Primal", "Dual");
			}
		else if (icontrol[59] < nextout)
			return;
		nextout += bdisplay;
		g_fmtp(buf, rcontrol[17], 8);
		printf("%4ld  %-14s", icontrol[59], buf);
		g_fmtp(buf, rcontrol[35], 8);
		t = rcontrol[17] - rcontrol[35];
		if (negobj1)
			t = -t;
		printf("  %-14s % .2e\n", buf, t);
		return;
	  case 101:
		solving = 1;
		return;
	  case 87:
	  case 102:
	  case 105:
		irget();
		best = rcontrol[22];
		if (solving && best - rcontrol[27] < mingap) {
			mingapped = icontrol[39];
			icontrol[39] = icontrol[36];
			irset();	/* force a stop */
			}
		if (mnum != 105) {
			if (mipdisplay == 1 || icontrol[36] < nextout)
				return;
			nextout += mipfreq;
			break;
			}
		intsols++;
		if (!mipdisplay || outlev >= 3)
			return;
		g_fmtp(buf, negobj1 ? -best : best, 8);
		}
	if (needhead) {
		printf("\n   Nodes  Integer  Simplex\n\
Searched   Soln's    Iters  Best%17s%14s\n","Bound","Gap");
		needhead = 0;
		}
	printf("%8ld%9ld %8ld  ", icontrol[36], intsols, icontrol[3]);
	if (solving && icontrol[41]) {
		g_fmtp(bufb, negobj1 ? -rcontrol[27]: rcontrol[27], 8);
		printf("%-16s%-16s%.1e\n", buf, bufb, best - rcontrol[27]);
		}
	else
		printf("%s\n", buf);
	}

 static void
getversion(VOID)
{
	static fint I0=0, I1=1, I2=2, I258=258, I259=259, I6=6, I7=7, Im1=-1;
	fint rc;
	char *s, *t;

	ekkrgcb(&I7, (Cbfunc*)ekkmsgu);
	dspace = (double*)Malloc(12000*sizeof(double));
	memset(dspace, 0, 12000*sizeof(double));
#ifdef OSL_V3
	/* Bypass bug with ekkmset: the first call always complains of an invalid */
	/* message number and chatters about it on stdout. */
	WantPrintf = WantPrintf1;
	ekkmset_(&rc, dspace, &I6, &I0, &Im1, &I0, &I2, &I259, &I2);
	WantPrintf = 1;
#endif
	ekkmset_(&rc, dspace, &I6, &I0, &Im1, &I0, &I2, &I259, &I2);
	chk("ekkmset");
	ekkinit_(&rc, dspace);
	chk("ekkinit");
	ekkmset_(&rc, dspace, &I6, &I0, &Im1, &I0, &I1, &I258, &I2);
	chk("ekkmset");
	free(dspace);
	t = versbuf + sizeof(versbuf) - 2;
	*t = ' ';
	for(s = versbuf; *s != ' '; s++);
	if (rc = s - versbuf) {
		*s = 0;
		Oinfo.bsname = (char*)Malloc(rc + 5);
		strcpy(Oinfo.bsname, "OSL ");
		strcpy(Oinfo.bsname+4, versbuf);
		if (t > s) {
			while(*--t == ' ');
			if (t <= s)
				return;
			*++t = '\n';
			*++t = 0;
			*s = ' ';
			s = Oinfo.version = (char*)Malloc(t - versbuf + 10);
			strcpy(s, "AMPL/OSL ");
			strcpy(s+9, versbuf);
			}
		}
	}

 static int
ircset(VOID)
{
	int i, rv;


	for(i = rv = 0; i < Icontrolsize; i++)
		if (icset[i+1]) {
			icontrol[i] = ic1[i];
			++rv;
			}
	for(i = 0; i < Rcontrolsize; i++)
		if (rcset[i+1]) {
			rcontrol[i] = rc1[i];
			++rv;
			}
	return rv;
	}

 static int
#ifdef KR_headers
compar(a, b, v) char *a, *b, *v;
#else
compar(const void *a, const void *b, void *v)
#endif
{
	int *p = (int*)v;
	return p[*(int*)a] - p[*(int*)b];
	}

 static fint
#ifdef KR_headers
priorities(nint_p, intnums_p, type_p,
	priority_p, ntotinfo_p, setindex_p,
	sets_p, dlnpc_p, dunpc_p)
	fint *nint_p, **intnums_p, **type_p, **priority_p, *ntotinfo_p;
	fint **setindex_p, **sets_p; double **dlnpc_p, **dunpc_p;
#else
priorities(fint *nint_p, fint **intnums_p, fint **type_p,
	fint **priority_p, fint *ntotinfo_p, fint **setindex_p,
	fint **sets_p, double **dlnpc_p, double **dunpc_p)
#endif
{
	SufDesc *dp, *dd;
	char *sostype;
	int i, inlv, isets, isos, j, k, n, nint, nlv, nreal, nsos, nsosnz;
	int copri[2], *dir, *num, *p, *pri, *q, *qi;
	int *sosbeg, *sosind, *sospri, *start;
	fint *intnums, ni, nsets, ntotinfo, pmax, *priority;
	fint *s0, *setindex, *sets, *type;
	double *dlnpc, down, *dunpc, up;
	real *sosref;

	nlv = nlvo;
	inlv = nlv - nlvoi;
	nreal = n_var - (nbv + niv);
	i = ASL_suf_sos_explict_free;
	if (!sos)
		i |= ASL_suf_sos_ignore_sosno;
	if (!sos2)
		i |= ASL_suf_sos_ignore_amplsos;
	copri[0] = (int)objpri;
	copri[1] = (int)conpri;
	isos = 0;	/* number of integers in SOS sets */
	if (nsos = suf_sos(i, &nsosnz, &sostype, &sospri, copri,
			&sosbeg, &sosind, &sosref)) {
		nreal = n_var - (nbv + niv);
		n = sosbeg[nsos] - 1;
		for(i = 0; i < n; )
			if ((j = sosind[i++]) > nreal
			 || j <= nlv && j > inlv)
				isos++;
		}
	isets = 0; /* sets of integers with the same priority */
	/* suf_sos may have changed nbv */
	if (nint = nbv + niv + nlvoi) {
		dp = suf_get("priority", ASL_Sufkind_var);
		if (p = dp->u.i) {
 have_p:
			q = qi = (int*)Malloc(nint*sizeof(int));
			nreal = n_var - (nbv + niv);
			n = nlv;
			for(i = inlv; i < n; i++)
				*qi++ = i;
			n = n_var;
			for(i = nreal; i < n; i++)
				*qi++ = i;
			n = 0;
			for(i = 0; i < nint; i++) {
				if ((j = p[i]) < 0) {
					n++;
					p[i] = 0;
					}
				else if (j > MAX_PRIORITY) {
					n++;
					p[i] = MAX_PRIORITY;
					}
				}
			if (n)
				fprintf(Stderr,
				 "%d priorit%s projected onto [0, %ld].\n",
					n, n > 1 ? "ies" : "y",
					(long)MAX_PRIORITY);
			if (!isets) {
				if (k > 1)
					qsortv(q, k, sizeof(int), compar, p);
				n = -1;
				for(i = 0; i < k; i++) {
					if (n != (j = p[q[i]])) {
						n = j;
						isets++;
						}
					}
				}
			}
		else {
			p = (int*)M1zapalloc(n_var*sizeof(int));
			n = mip_pri(&start, &num, &pri, (long)MAX_PRIORITY);
			if (n > 0) {
				isets = n;
				while(n-- > 0) {
					i = start[n];
					j = i + num[n];
					k = pri[n];
					while(i < j)
						p[i++] = k;
					}
				}
			else
				isets = 1;
			goto have_p;
			}
		}
	if (!(nsets = nsos + isets))
		return 0;
	*nint_p = ntotinfo = nint;
	ni = 0;
	if (nsos) {
		ntotinfo += nsosnz;
		*nint_p = ni = ntotinfo - isos;
		if (!isos)
			ni = 0;
		}
	*ntotinfo_p = ntotinfo;
	*dlnpc_p = dlnpc = (double*)
		M1alloc((ni+3*nsets+1+ntotinfo)*sizeof(fint)
			+ 2*sizeof(double)*ntotinfo);
	*dunpc_p = dunpc = dlnpc + ntotinfo;
	*intnums_p = intnums = (fint*)(dunpc + ntotinfo);
	*sets_p = sets = s0 = intnums + ni;
	*type_p = type = sets + ntotinfo;
	*priority_p = priority = type + nsets;
	*setindex_p = setindex = priority + nsets;
	if (ni) {
		n = nlv;
		for(i = inlv; i++ < n; )
			*intnums++ = i;
		n = n_var;
		for(i = nreal; i++ < n; )
			*intnums++ = i;
		n = sosbeg[nsos] - 1;
		for(i = 0; i < n; )
			if ((j = sosind[i++]) <= inlv
			 || j > nlv && j <= nreal)
				*intnums++ = j + 1;
		}
	down = up = .001;
	if (nint) {
		dd = suf_get("direction", ASL_Sufkind_var);
		dir = dd->u.i;
		n = -1;
		for(i = 0; i < nint; i++) {
			k = q[i];
			*sets++ = k + 1;
			if (n != (j = p[k = q[i]])) {
				*priority++ = n = j;
				*setindex++ = sets - s0;
				*type++ = 4;
				}
			if (dir) {
				down = up = .001;
				if (k = dir[k])
					if (k > 0)
						down = 1.001;
					else
						up = 1.001;
				}
			*dlnpc++ = down;
			*dunpc++ = up;
			}
		free(q);
		}
	j = n = 0;
	for(i = 0; i < nsos; ) {
		*setindex++ = sets - s0 + 1;
		*type++ = sostype[i] - '0';
		if ((k = sospri[i]) < 0) {
			k = 0;
			n++;
			}
		else if (k > MAX_PRIORITY) {
			k = MAX_PRIORITY;
			n++;
			}
		*priority++ = k;
		k = sosbeg[++i] - 1;
		while(j < k) {
			*sets++ = sosind[j];
			*dlnpc++ = .001;
			*dunpc++ = sosref[j++];
			}
		}
	*setindex = sets - s0 + 1;
	pmax = 1000;
	for(s0 = *priority_p; s0 < priority; s0++)
		if (pmax < *s0)
			pmax = *s0;
	pmax++;
	for(s0 = *priority_p; s0 < priority; s0++)
		*s0 = pmax - *s0;
	return nsets;
	}

 static void
#ifdef KR_headers
badstat(what, i, j) char *what; int i; int j;
#else
badstat(char *what, int i, int j)
#endif
{
	fprintf(Stderr, "OSL driver: %s[%d] = %d\n", what, i, j);
	}

 static void
#ifdef KR_headers
suppressed(bad, what) int bad; char *what;
#else
suppressed(int bad, char *what)
#endif
{
	fprintf(Stderr,
		"OSL driver: %d messages about bad %s values suppressed.\n",
			bad-1, what);
	}

#define BASIC	0x80000000
#define UPPER	0x40000000
#define LOWER	0x20000000
#define FREE	0x60000000
#define FIXED	0

 static void
#ifdef KR_headers
stat_map_in(os, stat, n, what) ULong *os; int *stat, n; char *what;
#else
stat_map_in(ULong *os, int *stat, int n, char *what)
#endif
{
	int bad, i, i1, j, j1;
	static ULong map[7] = { FREE, BASIC, FREE, LOWER, UPPER, FIXED, FREE };

	for(i = bad = 0; i < n; i++)
		if ((j = stat[i]) >= 0 && j <= 6)
			os[i] = map[j];
		else {
			os[i] = FREE;
			i1 = i;
			j1 = j;
			if (!bad++)
				badstat(what, i, j);
			}
	if (bad > 1)
		if (bad == 2)
			badstat(what, i1, j1);
		else
			suppressed(bad, what);
	}

 static void
#ifdef KR_headers
stat_map_out(os, stat, n) ULong *os; int *stat, n;
#else
stat_map_out(ULong *os, int *stat, int n)
#endif
{
	ULong u;
	int i;
	static int map[4] = { 5, 3, 4, 6 };

	for(i = 0; i < n; i++)
		stat[i] = ((u = os[i]) & BASIC) ? 1 : map[(u >> 29) & 3];
	}

 static void
#ifdef KR_headers
get_statuses(dspace, sti) double *dspace; Statinfo *sti;
#else
get_statuses(double *dspace, Statinfo *sti)
#endif
{
	SufDesc *sd;
	ULong *cs, *rs;
	fint cv[10], rc;
	int *cstat, havestats, i, j, m, n, *rstat;
	real *L, *U;
	static fint I10 = 10;

	if (!(statuses & 1))
		return;
	havestats = 0;
	sd = sti->csd;
	n = n_var;
	if (sd->kind & ASL_Sufkind_input)
		havestats = 1;
	sd = sti->rsd;
	m = n_con;
	if (sd->kind & ASL_Sufkind_input)
		havestats = 1;
	if (!havestats || !mipststat && niv + nbv)
		return;
	ekknget_(&rc, dspace, cv, &I10);
	chk("ekknget");
	pricing = 0;
	cs = (ULong*)dspace + (cv[9] - 1);
	rs = (ULong*)dspace + (cv[4] - 1);
	cstat = sti->cstat;
	rstat = sti->rstat;
	if (primaldual < 2) {
		if (sti->objadj) {
			sti->cstat[n++] = 1;
			rs[m] = FREE;
			}
		stat_map_in(cs, cstat, n, "cstat");
		stat_map_in(rs, rstat, m, "rstat");
		return;
		}
	m = n_con;
	n = n_var;
	L = LUrhs;
	U = Urhsx;
	for(i = 0; i < m; i++) {
		if (L[i] <= osl_neginfinity)
			cs[i] = rstat[i] == 1 ? LOWER : BASIC;
		else if (U[i] >= osl_infinity)
			cs[i] = rstat[i] == 1 ? UPPER : BASIC;
		else if (L[i] == U[i])
			cs[i] = BASIC;
		else
			cs[i] = rstat[i] == 3 ? BASIC : LOWER;
		}
	if (sti->objadj) {
		cs[i++] = BASIC;
		rs[n] = FREE;
		}
	if (nranges) {
		for(j = 0; j < m; j++)
			if (L[j] < U[j]
			 && L[j] > osl_neginfinity
			 && U[j] < osl_infinity)
				cs[i++] = rstat[j] == 4
					? BASIC : UPPER;
		}
	L = LUv;
	U = Uvx;
	for(j = 0; j < n; j++) {
		rs[j] = LOWER;
		if (U[j] < osl_infinity) {
			if (U[j])
				cs[i++] = cstat[j] == 4
					? BASIC : LOWER;
			else
				rs[j] = cstat[j] == 4
					? BASIC : LOWER;
			if (L[j] > osl_neginfinity)
				goto finite_lb;
			}
		else if (L[j] > osl_neginfinity) {
 finite_lb:
			if (L[j])
				cs[i++] = cstat[j] == 3
					? BASIC : LOWER;
			else
				rs[j] = cstat[j] == 3 ?
					BASIC : UPPER;
			}
		else if (cstat[j] != 1)
			rs[j] = BASIC;
		}
	}

 static void
cud_cadjust(Void)
{
	int i, m;
	real *L, *U;

	m = n_con;
	L = LUrhs;
	U = Urhsx;
	for(i = 0; i < m; i++)
		if (L[i] <= osl_neginfinity)
			L[i] = U[i];
	}

 static char *
#ifdef KR_headers
send_statuses(dspace, sti) double *dspace; Statinfo *sti;
#else
send_statuses(double *dspace, Statinfo *sti)
#endif
{
	ULong *cs, *rs;
	fint cv[48], rc;
	int i, j, k, m, n;
	int *cstat, *rstat;
	real *L, *U;
	static fint I1 = 1, I10 = 10, I2 = 2, I32 = 32, I48 = 48;

	if (!(Oinfo.wantsol & 1) && !amplflag)
		return 0;
	if (!havebasis || !(statuses & 2)) {
		sti->csd->kind &= ~ASL_Sufkind_output;
		sti->rsd->kind &= ~ASL_Sufkind_output;
		return 0;
		}
	ekknget_(&rc, dspace, cv, &I10);
	chk("ekknget");
	cs = (ULong*)dspace + (cv[9] - 1);
	rs = (ULong*)dspace + (cv[4] - 1);
	cstat = sti->cstat;
	rstat = sti->rstat;
	memset(cstat, 0, n_var*sizeof(int));
	memset(rstat, 0, n_con*sizeof(int));
	if (primaldual < 2) {
		stat_map_out(cs, cstat, n_var);
		stat_map_out(rs, rstat, n_con);
		equ_adjust(cstat, rstat);
		if (sensitivity && solve_result_num < 200) {
			ekksobj_(&rc, dspace);
			if (rc)
				return "Sensitivities not";
			ekknget_(&rc, dspace, cv, &I32);
			chk("ekknget");
			L = dspace + cv[31] - 1;
			U = dspace + cv[30] - 1;
			suf_rput("down", ASL_Sufkind_var, L);
			suf_rput("up", ASL_Sufkind_var, U);
			suf_rput("current", ASL_Sufkind_var, sti->c);
			ekksbnd_(&rc, dspace, &I1, &I2);
			if (rc)
				return "Variable sensitivity constraints";
			ekknget_(&rc, dspace, cv, &I48);
			chk("ekknget");
			L = dspace + cv[47] - 1;
			U = dspace + cv[46] - 1;
			cud_cadjust();
			suf_rput("down", ASL_Sufkind_con, L);
			suf_rput("up", ASL_Sufkind_con, U);
			suf_rput("current", ASL_Sufkind_con, LUrhs);
			return "Sensitivities";
			}
		return 0;
		}
	m = n_con;
	n = n_var;
	L = LUrhs;
	U = Urhsx;
	for(i = 0; i < m; i++) {
		if (L[i] <= osl_neginfinity)
			rstat[i] = cs[i] & LOWER ? 1 : 4;
		else if (U[i] >= osl_infinity)
			rstat[i] = cs[i] & UPPER ? 1 : 3;
		else if (L[i] == U[i])
			rstat[i] = 5;
		else
			rstat[i] = cs[i] & BASIC ? 3 : 1;
		}
	i += sti->objadj;
	if (nranges) {
		for(j = 0; j < m; j++)
			if (L[j] < U[j]
			 && L[j] > osl_neginfinity
			 && U[j] < osl_infinity
			 && cs[i++] & BASIC)
				rstat[j] = 4;
		}
	L = LUv;
	U = Uvx;
	for(j = 0; j < n; j++) {
		k = 1;
		if (U[j] < osl_infinity) {
			if (U[j]) {
				if (cs[i++] & BASIC)
					k = 4;
				}
			else {
				if (rs[j] & BASIC)
					k = 4;
				}
			if (L[j] > osl_neginfinity)
				goto finite_lb;
			}
		else if (L[j] > osl_neginfinity) {
 finite_lb:
			if (L[j]) {
				if (cs[i++] & BASIC)
					k = 3;
				}
			else {
				if (rs[j] & BASIC)
					k = 3;
				}
			}
		else {
			if (rs[j] & BASIC)
				k = 6;
			}
		cstat[j] = k;
		}
	equ_adjust(cstat, rstat);
	return 0;
	}

 void
#ifdef KR_headers
amplin(stub, av, sti) char *stub, **av; Statinfo *sti;
#else
amplin(char *stub, char **av, Statinfo *sti)
#endif
{
	FILE *nl;
	char *ns1;
	double *a, *a0, *a1, *a2, *c, *c1, *clb, *cub, *cub1;
	double *delsq, *dlnpc, *dunpc, *lb, *lb1, *lu, *ub, *ub1, *x0, **xpi0;
	double os, t;
	fint *colq, *ia, *ia0, *ia1, *ia2, *intnums, *ja, *ja0, *ja1;
	fint *priority, *rowq, *setindex, *sets, *type;
	fint L[7], MXROW, MXCOL;
	fint i, j, k, m, m0;
	fint n, nblock, nextra, nint, nmodels, nnz;
	fint nsets, ntotinfo, nwords, nxpi0, nz, nzextra, rc;
	int neednames;
	ograd *og;
	static fint I0 = 0, I1 = 1, I2 = 2, I256 = 256, I2999 = 2999;
	static fint I48 = 48, I50 = 50, I7 = 7, I8 = 8;

	nl = jacdim(stub, &M, &N, &NO, &NZ, &MXROW, &MXCOL, (fint)strlen(stub));

	if (N <= 0) {
		fprintf(Stderr, "%s has no variables\n", filename);
		fexit(4);
		}
	*stub_end = 0;
	sti->objadj = 0;

	/* Allow space for adding a row to adjust the objective value. */
	m = M + 2;
	n = N + 1;
	nz = NZ + 1;

	/* Allow space for adding a row to network problems */
	/* to make ekknslv happy: add one super source/sink */
	/* if presolve omitted some sources and sinks.	    */

	if (NZ > 2*N)
		netalg = 0;
	else {
		nz += 2*n - nz;
		m++;
		}

	sti->cstat = (int*)M1zapalloc((m+n)*sizeof(int));
	sti->rstat = sti->cstat + n;
	sti->csd = suf_iput("sstatus", ASL_Sufkind_var, sti->cstat);
	sti->rsd = suf_iput("sstatus", ASL_Sufkind_con, sti->rstat);

	LUrhs = clb = (double *)M1alloc((2*m + 3*n)*sizeof(double));
	Urhsx = cub = LUrhs + m;
	c = cub + m;
	LUv = lb = c + n;
	Uvx = ub = LUv + n;

	A_vals = a = (double *)Malloc(nz*sizeof(double)
					+ (nz + n + 1)*sizeof(fint));
	ia = (fint *)(A_vals + nz);
	ja = ia + nz;
	A_rownos = (int *)ia;
	A_colstarts = (int *)ja;
	Fortran = 1;
	want_xpi0 = 3;
	want_deriv = 0;

	qp_read(nl,0);

	lu = lb + N;
	ub += N;
	while(lu > lb) {
		if (*--lu <= negInfinity)
			*lu = osl_neginfinity;
		if (*--ub >= Infinity)
			*ub = osl_infinity;
		}
	lu = clb + M;
	cub += M;
	while(lu > clb) {
		if (*--lu <= negInfinity)
			*lu = osl_neginfinity;
		if (*--cub >= Infinity)
			*cub = osl_infinity;
		}

	if (M > 0)
		rowinc = 100. / M;
	if (rowinc < 0.25)
		rowinc = 0.25;

	if (getopts(av, &Oinfo))
		exit(1);
	if (need_nl && Nchatter) {
		need_nl = 0;
		printf("\n");
		}
	nelq = qpcheck(&rowq, &colq, &delsq);

	nsets = 0;
	nint = nbv + niv + nlvoi;
	if (relax)
		relax = nint;
	else {
		nsets = priorities(&nint, &intnums, &type, &priority,
				&ntotinfo, &setindex, &sets, &dlnpc, &dunpc);
		/* suf_sos() may have changed the following... */
		M = n_con;
		N = n_var;
		NZ = nzc;
		}

	nmodels = 1;
	if (sti->nsets = nsets) {
		primaldual = 1;
		netalg = 0;
		if (binconv || barrier >= 0)
			nmodels = 2;
		}
	if (nelq) {
		primaldual = 1;
		netalg = 0;
		simplify = -1;
		}

	if (primaldual == 1) {
		xpi0 = &X0;
		nxpi0 = N;
		}
	else {
		xpi0 = &pi0;
		nxpi0 = M;
		}

	binig = neednames = 0;
	if (nint) {
		binig = nbv;
		if (nlvoi) {
			lb = LUv;
			ub = Uvx;
			j = nlvo;
			for(i = j - nlvoi; i < j; i++)
				if (lb[i] == 0. && ub[i] == 1.)
					binig++;
			}
		if (binconv)
			if (binig == nint)
				binconv = 0;
			else
				neednames = 1;
		if (X0)
			neednames = 1;
		else
			binig = 0;
		}

	if (startbasis)
		neednames = 1;
	else if (!crash)
		switch(simplexinit) {
			case 0: passx0 = 0; break;
			case 2: pricing = 3; /* no break */
			case 1: passx0 = 1;
			}

	/* check for pure network */

	if (barrier >= 0)
		netalg = 0;
	else if (netalg)
		netalg = primaldual == 2 ? 0 : netcheck(a, ia, ja);

	nwords = 10000 + nmodels*1000 + 10*(M + N + NZ + 3)
		+ dspinc + nelq*dspqfac;
	/* OSL 1.2 sometimes faulted without nelq*dspqfac. */
	m = M;
	n = N;
	if (nsets) {
		nwords += 40000 + 10*nmodels*rowinc*(M + NZ)
				+ 3*(nmodels*nint + nsets);
		}
	if (dspinc)
		printf("Allocating %ld double words for dspace.\n", (long)nwords);
	dspace = (double *)Malloc(nwords*sizeof(double));
	memset((char *)dspace, 0, nwords*sizeof(double));

	msgzap(6L, 259L);
	ekkinit_(&rc, dspace);
	chk("ekkinit");
#ifdef msgDebug
	if (mipdisplay > 3)
		msgcatch2(1L, 8999L);
	else
#endif
		{
		msgcatch2(7L,47L);
		msgcatch2(51L, 258L);
		}
	msgzap1(81L);
	msgzap1(83L);
	msgcatch(85L);
	if (outlev & 1) {
		ekkmset_(&rc,dspace,&I1,&I0,&I0,&I0,&I0,&I2999,&I1);
		chk("ekkmset");
		}
	if (outlev < 3) {
		msgzap1(3001L);	/* msg about unbounded problem.	*/
		msgzap1(3048L);	/*  "	 "	 "	  "	*/
		msgzap1(3069L); /* "All 0-1 variables are satisfied" */
		msgzap(280L, 287L);
		msgzap1(3060L);
		if (outlev <= 0)
			msgzap(1L, 243L);
		else {
			msgzap1(  2L);
			msgzap1(  6L);
			msgzap(  16L,   20L);
			msgzap1( 23L);
			msgzap(  31L,	35L);
			msgzap1( 37L);
			msgzap(  75L,   78L);
			msgzap1( 82L);
			msgzap1( 84L);
			msgzap(  87L,  171L);
			msgzap1(174L);
			msgzap( 183L,  196L);
			msgzap( 198L,  243L);
			msgzap(3070L, 3070L);
			if (!icset[Ilogfreq.isub]) {
				msgzap1(   1L);
				msgzap1(  38L);
				msgzap1(  57L);
				msgzap1(3000L);
				}
			}
		}

	WantPrintf = 0; /* suppress license chatter; use "osl -v" to see it */
	ekkdsca_(&rc, dspace, &nwords, &nmodels);
	WantPrintf = 1;
	chk("ekkdsca");
	msgcatch1(7003L);
	nblock = nelq ? 5 : 1;
	ekkdscm_(&rc, dspace, &I1, &nblock);
	chk("ekkdscm");

	if (nint && rowinc > 0) {
		icset[9] = 1;
		if ((i = M) < 12)
			i = 12;
		ic1[8] = i + rowinc*i;
		}

	irget();
	icontrol[24] = 1;	/* msg no's at right */
	/* song and dance because we can only call ekkdsca once */
	if (prestrat != 1 || branch) {
		ic1[51] = prestrat | branch >> 2;
		icset[52] = 1;
		}
	nircset = ircset();
	if (icset[9]) {
		--nircset;
		icset[9] = 0;
		/* This otherwise causes trouble at the ircset() call below. */
		}

	if (trace) {
		ekkmset_(&rc,dspace,&I48,&I0,&I256,&I0,&I0,&I50,&trace);
		chk("ekkmset");
		}
	else
		msgzap1(53L);

	nz = NZ;
	memset(c, 0, N*sizeof(double));
	if (obj_no >= 0 && obj_no < n_obj) {
		for(og = Ograd[obj_no]; og; og = og->next)
			c[og->varno] = og->coef;
		if (need_minmax)
			objsen = objtype[obj_no] ? -1 : 1;
		if (nelq && objsen <= 0.)
			if (objsen < 0.) {
				objsen = -objsen;
				negobj = 1;
				fprintf(Stderr,
			"Maximizing a quadratic objective function; OSL\n\
will minimize the negative of the objective instead.\n");
				for(og = Ograd[obj_no]; og; og = og->next)
					c[og->varno] = -og->coef;
				for(i = 0; i < nelq; i++)
					delsq[i] = -delsq[i];
				}
			else
				nelq = 0;
		obj_adj = objconst(obj_no);
		if (negobj)
			obj_adj = -obj_adj;
		if (obj_adj || !m) {
			sti->objadj = 1;
			clb[m] = cub[m] = obj_adj;
			lb[n] = osl_neginfinity;
			ub[n] = osl_infinity;
			c[n++] = 1;
			ia[nz] = ++m;
			a[nz++] = 1;
			ja[n] = nz + 1;
			}
		}
	else
		objsen = 0;
	negobj1 = objsen < 0.;	/* objsen may be wrong in ekkmsgu_ */
	if (primaldual == 2)
		rc1[2] = rcontrol[2] = -rcontrol[2];
	if ((rc1[2] = objsen) != 1.)
		rcset[3] = 1;	/* for fortran */
	if (rc1[44])
		/* kludge around OSL rprintcpu initialization bug */
		msgcatch(197L);
	irset();

	nonlin(nlc, "nonlinear constraints");
	nonlin(plterms, "piecewise-linear terms");

	if (netalg) {
		ns1 = netsign - 1;
		if (i = 2*n - nz) {
			clb[m] = osl_neginfinity;
			netsign[m] = 2;
			cub[m++] = osl_infinity;
			a1 = a + nz;
			ia1 = ia + nz;
			nz += i;
			a2 = a + nz;
			ia2 = ia + nz;
			ja1 = ja + n;
			for(;;) {
				j = *ja1 - 1;
				*ja1 += i;
				*--a2 = *--a1;
				*--ia2 = *--ia1;
				if (*--ja1 == j) {
					*--a2 = ns1[*ia1] == 2 ? -*a1 : *a1;
					*--ia2 = m;
					if (!--i)
						break;
					}
				else {
					*--a2 = *--a1;
					*--ia2 = *--ia1;
					}
				}
			}
		for(i = 0; i < m; i++)
			if (netsign[i] == 2) {
				t = clb[i];
				clb[i] = -cub[i];
				cub[i] = -t;
				}
		a2 = a + nz;
		ia2 = ia + nz;
		while(a2 > a) {
			--a2;
			if (ns1[*--ia2] == 2)
				*a2 = -*a2;
			}
		}

	if (nint) {
		if (relax) {
			printf("ignoring integrality of %ld variable%s\n",
				(long)nint, nint > 1 ? "s" : "");
			}
		}
	if (primaldual == 2) {
		lu = lb + N;
		ub += N;
		nextra = 0;
		while(lu > lb) {
			if (*--lu > osl_neginfinity && *lu)
				nextra++;
			if (*--ub < osl_infinity && *ub)
				nextra++;
			}
		nzextra = nextra;
		ja0 = ja;
		ja = (fint *)Malloc(m*sizeof(fint));
		memset((char *)ja, 0, m*sizeof(fint));
		ia1 = ia + nz;
		while(ia1 > ia)
			ja[*--ia1 - 1]++;
		lu = clb + M;
		cub += M;
		while(lu > clb) {
			--cub;
			if (*--lu > osl_neginfinity
			 && *cub < osl_infinity
			 && *lu < *cub) {
				nextra++;
				i = lu - clb;
				nzextra += ja[i];
				}
			}
		free((char *)ja);
		m0 = m;
		i = m + nextra;
		m = n;
		n = i;
		x0 = *xpi0;
#if X0_IN_DSPACE
		if (x0) {
			if (nranges) {
				i = m0 + nranges;
				a0 = (real*)M1zapalloc(i*sizeof(real));
				memcpy(a0, x0, m0*sizeof(real));
				/* x0 was allocated by M1zapalloc */
				*xpi0 = x0 = a0;
				}
			if (t = -objsen) {
				if (t > 0.)
					t = -t;
				if (negobj)
					t = -t;
				if (t != 1.)
					for(i = 0; i < nxpi0; i++)
						x0[i] *= t;
				}
			}
#endif
		if (os = rcontrol[2])
			os = os > 0 ? 1. : -1.;
		nnz = nz + nextra + nzextra;
		a0 = a;
		a = (double *)Malloc((nnz + m + 3*n)*sizeof(double)
					+ (nnz + n + 1)*sizeof(fint));
		cub1 = a + nnz;
		c1 = cub1 + m;
		lb1 = c1 + n;
		ub1 = lb1 + n;
		ia0 = ia;
		ia = (fint *)(ub1 + n);
		ja = ia + nnz;
		memset((char *)ja, 0, m0*sizeof(fint));
		ia1 = ia0 + nz;
		while(ia1 > ia0)
			ja[*--ia1 - 1]++;
		for(i = j = 0; i < m0; i++)
			j = ja[i] += j;
		ja[m0] = nz;
		a1 = a0 + nz;
		ia1 = ia0 + nz;
		j = m + 1;	/* m = old n */
		while(--j >= 1) {
			ia2 = ia0 + (ja0[j-1] - 1);
			while(ia1 > ia2) {
				ia[i = --ja[*--ia1 - 1]] = j;
				a[i] = *--a1;
				}
			}
		n = m0;
		for(i = 0; i < m0; i++) {
			lb1[i] = 0.;
			ub1[i] = osl_infinity;
			if (clb[i] <= osl_neginfinity)
				c1[i] = cub[i];
			else if (cub[i] >= osl_infinity) {
				c1[i] = clb[i];
				lb1[i] = osl_neginfinity;
				ub1[i] = 0.;
				}
			else {
				c1[i] = cub[i];
				if (clb[i] == cub[i])
					lb1[i] = osl_neginfinity;
				else {
					lb1[n] = osl_neginfinity;
					ub1[n] = 0.;
					c1[n] = clb[i];
					if (x0)
						if (x0[i] >= 0.)
							x0[n] = 0.;
						else {
							x0[n] = x0[i];
							x0[i] = 0.;
							}
					ja[n++] = nz;
					j = ja[i];
					k = ja[i+1] - j;
					memcpy((char *)(a+nz), (char *)(a+j),
						k*sizeof(double));
					memcpy((char *)(ia+nz), (char *)(ia+j),
						k*sizeof(fint));
					nz += k;
					}
				}
			}
		clb = c;
		cub = cub1;
		c = c1;
		for(i = 0; i < m; i++) {
			cub[i] = clb[i] *= os;
			lu = ub + i;
			if (*lu < osl_infinity)
				if (*lu) {
					lb1[n] = 0.;
					ub1[n] = osl_infinity;
					c1[n] = *lu;
					ja[n++] = nz;
					a[nz] = 1.;
					ia[nz++] = i + 1;
					}
				else
					clb[i] = osl_neginfinity;
			lu = lb + i;
			if (*lu > osl_neginfinity)
				if (*lu) {
					lb1[n] = 0.;
					ub1[n] = osl_infinity;
					c1[n] = -*lu;
					ja[n++] = nz;
					a[nz] = -1.;
					ia[nz++] = i + 1;
					}
				else
					cub[i] = osl_infinity;
			}
		ja[n] = nz;
		for(i = 0; i <= n; i++)
			ja[i]++;
		if (os < 0)
			for(i = 0; i < n; i++)
				c[i] = -c[i];
		lb = lb1;
		ub = ub1;
		free((char *)A_vals);
		}
#ifndef noDEBUG
	if (fort) {
		fprintf(fort, "\tprogram bugsho\n\n");
		fprintf(fort,
			"\tinteger i, L(7), xoff\n\tdouble precision f\n\n");
		fprintf(fort, "\tdouble precision dspace(%ld), rctl(%d)\n",
			(long)nwords, Rcontrolsize);
		fprintf(fort, "\tinteger ia(%ld), ictl(%d), ja(%ld)\n",
			(long)nz, Icontrolsize, n+1);
		fprintf(fort, "\tinteger m, n, nwords, nz, rc\n");
		fprintf(fort,
		  "\tdouble precision a(%ld), c(%ld), clb(%ld), cub(%ld)\n",
				(long)nz, (long)n, (long)m, (long)m);
		fprintf(fort, "\tdouble precision lb(%ld), ub(%ld)\n",
			(long)n, (long)n);
		if (nint) {
			fprintf(fort, "\tinteger nint, nsets\n");
			fprintf(fort,
	"\tinteger intnms(%ld), prirty(%ld), setind(%ld), type(%ld)\n",
				(long)nint, (long)nsets, (long)nsets+1,
				(long)nsets);
			fprintf(fort,
				"\tdouble precision dlnpc(%ld), dunpc(%ld)\n",
				(long)nint, (long)nint);
			}
#if X0_IN_DSPACE
		if (*xpi0)
			fprintf(fort, "\tdouble precision x0(%ld)\n", (long)n);
#endif
		if (nelq)
			fprintf(fort, "\tinteger nelq, rowq(%ld), colq(%ld)\n\
	double precision delsq(%ld)\n",
				(long)nelq, (long)n+1, (long)nelq);
		fprintf(fort,"\tdata m/%ld/, n/%ld/, nwords/%ld/, nz/%ld/\n",
			(long)m, (long)n, (long)nwords, (long)nz);
		for(i = 0; i <= n; i++)
			fprintf(fort, "\tdata ja(%ld)/%ld/\n", (long)i+1,
				(long)ja[i]);
		for(i = 0; i < nz; i++)
			fprintf(fort, "\tdata ia(%ld)/%ld/, a(%ld)/%.16g/\n",
				(long)i+1, (long)ia[i], (long)i+1, a[i]);
		for(i = 0; i < m; i++)
			fprintf(fort,
				"\tdata clb(%ld)/%.16g/, cub(%ld)/%.16g/\n",
				(long)i+1, clb[i], (long)i+1, cub[i]);
		for(i = 0; i < n; i++)
			fprintf(fort,
				"\tdata lb(%ld)/%.16g/, ub(%ld)/%.16g/\n",
				(long)i+1, lb[i], (long)i+1, ub[i]);
		for(i = 0; i < n; i++)
			fprintf(fort, "\tdata c(%ld)/%.16g/\n", (long)i+1, c[i]);

		if (nint) {
		  fprintf(fort, "\tdata nint/%ld/, nsets/%ld/\n",
			(long)nint, (long)nsets);
		  for(i = 0; i < nsets; i++)
			fprintf(fort,
		"\tdata setind(%ld)/%ld/, prirty(%ld)/%ld/, type(%ld)/%ld/\n",
				(long)i+1, (long)setindex[i], (long)i+1,
				(long)priority[i],
				(long)i+1, (long)type[i]);
		  fprintf(fort, "\tdata setind(%ld)/%ld/, dlnpc/%ld*0.0001/\n",
			(long)nsets+1, (long)setindex[nsets], (long)nint);
		  for(i = 0; i < nint; i++)
			fprintf(fort,
				"\tdata intnms(%ld)/%ld/, dunpc(%ld)/%.16g/\n",
				(long)i+1, (long)intnums[i], (long)i+1, dunpc[i]);
		  }
		if (nelq) {
			fprintf(fort, "\tdata nelq/%ld/\n", (long)nelq);
			for(i = 0; i <= n; i++)
				fprintf(fort, "\tdata colq(%ld)/%ld/\n",
					(long)i+1, (long)colq[i]);
			for(i = 0; i < nelq; i++)
				fprintf(fort,
				  "\tdata rowq(%ld), delsq(%ld)/%ld, %.16g/\n",
					(long)i+1, (long)i+1, (long)rowq[i],
					delsq[i]);
			}
#if X0_IN_DSPACE
		if (x0 = *xpi0)
			for(i = 0; i < nxpi0; i++)
				fprintf(fort, "\tdata x0(%ld)/%.16g/\n",
					(long)i+1, x0[i]);
#endif
		fprintf(fort, "\n\tcall ekkdsca(rc, dspace, nwords, %ld)\n",
			(long)nmodels);
		fprintf(fort, "\tcall ekkdscm(rc, dspace, 1, %ld)\n",
			(long)nblock);
		for(i = 1; i <= Icontrolsize; i++)
			if (icset[i]) {
				fprintf(fort,
			"\tcall ekkiget(rc, dspace, ictl, %d)\n",
					Icontrolsize);
				for(; i <= Icontrolsize; i++)
					if (icset[i])
						fprintf(fort,
			"\tictl(%ld) = %ld\n", (long)i, (long)ic1[i-1]);
				fprintf(fort,
			"\tcall ekkiset(rc, dspace, ictl, %d)\n",
					Icontrolsize);
				break;
				}
		for(i = 1; i <= Rcontrolsize; i++)
			if (rcset[i]) {
				fprintf(fort,
			"\tcall ekkrget(rc, dspace, rctl, %d)\n",
					Rcontrolsize);
				for(; i <= Rcontrolsize; i++)
					if (rcset[i])
						fprintf(fort,
			"\trctl(%ld) = %.15g\n", (long)i, rc1[i-1]);
				fprintf(fort,
			"\tcall ekkrset(rc, dspace, rctl, %d)\n",
					Rcontrolsize);
				break;
				}
		fprintf(fort,
	"\tcall ekklmdl(rc, dspace, 2, m, n, nz, c, clb, cub, lb, ub,\n");
		fprintf(fort, "     1\t\t\tia, ja, a)\n");
		}
#endif
	sti->c = c;
	ekklmdl_(&rc, dspace, &I2, &m, &n, &nz, c, clb, cub, lb, ub,
		ia, ja, a);
	chk("ekklmdl");
	if (nelq) {
		msgcatch(3051L);
		msgcatch(7060L);
		ekkqmdl_(&rc, dspace, &I2, &nelq, rowq, colq, delsq);
		if (rc != 144 && rc != 114)
			chk("ekkqmdl");
#ifndef noDEBUG
		if (fort)
		  fprintf(fort,
		    "\tcall ekkqmdl(rc, dspace, 2, nelq, rowq, colq, delsq)\n");
#endif
		}
	get_statuses(dspace, sti);
#if X0_IN_DSPACE /* worked in OSL 1.2, may cause faults in OSL 2 */
	if (x0 = *xpi0) {
		ekknget_(&rc, dspace, L, &I7);
		chk("ekknget");
		if ((i = L[6] - 1) >= 0)
			memcpy((char *)(dspace + i), x0, nxpi0*sizeof(double));
#ifndef noDEBUG
		if (fort) {
		  fprintf(fort, "\tcall ekknget(rc, dspace, L, 7)\n");
		  fprintf(fort, "\tdo 1 i = 1, %d\n", n_var);
		  fprintf(fort, " 1\t\tdspace(L(7)+i-1) = x0(i)\n");
		  }
#endif
		}
#endif /* X0_IN_DSPACE */

	if (nint) {
#ifndef noDEBUG
		if (fort) {
			fprintf(fort,
	"\tcall ekkimdl(rc, dspace, nint, intnms, nsets, type, prirty,\n\
     1			nint, setind, intnms, dlnpc, dunpc)\n");
			}
#endif
		ekkimdl_(&rc, dspace, &nint, intnums, &nsets, type,
			priority, &ntotinfo, setindex, sets, dlnpc, dunpc);
		chk("ekkimdl");
		}

	if (outlev >= 3)
		msgzap1(6L);	/* page headings */

#ifndef noDEBUG
	if (fort)
		fflush(fort);
#endif
	if (mpsout) {
		ekkbcdo_(&rc, dspace, &I8, &I1, &mpsout);
		chk("ekkbcdo");
		ftncls_(&I8);
		}
	if (neednames) {
		/* assume f2c calling conventions for character variables */
		ekkname_(&rc, dspace, &I0, 0, &I1, &I0, 0, &I1, &I1);
		chk("ekkname");
		if (startbasis) {
			ekkbasi_(&rc, dspace, &startbasis);
			chk("ekkbasi");
			ftncls_(&startbasis);
			}
		}
	fflush(stdout);
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

 static int
iiscopy(fint i, fint j, int t, int *iis, int n, int *nbad, fint *table)
{
	int b, g, k;

	g = j - i + 1;
	if (!i || g <= 0)
		return 0;
	--iis;
	g = j - i + 1;
	for(b = 0; i <= j; i++)
		if ((k = (int)table[i]) < 0 || k > n)
			b++;
		else
			iis[k] = t;
	*nbad += b;
	return g - b;
	}

 static int
#ifdef KR_headers
send_iis(dspace, hbuf) double *dspace; char *hbuf;
#else
send_iis(double *dspace, char *hbuf)
#endif
{
	int *ci, i, j, nc, nv, *vi;
	fint index[9], rc, *table;
	static fint mask = 8, output = 2;

	if (!iisfind)
		return 0;
	index[0] = M + N;
	table = (fint*)M1alloc(index[0]*sizeof(fint));
	ekknfes_(&rc, dspace, &mask, &output, index, table);
	if (rc != 100)
		return Sprintf(hbuf, "\nReturn %ld from ekknfes.", (long)rc);
	ci = (int*)M1zapalloc((n_var+n_con)*sizeof(int));
	vi = ci + n_con;
	j = 0;
	--table;
	nc =	  iiscopy(index[1], index[2], 1, ci, n_con, &j, table)
		+ iiscopy(index[3], index[4], 3, ci, n_con, &j, table);
	nv =	  iiscopy(index[5], index[6], 1, vi, n_var, &j, table)
		+ iiscopy(index[7], index[8], 3, vi, n_var, &j, table);
	suf_iput("iis", ASL_Sufkind_con, ci);
	suf_iput("iis", ASL_Sufkind_var, vi);
	i = Sprintf(hbuf,
		"\nReturning %siis of %d variables and %d constraints.",
		j == 0 ? "" : "partial ", nv, nc);
	j += nc + nv;
	if (j != index[0])
		i += Sprintf(hbuf+i, "\n?? Got %d rather than %.d iis members.",
			j, (long)index[0]);
	return i;
	}

 static int
#ifdef KR_headers
send_ray(dspace, hbuf, wb) double *dspace; char *hbuf, *wb;
#else
send_ray(double *dspace, char *hbuf, char *wb)
#endif
{
	fint cv[30], rc;
	static fint I30 = 30;

	if (primaldual == 2)
		return 0;
	ekknget_(&rc, dspace, cv, &I30);
	chk("ekknget");
	suf_rput("unbdd", ASL_Sufkind_var, dspace + cv[29] - 1);
	return Sprintf(hbuf, "\n_var.unbdd %sreturned.", wb);
	}

 void
#ifdef KR_headers
amplout(mpre_trouble, sti) int mpre_trouble; Statinfo *sti;
#else
amplout(int mpre_trouble, Statinfo *sti)
#endif
{
	char buf[32], hbuf[480], *intfmt, *s, *simplex, *wb, *whatlim;
	fint L[7], bit, sit, sitbug;
	double *l, *le, *u, *x, *y, *y1, *ye, *z, *z1;
	double *dualobj, *primalobj, t;
	fint j, rc, status;
	int i, nint;
	typedef struct { char *msg; int code; int wantobj; } Sol_info;
	Sol_info *SI;
	static Sol_info solinfo[] = {
	 { "optimal solution",				000, 1 },
	 { "primal infeasible",				200, 1 },
	 { "primal unbounded",				300, 0 },
	 { "iteration limit",				400, 1 },
	 { "couldn't get feasible",			203, 1 },
	 { "solution limit",				101, 1 },
	 { "ran out of space",				500, 0 },
	 { "status unknown",				501, 1 },
	 { "bug!",					502, 0 },
	 { "iteration limit in mpre (osl's mixed-integer presolve)", 401, 0 },
	 { "best MIP solution so far restored",		101, 1 },
	 { "failed to restore best MIP solution",	503, 1 },
	 { "optimal (?) solution",			100, 1 },
	 { "relaxed LP is infeasible",			201, 0 },
	 { "relaxed LP is unbounded",			301, 0 },
	 { "dual infeasible",				202, 1 },
	 { "dual unbounded",				302, 0 }
		};
	static char *limname[] = { "iteration", "node", "gap" };
	static fint I0 = 0, I1 = 1;

	if (endbasis) {
		ekkbaso_(&rc, dspace, &endbasis, &I0);
		ftncls_(&endbasis);
		chk("ekkbaso");
		}
	irget();
	status = icontrol[46];
	sit = icontrol[3];
	sitbug = 0;
	if (dualused || method == 2) {
		simplex = "dual simplex";
		if (status == 2) {
			/* Cope with OSL BUG: dual simplex gives status 2 */
			/* for both infeasible and unbounded problems.	  */
			/* Run the primal simplex to tell them apart.	  */
			printf("Rerunning with primal simplex because of dual-simplex status bug.\n");
			ekksslv_(&rc, dspace, &I1, &I0);
			if (rc != 110 && rc != 100)
				chk("ekksslv");
			irget();
			status = icontrol[46];
			sitbug = icontrol[3] - sit;
			}
		}
	else
		simplex = "simplex";
	if (status > 6)
		status = 8;
	else if (status < 0)
		status = 7;
	else if (primaldual == 2 && (status == 1 || status == 2))
		status += 14;
	whatlim = 0;
	switch(mpre_trouble) {
		case 1: status = 9; break;
		case 2:
		case 3:
		case 4:	whatlim = limname[mpre_trouble-2];
			status = status == 0 ? 10 : 11;
		}

	if (!status) {
		if (intrestore_bug)
			status = 11;
		else if (nonconvex)
			status = 12;
		}
	else if (relaxed_infeas)
		status = relaxed_infeas + 12;
	SI = solinfo + status;
	solve_result_num = SI->code;
	i = Sprintf(hbuf, "%s: %s", Oinfo.bsname, SI->msg);
	if (whatlim)
		i += Sprintf(hbuf+i, "\n\tafter %s limit", whatlim);
	nint = nbv + niv + nlvoi;
	if (status == 3 && nint && icontrol[41] <= 0)
		i += Sprintf(hbuf+i, " without finding any feasible\n\
	integer solutions");
	x = y = 0;
	if (SI->wantobj) {
		primalobj = rcontrol + 17;
		if (barrier >= 0 && !nint) {
			if (primaldual == 2) {
				dualobj = primalobj;
				primalobj = rcontrol + 35;
				}
			else
				dualobj = rcontrol + 35 ;
			if (primalobj) {
				g_fmtop(buf, *primalobj);
				i += Sprintf(hbuf+i,
					"\nprimal objective %s", buf);
				}
			if (dualobj) {
				g_fmtop(buf, *dualobj);
				i += Sprintf(hbuf+i,
					"\n  dual objective %s", buf);
				}
			if (icontrol[59] < icontrol[3])
				i += Sprintf(hbuf+i,
					" (before final %s iterations)",
					simplex);
			}
		else {
			g_fmtop(buf, negobj ? -*primalobj : *primalobj);
			i += Sprintf(hbuf+i, status == 11
				? ";\nobjective %s before rounding to integer"
				: "; objective %s", buf);
			}

		j = 7;
		ekknget_(&rc, dspace, L, &j);
		chk("ekknget");
		if ((j = L[6]) > 0) {
			x = dspace + j - 1;
			if (nint) {
				/* round to integer */
				if (!relax) {
					y1 = x + n_var - nint;
					if (nshift) {
						round(y1, nbv - shift);
						round(y1 + nbv, niv);
						}
					else
						round(y1, nint);
					}
				}
			}
		if ((j = L[3]) > 0) {
			y = dspace + j - 1;
			if (t = objsen) {
				t = -1. / t;
				/* "nelq ||" patches OSL bug */
				if (nelq || primaldual == 2 && t > 0.)
					t = -t;
				if (negobj)
					t = -t;
				if (t != 1.) {
					ye = y + (primaldual == 2 ? N : M);
					for(y1 = y; y1 < ye; y1++)
						*y1 *= t;
					}
				}
			}
		}
	intfmt = "\n";
	if (barrier >= 0) {
		bit = icontrol[59];
		i += Sprintf(hbuf+i, "\n%ld barrier iterations", (long)bit);
		if (sit -= bit)
			i += Sprintf(hbuf+i, ", %ld %s iterations", (long)sit,
				simplex);
		}
	else if (nelq)
		i += Sprintf(hbuf+i, "\n%ld quadratic-programming iterations",
				(long)sit);
	else {
		i += Sprintf(hbuf+i,  "\n%ld %s%s iterations", (long)sit,
				netalg > 0 ? "network " : "", simplex);
		intfmt = "; ";
		}
	if (nint)
		i += Sprintf(hbuf+i, "%s%ld branch-and-bound nodes",
				intfmt, (long)icontrol[36]);
	if (intsols > 1)
		i += Sprintf(hbuf+i, "\n%ld feasible integer solutions found",
			(long)intsols);
	if (sitbug)
		i += Sprintf(hbuf+i,
		 "\n%ld primal simplex iterations to fix dual simplex status bug",
			(long)sitbug);
	if (primaldual == 2) {
		z = x;
		x = y;
		if (y = z) {
			z1 = z + n_con;
			l = LUrhs;
			u = Urhsx;
			/* fix dual variables */
			for(le = l + n_con; l < le; l++, u++, z++)
				if (*u > *l
				 && *l > osl_neginfinity
				 && *u < osl_infinity) {
					if (*z < -*z1)
						*z = *z1;
					z1++;
					}
			}
		if (rcontrol[2] < 0. && y)
			for(z1 = y + n_con; z1 > y; ) {
				--z1;
				*z1 = -*z1;
				}
		}
	else if (netsign) {
		if (y)
			for(j = 0; j < M; j++)
				if (netsign[j] == 2)
					y[j] = -y[j];
		}
	if (solve_result_num < 500) {
		wb = (Oinfo.wantsol || amplflag) ? "" : "would be ";
		if (s = send_statuses(dspace, sti))
			i += Sprintf(hbuf+i, "\n%s %sreturned.", s, wb);
		if (solve_result_num >= 200) {
			if (solve_result_num < 300)
				i += send_iis(dspace, hbuf+i);
			else if (solve_result_num < 400)
				i += send_ray(dspace, hbuf+i, wb);
			}
		}
	write_sol(hbuf, x, y, &Oinfo);
	}

 static double Times[4];

 static void
show_times(Void)
{
	int i;

	Times[3] = xectim_();
	for(i = 1; i <= 2; i++)
	    if (timing & i) {
		fprintf(i == 1 ? stdout : Stderr,
		"\nTimes (seconds):\nInput =  %g\nSolve =  %g\nOutput = %g\n",
			Times[1] - Times[0], Times[2] - Times[1],
			Times[3] - Times[2]);
		}
	}

 static double*
#ifdef KR_headers
kludge1(pnk, dspace) fint *pnk; double *dspace;
#else
kludge1(fint *pnk, double *dspace)
#endif
{
	double *sol0, *z, *z1;
	fint L[28], nk, rc, *x, *y;
	static fint I28 = 28;

	irget();
	ekknget_(&rc, dspace, L, &I28);
	chk("ekknget");
	nk = *pnk = icontrol[42];
	if (nk <= 0)
		return 0;
	z = z1 = (double *)Malloc(nk*(sizeof(double) + sizeof(fint)));
	y = (fint *)(z + nk);
	x = (fint *)dspace + (L[27] - 1);
	sol0 = dspace + (L[6] - 2);
	for(; nk > 0; --nk, x += 4)
		*z1++ = sol0[*y++ = *x];
	return z;
	}

 static void
#ifdef KR_headers
kludge2(nk, z, dspace) fint nk; double *z, *dspace;
#else
kludge2(fint nk, double *z, double *dspace)
#endif
{
	fint L[8], rc, *x;
	double *lb0, *ub0, *z0;
	static fint I8 = 8;

	irget();
	ekknget_(&rc, dspace, L, &I8);
	chk("ekknget");

	lb0 = dspace + (L[5] - 2);
	ub0 = dspace + (L[7] - 2);
	z0 = z;
	for(x = (fint *)(z + nk); --nk >= 0; x++)
		lb0[*x] = ub0[*x] = *z++;
	free(z0);
	}

 static void
#ifdef KR_headers
miplim_msg(aok) int aok;
#else
miplim_msg(int aok)
#endif
{
	switch(aok) {
	  case 2:
		printf("Simplex iteration limit (maxiter %ld)",
			(long)icontrol[4]);
		break;
	  case 3:
		printf("Node limit (maxnodes %ld)", (long)icontrol[39]);
		break;
	  case 4:
		printf("Mingap (%g) reached", mingap);
		}
	printf(icontrol[4] - icontrol[3] < iter_inc
		? "%s;\nallowing iter_inc = %ld\n\
more iterations to recover the best one.\n" : "%s.\n",
		" during branch-and-bound", (long)iter_inc);
	}

 static void
mb_opn(VOID)
{
	int k, pid;
	char *s, *td;

	k = (td = getenv("TMPDIR")) ? strlen(td) + 20 : 20;
	pid = getpid();
	munit.name = s = (char *)M1alloc(2*k);
	bunit.name = s + k;
	atexit(rmtemps);
	sigcatch_ASL();
	opn(&munit, td, pid);
	opn(&bunit, td, pid);
	}

 static void OSLCALLBACK
#ifdef KR_headers
qpslv(dspace, mspace, dobjval, jrtcod, jtype) double *dspace; fint *mspace; double *dobjval; fint *jrtcod; fint *jtype;
#else
qpslv(double *dspace, fint *mspace, double *dobjval, fint *jrtcod, fint *jtype)
#endif
{
	fint init = 1, rc = 0;
	Not_Used(mspace);

	if (*jtype == 2)
		init = 0;
	ekkqslv_(&rc, dspace, &qpnodealg, &init);
	irget();
	*jrtcod = icontrol[46];
	*dobjval = rcontrol[17];
	}

 static void OSLCALLBACK
#ifdef KR_headers
mipbslv(dspace, mspace, dobjval, jrtcod, jtype) double *dspace; fint *mspace; double *dobjval; fint *jrtcod; fint *jtype;
#else
mipbslv(double *dspace, fint *mspace, double *dobjval, fint *jrtcod, fint *jtype)
#endif
{
	static fint I3 = 3;
	fint init = 1, rc = 0;
	Not_Used(mspace);

	if (*jtype == 2)
		init = 0;
	ekkbslv_(&rc, dspace, &I3, &init);
	havebasis = 0;
	irget();
	*jrtcod = icontrol[46];
	*dobjval = rcontrol[17];
	}

 static int OSLCALLBACK
#ifdef KR_headers
heuu(dsp, mstatbin, mtobin, mfix, nbin, nfix1, nfix2, jtype)
	double *dsp;
	fint *mstatbin, *mtobin, *mfix, *nbin, *nfix1, *nfix2, *jtype;
#else
heuu(double *dsp, fint *mstatbin, fint *mtobin, fint *mfix, fint *nbin,
	fint *nfix1, fint *nfix2, fint *jtype)
#endif
{
	fint i, j, k, n;
	double t, *x0;

	if (*jtype != 1 || !binig1)
		return 0;
	binig1 = 0;
	n = N;
	x0 = X0 - 1;
	i = k = 0;
	while(i < n)
		if (j = mtobin[i++]) {
			if (!(t = x0[j]))
				mfix[k++] = -i;
			else if (t == 1.)
				mfix[k++] = i;
			}
	*nfix1 = k;
	*nfix2 = 0;
	return 0;
	}

 static char iis_table[] = "\n\
0	non	not in the iis\n\
1	low	at lower bound\n\
2	fix	fixed\n\
3	upp	at upper bound\n";

 static SufDecl
suftab[] = {
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

 int
#ifdef KR_headers
main(argc, argv) int argc; char **argv;
#else
main(int argc, char **argv)
#endif
{
	Statinfo sti;
	int aok, nint;
	fint balg, nk, rc;
	double *z;
	static fint I0 = 0, I1 = 1, I2 = 2, I7 = 7;
	static fint CBHEU = 5, CBSLV = 10, IBMPRS = 9;
	static fint netinitv[] = { 3, 1 };
	char *stub;
	Cbfunc *slvu;

	/* Now we get to work... */

	Times[0] = xectim_();

	ASL_alloc(ASL_read_fg);
	if ((stub = argv[1]) && !strcmp(stub, "-v"))
		WantPrintf1 = 1;
	getversion();
	if (!(stub = getstub(&argv, &Oinfo)))
		usage_ASL(&Oinfo,1);

	suf_declare(suftab, sizeof(suftab)/sizeof(SufDecl));
	amplin(stub, argv, &sti);

	Times[1] = xectim_();
	msgcatch1(101L);
	msgcatch1(7034L);

	if (nint = nbv + niv + nlvoi) {
		if (binconv)
			if (nelq)
				binconv = 0;
			else {
				ekkbmpr_(&rc, dspace, &IBMPRS, &binconv);
				chk("ekkbmpr");
				if (fort)
					fprintf(fort,
				"\n\tcall ekkbmpr(rc, dspace, %ld, %ld)\n",
						(long)IBMPRS, (long)binconv);
				if (nircset) {
					/* ekkbmpr may have screwed up */
					/* some nondefault values, so  */
					/* set them again to be safe.  */
					irget();
					ircset();
					irset();
					}
				}
		}
	else
		binconv = 0;
	if (simplify >= 0) {
		ekkprsl_(&rc, dspace, &I7, &simplify);
		chk("ekkprsl");
#ifndef noDEBUG
		if (fort)
			fprintf(fort, "\tcall ekkprsl(rc, dspace, 7, %ld)\n",
				(long)simplify);
#endif
		}
	if (netalg > 0) {
		if (simplify >= 0) {
			ekkpssl_(&rc, dspace, &I7);
			chk("ekkpssl");
			}
		ekknslv_(&rc, dspace, &netalg, netinitv + netinit);
#ifndef noDEBUG
		if (fort)
			fprintf(fort, "\tcall ekknslv(rc, dspace, %ld, %ld)\n",
				(long)netalg, (long)netinitv[netinit]);
#endif
		goto solved;
		}
	if (scale) {
		ekkscal_(&rc, dspace);
		chk("ekkscal");
#ifndef noDEBUG
		if (fort)
			fprintf(fort, "\tcall ekkscal(rc, dspace)\n");
#endif
		}

	if (crash > 0) {
		if (crash < 5) {
			ekkcrsh_(&rc, dspace, &crash);
			chk("ekkcrsh");
#ifndef noDEBUG
			if (fort)
				fprintf(fort,
					"\tcall ekkcrsh(rc, dspace, %ld)\n",
					(long)crash);
#endif
			}
		else if (nint > 0 || nelq > 0) {
			balg = crash == 5 ? 0 : -1;
			ekkbslv_(&rc, dspace, &balg, &I2);
			/* rc == 100 means infeasible problem */
			/* rc == 110 means there was a warning */
			/* about ill-conditioning */
			switch(rc) {
				case 100:
					relaxed_infeas = 1;
					goto unsimp;
				case 110:
					break;
				default:
					chk("ekkbslv");
				}
			}
		}

	if (method == 2 && pricing > 1)
		pricing = 1;
	slvu = 0;
	if (sti.nsets) {
		if (nelq) {
			qpnodealg = barrier >= 0 ? -1 : 1;
			slvu = qpslv;
			/*goto nosimp;*/
			}
		else if (barrier >= 0) {
			slvu = mipbslv;
			ekknwmt_(&rc, dspace, &I2);
			chk("ekknwmt");
			if (fort)
				fprintf(fort, "\tcall ekknwmt(rc, dspace, 2)\n");
			}
		msgcatch1(105L);	/* good solution found */
		if (mipdisplay && outlev < 3) {
			msgcatch(87L);	/* "good" node */
			msgcatch(102L);	/* infeasible node */
			}
		if (slvu)
			ekkrgcb(&CBSLV, slvu);
		if (binig1 = binig)
			ekkrgcb(&CBHEU, (Cbfunc*)heuu);
		if (binsimp && nbv) {
			nk = 3 - binsimp;
			ekkmpre_(&rc, dspace, &nk);
			if (rc != 113)
				chk("ekkmpre");
#ifndef noDEBUG
			if (fort)
				fprintf(fort,
				"\tcall ekkmpre(rc, dspace, %ld)\n", (long)nk);
#endif
			/* check for exceeded iteration limit */
			irget();
			if (icontrol[3] >= icontrol[4]) {
				amplout(1, &sti);
				goto done;
				}
			}
 /*nosimp:*/
		if (mbfile)
			mb_opn();
		need_rm = 1;
		msgcatch1(3011L); /* catch OSL memory-overwrite bug */
		msgcatch1(3L);	  /* catch OSL infeasibility message */
		ekkmslv_(&rc, dspace, &I1, munit.unit + mbfile,
			bunit.unit + mbfile);
		if (binig)
			ekkclcb(&CBHEU);
		if (slvu)
			ekkclcb(&CBSLV);

/*** ekkmslv gives rc = 123 (but no error message and a correct solution)
 *** on the following:

	set I := 1..4; var b{I} binary;
	minimize zot: sum{i in I} i*b[i];
	suffix sosno; suffix ref;
	let{i in I} b[i].sosno := 1;
	let{i in I} b[i].ref := i;
	s.t. zap: sum{i in I} b[i] >= 1;

*** so we ignore return 123 from ekkmslv... */

		if (rc != 100 && rc != 113 && rc != 123)
			chk("ekkmslv");
		msgzap1(3L);
#ifndef noDEBUG
		if (fort)
			fprintf(fort,
			"\tcall ekkmslv(rc, dspace, 1, %ld, %ld)\n",
				(long)munit.unit[mbfile],
				(long)bunit.unit[mbfile]);
#endif
		irget();
		if (icontrol[46] == 3 && icontrol[41] > 0 && iter_inc > 0) {
			/* Try to recover the best integer solution */
			aok = mingapped ? 4 : icontrol[3] >= icontrol[4] ? 2 : 3;
			if (outlev >= 1)
				miplim_msg(aok);
			if (icontrol[4] - icontrol[3] < iter_inc)
				icontrol[4] += iter_inc;
			irset();
			ekksslv_(&rc, dspace, &I1, &I2);
			if (rc != 110 && rc != 100)
				chk("ekksslv");
			}
		rmtemps();
		}
	else if (nelq) {
		msgcatch(3063L);
		msgcatch(7068L);
		msgzap1(3035L);
		balg = barrier >= 0 ? -1 : 1;
		ekkqslv_(&rc, dspace, &balg, startbasis ? &I0 : &I1);
		switch(rc) {
			case 100:
			case 110:
			case 114:
			case 132:
				break;
			default:
				chk("ekkqslv");
			}
#ifndef noDEBUG
		if (fort)
			fprintf(fort, "\tcall ekkqslv(rc, dspace, %ld, %d)\n",
				(long)balg, startbasis ? 0 : 1);
#endif
		}
	else if (barrier >= 0) {
		if (bdisplay > 0) {
			mipfreq = bdisplay;
			msgcatch(153L);
			msgcatch(204L);
			}
		if (nint)
			bs_switch = 2;
		if ((balg = barrier) == 4)
			balg = -1;
		ekkbslv_(&rc, dspace, &balg, &bs_switch);
		havebasis = bs_switch >= 2 && bs_switch <= 4;
		/* rc == 100 means infeasible problem */
		/* rc == 110 means there was a warning about ill-conditioning */
		switch(rc) {
			case 100:
				if (nint) {
					irget();
					relaxed_infeas = (int)icontrol[46];
					if (relaxed_infeas < 0
					 || relaxed_infeas > 2)
						relaxed_infeas = 0;
					}
			case 110:
				break;
			default:
				chk("ekkbslv");
			}
#ifndef noDEBUG
		if (fort)
			fprintf(fort, "\tcall ekkbslv(rc, dspace, %ld, %ld)\n",
				(long)balg, (long)bs_switch);
#endif
		}
	else {
		ekksslv_(&rc, dspace, &method, &pricing);
		switch(rc) {
			case 100:
			case 110:
				break;
			default:
				chk("ekksslv");
			}
#ifndef noDEBUG
		if (fort)
			fprintf(fort, "\tcall ekksslv(rc, dspace, %ld, %ld)\n",
				(long)method, (long)pricing);
#endif
		}

	aok = 0;

 unsimp:
	if (simplify >= 0) {
		if (nint)
			z = kludge1(&nk, dspace);
		ekkpssl_(&rc, dspace, &I7);
		chk("ekkpssl");
#ifndef noDEBUG
		if (fort)
			fprintf(fort, "\tcall ekkpssl(rc, dspace, 7)\n");
#endif
		if ((barrier < 0 || bs_switch) && !relaxed_infeas) {
			if (nint)
				kludge2(nk, z, dspace);
			ekksslv_(&rc, dspace, &method, &pricing);
			if (rc != 110 && rc != 100)
				chk("ekksslv");
#ifndef noDEBUG
			if (fort)
				fprintf(fort,
			"\tcall ekksslv(rc, dspace, %ld, %ld)\n",
					(long)method, (long)pricing);
#endif
			}
		}
	if (binconv) {
		if (fort) {
			fprintf(fort, "\n\tcall ekkbmps(rc, dspace, %ld)\n",
				(long)IBMPRS);
			fflush(fort);
			/* flush now in case ekkbmps faults... */
			}
		ekkbmps_(&rc, dspace, &IBMPRS);
		chk("ekkbmps");
		if (fort)
			fprintf(fort, "\n\tcall ekkbmps(rc, dspace, %ld)\n",
				(long)IBMPRS);
		}
 solved:
#ifndef noDEBUG
	if (fort) {
		fprintf(fort, "\n\twrite(*,*) 'rc = ', rc\n\
	call ekkiget(rc, dspace, ictl, %d)\n\
	call ekkrget(rc, dspace, rctl, %d)\n\
	write(*,*) 'Iprobstat = ', ictl(47)\n\
	write(*,*) 'Robjvalue = ', rctl(18)\n", Icontrolsize, Rcontrolsize);
		fprintf(fort, "\tcall ekknget(rc, dspace, L, 7)\n\
	if (L(7) .gt. 0) then\n\
		xoff = L(7) - 1\n\
		write(*,*) 'solution:'\n\
		f = 0\n\
		do 10 i = 1, n\n\
			write(*,*) i, dspace(i+xoff)\n\
			f = f + c(i)*dspace(i+xoff)\n\
 10			continue\n\
		write(*,*) 'c**T * x  = ', f\n\
		endif\n");
		fprintf(fort, "\tend\n");
		}
#endif
	Times[2] = xectim_();

	amplout(aok, &sti);
 done:
	show_times();

	return 0;
	}

#undef fmt
#undef fprintf
#undef printf
#undef sprintf

 int
fprintf(FILE *f, const char *fmt, ...)
{
	int rc;
	va_list ap;
	if (f != stdout)
		goto accept;
	if (rc = WantPrintf) {
		Nchatter++;	/* Who knows when OSL will chatter? */
 accept:
		va_start(ap, fmt);
		rc = Vfprintf(f, fmt, ap);
		va_end(ap);
		}
	return rc;
	}

 int
printf(const char *fmt, ...)
{
	int rc;
	va_list ap;
	va_start(ap, fmt);
	rc = Vfprintf(stdout, fmt, ap);
	va_end(ap);
	Nchatter++;
	return rc;
	}

 int
sprintf(char *s, const char *fmt, ...)
{
	int rc;
	va_list ap;
	va_start(ap, fmt);
	rc = Vsprintf(s, fmt, ap);
	va_end(ap);
	return rc;
	}
