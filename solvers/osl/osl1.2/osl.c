/****************************************************************
Copyright (C) 1992-1997 Lucent Technologies
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

/* Ideally I would #include "ekkc.h",
   but I nearly barfed when I looked at it.  */

/* #define noDEBUG to omit the "fortran" keyword ("fortran=foo"
		makes osl write Fortran to file foo.  The Fortran
		currently omits the kludge needed for recovering
		the solution with "simplify" >= 0 when there are
		integer variables.) */

/* #define NETBUGfixed once EKKNSLV always gives dual values
		of signs consistent with the duals from EKSSLV.
		#defining NETBUGfixed removes the "netbug" keyword. */

/* #define MPREBUGfixed once EKKMPRE can be used with Special Ordered
		sets of Type 2 */

/* #define MAXSOLS_BUG_FIXED once Imaxsols doesn't have to be incremented
		when it is > 1 */

#ifdef NDP
#ifndef COPYRIGHT_NOTICE
#define COPYRIGHT_NOTICE \
	"\r\nPortions of this program (c) 1987, 1991 Microway, Inc.\r\n"
#endif
char Copyright_Notice[] = COPYRIGHT_NOTICE;
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

#ifdef Use_tolower
#include "ctype.h"
#define Tolower(x) tolower(x)
#define lc_init()
#else
static unsigned char lc[256];
 static void
lc_init(Void)
{
	register int i;
	register char *s;
	for(i = 0; i < 256; i++)
		lc[i] = i;
	for(s = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"; *s; s++)
		lc[*s] = *s + 'a' - 'A';
	}
#define Tolower(x) lc[x]
#endif

#ifdef SVS
#include "osl.h"
extern void _cinitf77(int,char**);
#else
#define _cinitf77(a,b) Not_Used(a)
#endif


#ifdef NDP
/* Stuff for NDP C */
int xargc;
char **xargv;
int zero = 0, one = 1;
extern void fmt_init(Void), rec_init(Void), rec_uninit(Void);
void (*initrec)(Void) = rec_init, (*uninitrec)(Void) = rec_uninit;
char *tempfiles[100] = {0};

/* End of NDP stuff (for now) */
#else
#define rec_uninit() /* nothing */
#include "osl.h"
#endif

extern void ekkbasi_ ANSI((fint*, double*, fint*));
extern void ekkbaso_ ANSI((fint*, double*, fint*, fint*));
extern void ekkbcdo_ ANSI((fint*, double*, fint*, fint*, fint*));
extern void ekkbslv_ ANSI((fint*, double*, fint*, fint*));
extern void ekkcrsh_ ANSI((fint*, double*, fint*));
extern void ekkdsca_ ANSI((fint*, double*, fint*, fint*));
extern void ekkdscm_ ANSI((fint*, double*, fint*, fint*));
extern void ekkiget_ ANSI((fint*, double*, fint*, fint*));
extern void ekkimdl_ ANSI((fint*, double*, fint*, fint*, fint*, fint*, fint*,
			fint*, fint*, fint*, double *, double *));
extern void ekkiset_ ANSI((fint*, double*, fint*, fint*));
extern void ekklmdl_ ANSI((fint *rc, double *dspace, fint *type, fint *nr,
			fint *nc, fint *nel, double *c, double *clb,
			double *cub, double *lb, double *ub,
			fint *ia, fint *ja, double *a));
extern void ekkmpre_ ANSI((fint*, double*, fint*));
extern void ekkmset_ ANSI((fint *rc, double *dspace, fint *strtnum,
			fint *maxalw, fint *maxprt, fint *trace,
			fint *usrexit, fint *endnum, fint *nonum));
extern void ekkmslv_ ANSI((fint*, double*, fint*, fint*, fint*));
extern void ekkname_ ANSI((fint*, double*, fint*, char*, fint*, fint*, char*,
			fint*, fint*, fint, fint));
extern void ekknslv_ ANSI((fint*, double*, fint*, fint*));
extern void ekknget_ ANSI((fint*, double*, fint*, fint*));
extern void ekkprsl_ ANSI((fint*, double*, fint*, fint*));
extern void ekkpssl_ ANSI((fint*, double*, fint*));
extern void ekkqmdl_ ANSI((fint*, double*, fint*, fint*, fint*, fint*, double*));
extern void ekkqslv_ ANSI((fint*, double*, fint*, fint*));
extern void ekkrget_ ANSI((fint*, double*, double*, fint*));
extern void ekkrset_ ANSI((fint*, double*, double*, fint*));
extern void ekkscal_ ANSI((fint*, double*));
extern void ekksslv_ ANSI((fint*, double*, fint*, fint*));

extern char *getenv ANSI((const char *));
extern void ftncls_ ANSI((fint*));
extern void ftnopn_ ANSI((fint*,char*,fint*,fint*,fint));

#ifdef NDP
#define ftncls_(x) /* */
#endif

fint convex_qp = 1, outlev = 1;
static double osl_infinity = 1e31, osl_neginfinity = -1e31;
static char *netsign;
static int nshift, shift;
static int intrestore_bug, relaxed_infeas;

#ifdef bunit_and_munit_in_TMPDIR
static char *tmpdir, *tmpstub;
static int pid;
#ifndef KR_headers
#include "unistd.h"	/* for getpid() */
#endif
#endif

#ifndef NETBUGfixed
static fint netbug;
#endif

 static fint M, N, NO, NZ;
#ifndef noDEBUG
 static FILE *fort;
#endif
 static fint	bfile = 1, iter_inc = 2000, mfile = 1, pricing = 2,
		prestrat = 1, sos2 = 1;
 static fint	branch, bs_switch, dspinc, method, mpsout, nelq, netinit,
		scale, timing, trace, wantns2;
 static fint bdisplay, endbasis, mipdisplay, simplexinit, startbasis;
 static fint mipfreq = 20;
 static double *dspace;
 static double mingap, rowinc;
 static fint dspqfac = 10.;
 static int need_minmax = 1;
 static int mingapped, need_rm, negobj, negobj1, relax, solving;
 static fint barrier = -1, binsimp = 1, crash = -1, dual_thresh = 32000,
		netalg = 1, simplify = -1;
 static fint	bunit[] = { 0, 3 },
		munit[] = { 0, 4 };
 static int intsols, nonconvex, primaldual;
#define Icontrolsize 61
#define Rcontrolsize 45
 static double	rcontrol[Rcontrolsize], rc1[Rcontrolsize];
 static fint	icontrol[Icontrolsize], ic1[Icontrolsize];
 static char	icset[Icontrolsize+1], rcset[Rcontrolsize+1];
 extern char *progname;
 static fint passx0 = 1;

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
	Barrier		= { -1, 3, &barrier },
	Bdisplay	= { 0, LONG_MAX, &bdisplay },
	Bfile		= { 0, 1, &bfile },
	Binsimp		= { 0, 2, &binsimp },
	Branch		= { 0, 15, &branch },
	Bs_switch	= { 0, 4, &bs_switch },
	Convex_qp	= { 0, 1, &convex_qp },
	Dspinc		= { 0, LONG_MAX, &dspinc },
	Dspqfac		= { 0, 1000, &dspqfac },
	Dualthresh	= { LONG_MIN, LONG_MAX, &dual_thresh },
#ifdef Use_basunits
	Endbasis	= { 9, 99, &endbasis },
#endif
	Iadjactype	= { 0, 1, &ic1[13], 14 },
	Icrashtype	= { 0, 4, &crash},
	Idensecol	= { 10, LONG_MAX, &ic1[15], 16 },
	Idevexmode	= { -3, 3, &ic1[16], 17 },
	Idroprowct	= { 1, 30, &ic1[18], 19 },
	Ifastits	= { LONG_MIN, LONG_MAX, &ic1[38], 39 },
	Iformntype	= { 0, 1, &ic1[14], 15 },
	Iheurpass	= { 0, LONG_MAX, &ic1[54], 55 },
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
	Mfile		= { 0, 1, &mfile },
#ifdef msgDebug
	Mipdisplay	= { 0, 3, &mipdisplay },
#else
	Mipdisplay	= { 0, 2, &mipdisplay },
#endif
	Mipfreq		= { 1, LONG_MAX, &mipfreq },
	Mpsout		= { 0, 2, &mpsout },
	Netalg		= { -1, 2, &netalg },
#ifndef NETBUGfixed
	Netbug		= { 0, 1, &netbug },
#endif
	Netinit		= { 0, 1, &netinit },
	Objno		= { 0, LONG_MAX, 0 },
	Outlev		= { 0, 4, &outlev },
	Passx0		= { 0, 1, &passx0 },
	Prestrat	= { 0, 3, &prestrat },
	Pricing		= { 1, 3, &pricing },
	Scale		= { 0, 1, &scale },
	Simplexinit	= { 0, 2, &simplexinit },
	Simplify	= { -1, 3, &simplify},
	Sos2		= { 0, 1, &sos2 },
#ifdef Use_basunits
	Startbasis	= { 9, 99, &startbasis },
#endif
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
rm(n) fint n;
#else
rm(fint n)
#endif
{
#ifdef bunit_and_munit_in_TMPDIR
	char *b, *bs, buf[16];

	if (b = tmpdir)
		bs = tmpstub;
	else
		b = bs = buf;
	Sprintf(bs, "fort%d.%ld", pid, n);
	remove(b);
#else
	char buf[16];
#ifdef WATCOM
	Sprintf(buf, "for%03ld", n);
#else
#ifdef HPUX
#define FORT_fmt n < 10 ? "ftn0%ld" : "ftn%ld"
#else
#define FORT_fmt "fort.%ld"
#endif
	Sprintf(buf, FORT_fmt, n);
#endif
	ftncls_(&n);
	remove(buf);
#endif
	}

 static void
rmtemps(Void)
{
	if (need_rm) {
		need_rm = 0;
		if (mfile)
			rm(munit[1]);
		if (bfile)
			rm(bunit[1]);
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
	rec_uninit();	/* NDP shutdown */
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

#ifdef bunit_and_munit_in_TMPDIR
 static void
#ifdef KR_headers
opn(n) fint n;
#else
opn(fint n)
#endif
{
	fint ierr, len;
	static fint I1 = 1;

	Sprintf(tmpstub, "fort%d.%ld", pid, n);
	len = strlen(tmpdir);
	ftnopn_(&n, tmpdir, &I1, &ierr, len);
	if (ierr)
		cantopen(tmpdir, "scratch");
	}
#endif

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
		fprintf(Stderr, /*(*/ "[%ld, %ld]).\n", it->tmin, it->tmax);
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

#ifdef Use_basunits
#define sf_bas sf_int
#else

 struct Basinfo {
	fint create;
	fint unit;
	fint *up;
	};
 typedef struct Basinfo Basinfo;
 static Basinfo
	Startbasis = { 0, 9, &startbasis },
	Endbasis = { 1, 10, &endbasis };
#if 0
	subroutine ftnopn(n, name, create, ierr)
	integer n, ierr
	character*(*) name
	logical create

	open(n,file=name,status='OLD',iostat=ierr)
	if (ierr .ne. 0 .and. create) then
		open(n,file=name,status='NEW',iostat=ierr)
		endif
	if (ierr .eq. 0) rewind n
	end
#endif

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
	ftnopn_(&bi->unit, buf, &bi->create, &ierr, len);
	if (ierr)
		cantopen(buf, kw->name);
	else
		*bi->up = bi->unit;
	return v;
	}
#endif

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
 static keyword keywds[] = {
	{ "adjactype",	sf_int,		VP &Iadjactype },
	{ "barrier",	sf_int,		VP &Barrier },
	{ "bb_bfile",	sf_int,		VP &Bfile },
	{ "bb_mfile",	sf_int,		VP &Mfile },
	{ "bbcutoff",	sf_dbl,		VP &Rbcutoff },
	{ "bbdispfreq",	sf_int,		VP &Mipfreq },
	{ "bbdisplay",	sf_int,		VP &Mipdisplay },
	{ "bbimprove",	sf_dbl,		VP &Rimprove },
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
	{ "dualthresh",	sf_int,		VP &Dualthresh },
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
	{ "ilinelen",	sf_int,		VP &Ilinelen },
	{ "imajorits",	sf_int,		VP &Imajorits },
	{ "imaxfactor",	sf_int,		VP &Imaxfactor },
	{ "imaxiter",	sf_int,		VP &Imaxiter },
	{ "imaxiterb",	sf_int,		VP &Imaxiterb },
	{ "imaxnodes",	sf_int,		VP &Imaxnodes },
	{ "imaxprojns",	sf_int,		VP &Imaxprojns },
	{ "imaxsols",	sf_int,		VP &Imaxsols },
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
	{ "mpsout",	sf_int,		VP &Mpsout },
	{ "mslv_bfile",	sf_int,		VP &Bfile },
	{ "mslv_mfile",	sf_int,		VP &Mfile },
	{ "mufactor",	sf_dbl,		VP &Rmufactor },
	{ "muinit",	sf_dbl,		VP &Rmuinit },
	{ "mulimit",	sf_dbl,		VP &Rmulimit },
	{ "mulinfac",	sf_dbl,		VP &Rmulinfac },
	{ "netalg",	sf_int,		VP &Netalg },
#ifndef NETBUGfixed
	{ "netbug",	sf_int,		VP &Netbug },
#endif
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
	{ "simplex",	sf_int,		VP &Method },
	{ "simplexinit",sf_int,		VP &Simplexinit },
	{ "simplify",	sf_int,		VP &Simplify },
	{ "sos2",	sf_int,		VP &Sos2 },
	{ "startbasis",	sf_bas,		VP &Startbasis },
	{ "stderr",	sf_file,	VP &Stderr },
	{ "stepmult",	sf_dbl,		VP &Rstepmult },
	{ "target",	sf_dbl,		VP &Rtarget },
	{ "timing",	sf_int, 	VP &set_timing },
	{ "toldinf",	sf_dbl,		VP &Rtoldinf },
	{ "tolint",	sf_dbl,		VP &Rtolint },
	{ "tolpinf",	sf_dbl,		VP &Rtolpinf },
	{ "trace",	sf_int,		VP &Trace },
	{ "version",	Ver_val,	0 },
	{ "wantsol",	WS_val,		0, "write .sol file (without -AMPL)" }
	};

 extern char osl_vers[];

 static Option_Info Oinfo = {"osl", "OSL 1.2", "osl_options",
				keywds, nkeywds, 0, osl_vers+1};

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
		fprintf(Stderr, "\nReturn code %ld from %s\n", rc, who);
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

 static fint
#ifdef KR_headers
nsets2(ia,ja,np,m,setindex,intnums,nint1,a,clb,cub,dunpc,nsets1,lb,ub,c,objadj)
 fint *ia, *ja, *np, *m, *setindex, *intnums, *nint1, nsets1;
 double *a, *clb, *cub, *dunpc, *lb, *ub, *c;
 int objadj;
#else
nsets2(fint *ia, fint *ja, fint *np, fint *m, fint *setindex, fint *intnums,
	fint *nint1, double *a, double *clb, double *cub, double *dunpc,
	fint nsets1, double *lb, double *ub, double *c, int objadj)
#endif
{
	fint i, i1, j, jlim, k, k0, k1, L, L1, m0, mdec, n, n0, nint, rv;
	fint *ja0, *ja00, *ja1, *z;
	double *c1, *lb1, *ub1;
	static fint inc;

	--ia;
	--a;
	n = *np;
	ja00 = ja;
	ja0 = ja += n0 = ampl_options[3] - 1;
	nint = niv + nbv;
	jlim = n - n0 - nint - objadj;
	m0 = *m;
	k0 = ampl_options[4];
	if (setindex) {
		if (inc) {
			dunpc += inc;
			if (nsets1)
				setindex += nsets1 - 1;
			else
				nsets1 = 1;
			}
		else
			--setindex;
		intnums += inc;
		z = (fint *)Malloc((m0+objadj)*sizeof(fint));
		L = L1 = mdec = 0;
		while(L < k0)
			z[L++] = ++L1;
		if (!inc || jlim == nbv)
			msgzap1(3068L); /* spurious "no 0-1 variables" */
		}
	for(rv = 1;;) {
		k = ia[i = *ja];
		if (k <= k0)
			k = ia[++i];
		if (setindex) {
			j = *++ja;
			if (ia[j] <= k0)
				j++;
			ja1 = ja;
			k1 = k + 1;
			*++setindex = ja - ja0 + inc;
			*dunpc++ = ia[++i] == k1 ? -a[i] : 0.;
			*dunpc++ = ia[++j] == k1 ? -a[j] : 0.;
			for(;;) {
				if ((i1 = ia[i = *++ja]) <= k0)
					i1 = ia[++i];
				if (i1 != k)
					break;
				*dunpc++ = ia[++i] == k1 ? -a[i] : 0.;
				}
			k1 = k - 1;
			if (L < k1) {
				mdec += k1 - L;
				do z[L++] = 0;
					while(L < k1);
				}
			z[L++] = ++L1;
			z[L++] = ++L1;
			if (ja - ja0 >= jlim) {
				k = (ja - ja1) << 1;
				mdec += ++k;
				while(--k >= 0)
					z[L++] = 0;
				break;
				}
			}
		else
			for(;;) {
				if ((i1 = ia[i = *++ja]) <= k0)
					i1 = ia[++i];
				if (i1 != k)
					break;
				}
		rv++;
		if (ja - ja0 >= jlim)
			break;
		}
	if (setindex) {
		if (objadj)
			z[L++] = L1 + 1;
		*m -= mdec;
		ja1 = ja0 + nbv;
		L1 = 1;
		--z;
		for(ja = ja00; ja < ja1; ja++) {
			i = *ja;
			*ja = L1;
			for(j = ja[1]; i < j; i++)
				if (k = z[ia[i]]) {
					ia[L1] = k;
					a[L1++] = a[i];
					}
			}
		if (k1 = nshift = niv + objadj) {
			shift = jlim;
			n = n0 + nbv;
			ja1 += jlim;
			i = *ja1;
			lb += n;
			ub += n;
			c += n;
			lb1 = lb + jlim;
			ub1 = ub + jlim;
			c1 = c + jlim;
			while(--k1 >= 0) {
				*ja++ = L1;
				*lb++ = *lb1++;
				*ub++ = *ub1++;
				*c++ = *c1++;
				for(j = *++ja1; i < j; i++)
					if (k = z[ia[i]]) {
						ia[L1] = k;
						a[L1++] = a[i];
						}
				}
			}
		*ja = L1;
		--clb;
		--cub;
		++L;
		for(i = 1; i < L; i++)
			if (j = z[i]) {
				clb[j] = clb[i];
				cub[j] = cub[i];
				}
		free((char *)(z+1));
		*np -= jlim;
		n = n0 + jlim;
		while(n0 < n)
			*intnums++ = ++n0;
		}
	else {
		inc = nbv + niv - jlim;
		if (nint1)
			*nint1 = inc;
		if (rv > 1 && !inc)
			--rv;
		}
	return rv;
	}

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

 void
#ifdef KR_headers
amplin(stub, av) char *stub, **av;
#else
amplin(char *stub, char **av)
#endif
{
	fint L[7], MXROW, MXCOL;
	fint *colq, *ia, *ia0, *ia1, *ia2, *intnums, *ja, *ja0, *ja1,
		*priority, *rowq, *setindex, *type;
	int *num, *pri, *start;
	double *a, *a0, *a1, *a2, *c, *c1, *clb, *cub, *cub1,
		*delsq, *dlnpc, *dunpc, *lb, *lb1, *lu, *ub, *ub1, *x0, **xpi0;
	double os, t;
	ograd *og;
	fint Lextra, i, intoff, j, k, k1, k2, kk, m, m0,
		n, n0, nblock, nextra, nint, nint1, nnz,
		nsets, nsets0, nsets1, nsets1a,
		nwords, nz, nzextra, rc;
	FILE *nl;
	char *ns1, *z;
	int objadj;
	static fint I0 = 0, I1 = 1, I2 = 2, I243 = 243, I256 = 256,
			I48 = 48, I50 = 50, I7 = 7, I8 = 8;
	static double intinf = 2147483647.;

	nl = jacdim(stub, &M, &N, &NO, &NZ, &MXROW, &MXCOL, (fint)strlen(stub));

	if (N <= 0) {
		fprintf(Stderr, "%s has no variables\n", filename);
		fexit(4);
		}
	*stub_end = 0;
	objadj = 0;

	/* Allow space for adding a row to adjust the objective value. */
	m = M + 1;
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

	LUrhs = clb = (double *)Malloc((2*m + n)*sizeof(double));
	Urhsx = cub = LUrhs + m;
	c = cub + m;
	LUv = lb = (double *)Malloc((2*n + nz)*sizeof(double)
					+ (nz + n + 1)*sizeof(fint));
	Uvx = ub = LUv + n;
	A_vals = a = ub + n;
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
	if (need_nl) {
		need_nl = 0;	/* Who knows when OSL will chatter? */
		printf("\n");
		}
	nelq = qpcheck(&rowq, &colq, &delsq);

	if (nint = nint1 = nbv + niv) {
		if (ampl_options[3] && ampl_options[3] <= N - nint)
			wantns2 = sos2;
		if (!relax || wantns2)
			primaldual = 1;
		}
	if (!primaldual)
		primaldual = m - n > dual_thresh ? 2 : 1;

	if (i = nint + nlvbi + nlvci + nlvoi)
		if (relax)
			relax = i;
		else {
			if (nelq) {
				fprintf(Stderr,
				  "Sorry, osl can't handle integer QPs\n");
				fexit(1);
				}
			netalg = 0;
			}
	else
		relax = 0;
	if (nelq) {
		primaldual = 1;
		netalg = 0;
		barrier = -1;
		}

	if (primaldual == 1) {
		xpi0 = &X0;
		Free(&pi0);
		}
	else {
		xpi0 = &pi0;
		Free(&X0);
		}

	if (!startbasis && !crash)
		switch(simplexinit) {
			case 0: passx0 = 0; break;
			case 3: pricing = 3; /* no break */
			case 1: passx0 = 1;
			}
	if (!passx0)
		Free(xpi0);

	/* check for pure network */

	if (barrier >= 0)
		netalg = 0;
	else if (netalg)
		netalg = primaldual == 2 ? 0 : netcheck(a, ia, ja);

	nwords = 11000 + 10*(M + N + NZ + 3) + dspinc + nelq*dspqfac;
	/* kludge: OSL sometimes memory faults without nelq*dspqfac. */
	Lextra = 0;
	m = M;
	n = N;
	if (nint) {
		nwords += 40000 + 10*rowinc*(M + NZ);
		nsets1 = nsets1a = mip_pri(&start, &num, &pri, 1000L);
		nsets = nsets0 = wantns2
			? nsets2(ia,ja,&n,&m,0,0,&nint1,0,0,0,0,0,0,0,0,0) : 1L;
		if (nsets1) {
			j = n - nint + nint1;
			k = nint1;
			for(i = 0; i < nsets1; )
					k -= num[i++];
			nsets += nsets1;
			if (k > 0 && !relax)
				nsets1a++;
			else
				nsets--;
			}
		Lextra = nint*(2*sizeof(double) + sizeof(fint))
				+ (3*nsets+1)*sizeof(fint);
		}
	if (dspinc)
		printf("Allocating %ld double words for dspace.\n", nwords);
	dspace = (double *)Malloc(nwords*sizeof(double) + Lextra);
	memset((char *)dspace, 0, nwords*sizeof(double) + Lextra);

	msgzap(48L, 50L);
	msgzap1(81L);
	if (outlev & 1) {
		ekkmset_(&rc,dspace,&I1,&I0,&I0,&I0,&I0,&I243,&I1);
		chk("ekkmset");
		}
	if (outlev < 3) {
		msgzap1(3001L);	/* msg about unbounded problem.	*/
		msgzap1(3048L);	/*  "	 "	 "	  "	*/
		if (outlev <= 0)
			msgzap(1L, 243L);
		else {
			msgzap1(  2L);
			msgzap1(  6L);
			msgzap(  16L,   20L);
			msgzap1( 23L);
			msgzap1( 31);
			msgzap(  33L,	35L);
			msgzap1( 37L);
			msgzap(  75L,   78L);
			msgzap1( 82L);
			msgzap(  87L,  171L);
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

	ekkdsca_(&rc, dspace, &nwords, &I1);
	chk("ekkdsca");
	msgcatch1(7003L);
	nblock = nelq ? 5 : 1;
	ekkdscm_(&rc, dspace, &I1, &nblock);
	chk("ekkdscm");

	if (nint && rowinc > 0) {
		icset[9] = 1;
		ic1[8] = M + rowinc*M;
		}

	irget();
	icontrol[24] = 1;	/* msg no's at right */
	/* song and dance because we can only call ekkdsca once */
	if (prestrat != 1 || branch) {
		ic1[51] = prestrat | branch >> 2;
		icset[52] = 1;
		}
	for(i = 0; i < Icontrolsize; i++)
		if (icset[i+1])
			icontrol[i] = ic1[i];
	for(i = 0; i < Rcontrolsize; i++)
		if (rcset[i+1])
			rcontrol[i] = rc1[i];

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
			objadj = 1;
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
#ifndef MAXSOLS_BUG_FIXED
	if (i = icontrol[40])
		icontrol[40] = 2*i - 1;
#endif
	irset();

	nonlin(nlc, "nonlinear constraints");
	nonlin(plterms, "piecewise-linear terms");

	if (netalg) {
		ns1 = netsign - 1;
		if (i = 2*n - nz) {
			clb[m] = osl_neginfinity;
			/* osl gets confused by osl_infinity if we've added */
			/* a constraint to adjust the objective value. */
			netsign[m] = 2;
			cub[m++] = obj_adj ? intinf : osl_infinity;
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
			/* kludge around OSL bug: upper bound of +infinity */
			/* seems to be treated as 1 */
		double *ubi, *ube;
		ube = ub + N;
		for(ubi = ub + (N - nint); ubi < ube; ubi++)
			if (*ubi >= osl_infinity)
				*ubi = intinf;

		dlnpc = dspace + nwords;
		dunpc = dlnpc + nint;
		intnums = (fint *)(dunpc + nint);
		type = intnums + nint;
		priority = type + nsets;
		setindex = priority + nsets;
		intoff = N - nint + 1;
		for(i = 0; i < nint; i++) {
			dlnpc[i] = dunpc[i] = 0.0001;
			intnums[i] = i + intoff;
			}
		for(i = 0; i < nsets; i++)
			priority[i] = 1000;
		*type = 4;
		setindex[0] = 1;
		setindex[nsets] = nint + 1;
		if (wantns2) {
			Free(xpi0);	/*!! for now; should adjust, use */
			n0 = n;
			nsets2(ia, ja, &n, &m, setindex, intnums, &nint1, a,
				clb, cub, dunpc, nsets1a, lb, ub, c, objadj);
			nz = ja[n] - 1;
			if ((i = nint1 > 0) < nsets) {
#ifndef MPREBUGfixed
				binsimp = 0;	/* avoid OSL bug */
#endif
				do type[i] = 2;
					while(++i < nsets);
				}
			}
		if (relax) {
			nint1 = nint;
			if (wantns2) {
				if (nint1 -= nint = n0 - n) {
					--nsets;
					++type;
					}
				dlnpc += nint1;
				dunpc += nint1;
				intnums += nint1;
				setindex[nsets] = nint + 1;
				}
			else
				nint = nbv = niv = 0;
			if (nint1)
				printf(
				"ignoring integrality of %ld variable%s\n",
					nint1, nint1 > 1 ? "s" : "");
			}
		else if (nsets1) {
			z = (char *)Malloc(nint);
			memset(z, 0, nint);
			k2 = 0;
			--intoff;
			j = N - (shift + nshift);
			for(i = 0; i < nsets1; i++) {
				type[i] = 4;
				setindex[i] = k2 + 1;
				k1 = start[i];
				k = k1 - intoff;
				if (k1 > j)
					k1 -= shift;
				priority[i] = 1001 - pri[i];
				kk = num[i];
				while(--kk >= 0) {
					z[k++] = 1;
					intnums[k2++] = ++k1;
					}
				}
			if (nsets1a > nsets1) {
				type[i] = 4;
				setindex[i++] = k2 + 1;
				++intoff;
				for(k = 0; k < nint; k++)
					if (!z[k])
						intnums[k2++] = k + intoff;
				}
			setindex[i] = k2 + 1;
			free(z);
			}
		}
	if (primaldual == 2) {
		if ((x0 = *xpi0) && (t = -objsen)) {
			if (t > 0.)
				t = -t;
			if (negobj)
				t = -t;
			if (t != 1.)
				for(i = 0; i < m; i++)
					x0[i] *= t;
			}
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
		free((char *)LUv);
		}
#ifndef noDEBUG
	if (fort) {
		fprintf(fort, "\tprogram bugsho\n\n");
		fprintf(fort,
			"\tinteger i, L(7), xoff\n\tdouble precision f\n\n");
		fprintf(fort, "\tdouble precision dspace(%ld), rctl(%d)\n",
			nwords, Rcontrolsize);
		fprintf(fort, "\tinteger ia(%ld), ictl(%d), ja(%ld)\n",
			nz, Icontrolsize, n+1);
		fprintf(fort, "\tinteger m, n, nwords, nz, rc\n");
		fprintf(fort,
		  "\tdouble precision a(%ld), c(%ld), clb(%ld), cub(%ld)\n",
				nz, n, m, m);
		fprintf(fort, "\tdouble precision lb(%ld), ub(%ld)\n", n, n);
		if (nint) {
			fprintf(fort, "\tinteger nint, nsets\n");
			fprintf(fort,
	"\tinteger intnms(%ld), prirty(%ld), setind(%ld), type(%ld)\n",
				nint, nsets, nsets+1, nsets);
			fprintf(fort,
				"\tdouble precision dlnpc(%ld), dunpc(%ld)\n",
				nint, nint);
			}
		if (*xpi0)
			fprintf(fort, "\tdouble precision x0(%ld)\n", n);
		if (nelq)
			fprintf(fort, "\tinteger nelq, rowq(%ld), colq(%ld)\n\
	double precision delsq(%ld)\n",
				nelq, n+1, nelq);
		fprintf(fort,"\tdata m/%ld/, n/%ld/, nwords/%ld/, nz/%ld/\n",
			m, n, nwords, nz);
		for(i = 0; i <= n; i++)
			fprintf(fort, "\tdata ja(%ld)/%ld/\n", i+1, ja[i]);
		for(i = 0; i < nz; i++)
			fprintf(fort, "\tdata ia(%ld)/%ld/, a(%ld)/%.16g/\n",
				i+1, ia[i], i+1, a[i]);
		for(i = 0; i < m; i++)
			fprintf(fort,
				"\tdata clb(%ld)/%.16g/, cub(%ld)/%.16g/\n",
				i+1, clb[i], i+1, cub[i]);
		for(i = 0; i < n; i++)
			fprintf(fort,
				"\tdata lb(%ld)/%.16g/, ub(%ld)/%.16g/\n",
				i+1, lb[i], i+1, ub[i]);
		for(i = 0; i < n; i++)
			fprintf(fort, "\tdata c(%ld)/%.16g/\n", i+1, c[i]);

		if (nint) {
		  fprintf(fort, "\tdata nint/%ld/, nsets/%ld/\n",
			nint, nsets);
		  for(i = 0; i < nsets; i++)
			fprintf(fort,
		"\tdata setind(%ld)/%ld/, prirty(%ld)/%ld/, type(%ld)/%ld/\n",
				i+1, setindex[i], i+1, priority[i],
				i+1, type[i]);
		  fprintf(fort, "\tdata setind(%ld)/%ld/, dlnpc/%ld*0.0001/\n",
			nsets+1, setindex[nsets], nint);
		  for(i = 0; i < nint; i++)
			fprintf(fort,
				"\tdata intnms(%ld)/%ld/, dunpc(%ld)/%.16g/\n",
				i+1, intnums[i], i+1, dunpc[i]);
		  }
		if (nelq) {
			fprintf(fort, "\tdata nelq/%ld/\n", nelq);
			for(i = 0; i <= n; i++)
				fprintf(fort, "\tdata colq(%ld)/%ld/\n",
					i+1, colq[i]);
			for(i = 0; i < nelq; i++)
				fprintf(fort,
				  "\tdata rowq(%ld), delsq(%ld)/%ld, %.16g/\n",
					i+1, i+1, rowq[i], delsq[i]);
			}
		if (x0 = *xpi0)
			for(i = 0; i < n; i++)
				fprintf(fort, "\tdata x0(%ld)/%.16g/\n",
					i+1, x0[i]);

		fprintf(fort, "\n\tcall ekkdsca(rc, dspace, nwords, 1)\n");
		fprintf(fort, "\tcall ekkdscm(rc, dspace, 1, %ld)\n", nblock);
		for(i = 1; i <= Icontrolsize; i++)
			if (icset[i]) {
				fprintf(fort,
			"\tcall ekkiget(rc, dspace, ictl, %d)\n",
					Icontrolsize);
				for(i; i <= Icontrolsize; i++)
					if (icset[i])
						fprintf(fort,
			"\tictl(%ld) = %ld\n", i, ic1[i-1]);
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
				for(i; i <= Rcontrolsize; i++)
					if (rcset[i])
						fprintf(fort,
			"\trctl(%ld) = %.15g\n", i, rc1[i-1]);
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
	if (x0 = *xpi0) {
		ekknget_(&rc, dspace, L, &I7);
		chk("ekknget");
		if ((i = L[6] - 1) >= 0)
			memcpy((char *)(dspace + i), x0, n*sizeof(double));
#ifndef noDEBUG
		if (fort) {
		  fprintf(fort, "\tcall ekknget(rc, dspace, L, 7)\n");
		  fprintf(fort, "\tdo 1 i = 1, %d\n", n_var);
		  fprintf(fort, " 1\t\tdspace(L(7)+i-1) = x0(i)\n");
		  }
#endif
		}

	if (nint) {
#ifndef noDEBUG
		if (fort) {
			fprintf(fort,
	"\tcall ekkimdl(rc, dspace, nint, intnms, nsets, type, prirty,\n\
     1			nint, setind, intnms, dlnpc, dunpc)\n");
			}
#endif
		ekkimdl_(&rc, dspace, &nint, intnums, &nsets, type,
			priority, &nint, setindex, intnums, dlnpc, dunpc);
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
#ifndef NDP
		ftncls_(&I8);
/* with Fortran
	subroutine ftncls(n)
	close(n)
	end
*/
#endif
		}
	if (startbasis) {
		/* assume f2c calling conventions for character variables */
		ekkname_(&rc, dspace, &I0, 0, &I1, &I0, 0, &I1, &I1, 8L, 8L);
		chk("ekkname");
		ekkbasi_(&rc, dspace, &startbasis);
		chk("ekkbasi");
		ftncls_(&startbasis);
		}
	fflush(stdout);
	}

 static double *
#ifdef KR_headers
copyup(x) double *x;
#else
copyup(double *x)
#endif
{
	double *x0, *x1, *x2;

	x0 = x1 = (double *)Malloc(n_var*sizeof(double));
	x2 = x1 + (n_var - shift - nshift);
	while(x1 < x2)
		*x1++ = *x++;
	x2 += shift;
	while(x1 < x2)
		*x1++ = 0;
	x2 += nshift;
	while(x1 < x2)
		*x1++ = *x++;
	return x0;
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

 void
#ifdef KR_headers
amplout(mpre_trouble) int mpre_trouble;
#else
amplout(int mpre_trouble)
#endif
{
	char buf[32], hbuf[320], *intfmt,  *whatlim;
	fint L[7], bit, sit;
	double *l, *le, *u, *x, *y, *y1, *ye, *z, *z1;
	double *dualobj, *primalobj, t;
	fint j, rc, status;
	int i, nint;
	static char *statmsg[] = {
		"optimal solution",
		"primal infeasible",
		"primal unbounded",
		"iteration limit",
		"couldn't get feasible",
		"solution limit",
		"ran out of space",
		"status unknown",
		"bug!",
		"iteration limit in mpre (osl's mixed-integer presolve)",
		"best MIP solution so far restored",
		"failed to restore best MIP solution",
		"optimal (?) solution",
		"relaxed LP is infeasible or unbounded"
		};
	static char *limname[] = { "iteration", "node", "gap" };
	static char wantobj[] = {1,1,0,1,1,1,0,1,0,0,1,1,1,0};
	static fint I0 = 0;


	if (endbasis) {
		ekkbaso_(&rc, dspace, &endbasis, &I0);
		ftncls_(&endbasis);
		chk("ekkbaso");
		}
	irget();
	status = icontrol[46];
	if (status > 6)
		status = 8;
	else if (status < 0)
		status = 7;
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
		else if (barrier >= 0
			&& barrier != 1 && icontrol[59] >= icontrol[12])
			status = 3;	/* circumvent OSL bug */
		else if (nonconvex)
			status = 12;
		}
	else if (relaxed_infeas)
		status = 13;
	i = Sprintf(hbuf, "OSL 1.2: %s", statmsg[status]);
	if (whatlim)
		i += Sprintf(hbuf+i, "\n\tafter %s limit", whatlim);
	nint = nbv + niv;
	if (status == 3 && nint && icontrol[41] <= 0)
		i += Sprintf(hbuf+i, " without finding any feasible\n\
	integer solutions");
	x = y = 0;
	if (wantobj[status]) {
		primalobj = rcontrol + 17;
		if (barrier >= 0 && !nint) {
			if (primaldual == 2) {
				dualobj = primalobj;
				primalobj = rcontrol + 35;
				}
			else
				dualobj = rcontrol + 35 ;
			if (barrier == 1)
				if (primaldual == 2)
					primalobj = 0;
				else
					dualobj = 0;
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
					" (before final simplex iterations)");
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
				if (wantns2)
					x = copyup(x);
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
	sit = icontrol[3];
	if (barrier >= 0) {
		bit = icontrol[59];
		i += Sprintf(hbuf+i, "\n%ld barrier iterations", bit);
		if (sit -= bit)
			i += Sprintf(hbuf+i, ", %ld simplex iterations", sit);
		intfmt = "\n";
		}
	else if (nelq)
		i += Sprintf(hbuf+i, "\n%ld quadratic-programming iterations",
				sit);
	else {
		i += Sprintf(hbuf+i,  "\n%ld %ssimplex iterations", sit,
				netalg > 0 ? "network " : "");
		intfmt = "; ";
		}
	if (nint)
		i += Sprintf(hbuf+i, "%s%ld branch-and-bound nodes",
				intfmt, icontrol[36]);
	if (intsols > 1)
		Sprintf(hbuf+i, "\n%ld feasible integer solutions found",
			intsols);
	if (primaldual == 2) {
		z = x;
		x = y;
		y = z;
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
		if (rcontrol[2] < 0.)
			for(z1 = y + n_con; z1 > y; ) {
				--z1;
				*z1 = -*z1;
				}
		}
	else if (netsign) {
#ifdef NETBUGfixed
#define j2 2
#else
		/* Kludge around bug in OSL's computation of */
		/* dual values for network problems: negate the rows */
		int j2 = 2 - netbug;
#endif
		if (y)
			for(j = 0; j < M; j++)
				if (netsign[j] == j2)
					y[j] = -y[j];
#undef j2
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
	double *lb0, *ub0;
	static fint I8 = 8;

	irget();
	ekknget_(&rc, dspace, L, &I8);
	chk("ekknget");

	lb0 = dspace + (L[5] - 2);
	ub0 = dspace + (L[7] - 2);
	for(x = (fint *)(z + nk); --nk >= 0; x++)
		lb0[*x] = ub0[*x] = *z++;
	}

 void
#ifdef KR_headers
ekkmsgu_(dsp, Idsp, msgnum, nrvec, rvec, nivec, ivec, ncvec, cvec, cvec_len)
	double *dsp, *rvec;
	fint *Idsp, *msgnum, *nrvec, *nivec, *ivec, *ncvec, cvec_len;
	char *cvec;
#else
ekkmsgu_(double *dsp, fint *Idsp, fint *msgnum, fint *nrvec,
	 double *rvec, fint *nivec, fint *ivec, fint *ncvec,
	 char *cvec, fint cvec_len)
#endif
{
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
#ifdef msgDebug
    if (mipdisplay > 2) {
	fint i, n, n_cv, n_iv, n_rv;

	printf("ekkmsgu(%ld,%ld,%ld,%ld,%ld)\n", mnum, n_rv = *nrvec,
		n_iv = *nivec, n_cv = *ncvec, cvec_len);
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
					"Warning: ", not_convex, ivec[1]);
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
#ifdef KR_headers
miplim_msg(aok) int aok;
#else
miplim_msg(int aok)
#endif
{
	switch(aok) {
	  case 2:
		printf("Simplex iteration limit (maxiter %ld)", icontrol[4]);
		break;
	  case 3:
		printf("Node limit (maxnodes %ld)", icontrol[39]);
		break;
	  case 4:
		printf("Mingap (%g) reached", mingap);
		}
	printf(icontrol[4] - icontrol[3] < iter_inc
		? "%s;\nallowing iter_inc = %ld\n\
more iterations to recover the best one.\n" : "%s.\n",
		" during branch-and-bound", iter_inc);
	}

 int
#ifdef KR_headers
main(argc, argv) int argc; char **argv;
#else
main(int argc, char **argv)
#endif
{
	int aok, nint;
	fint nk, rc;
	double *z;
	static fint I0 = 0, I1 = 1, I2 = 2, I7 = 7;
	static fint netinitv[] = { 3, 1 };
	char *stub;

#ifdef NDP
	/* More NDP stuff */
	static char *xxargv[2];
	xargc = 1;	/* don't let Fortran see args */
	xxargv[0] = argv[0];
	xargv = xxargv;
	fmt_init();
	rec_init();
#else
	_cinitf77(argc, argv);
#endif

	/* Now we get to work... */

	Times[0] = xectim_();

	ASL_alloc(ASL_read_fg);
	if (!(stub = getstub(&argv, &Oinfo)))
		usage_ASL(&Oinfo,1);
	lc_init();
	amplin(stub, argv);

	Times[1] = xectim_();
	msgcatch1(101L);
	msgcatch1(7034L);
	if (simplify >= 0) {
		ekkprsl_(&rc, dspace, &I7, &simplify);
		chk("ekkprsl");
#ifndef noDEBUG
		if (fort)
			fprintf(fort, "\tcall ekkprsl(rc, dspace, 7, %ld)\n",
				simplify);
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
				netalg, netinitv[netinit]);
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
	if (crash > 0 && barrier < 0) {
		ekkcrsh_(&rc, dspace, &crash);
		chk("ekkcrsh");
#ifndef noDEBUG
		if (fort)
			fprintf(fort, "\tcall ekkcrsh(rc, dspace, %ld)\n",
				crash);
#endif
		}

	nint = nbv + niv;

	if (method == 2 && pricing > 1)
		pricing = 1;
	if (barrier >= 0) {
		if (bdisplay > 0) {
			mipfreq = bdisplay;
			msgcatch(153L);
			msgcatch(204L);
			}
		if (nint)
			bs_switch = 2;
		ekkbslv_(&rc, dspace, &barrier, &bs_switch);
		/* rc == 100 means infeasible problem */
		/* rc == 110 means there was a warning about ill-conditioning */
		switch(rc) {
			case 100:
				relaxed_infeas = 1;
			case 110:
				break;
			default:
				chk("ekkbslv");
			}
#ifndef noDEBUG
		if (fort)
			fprintf(fort, "\tcall ekkbslv(rc, dspace, %ld, %ld)\n",
				barrier, bs_switch);
#endif
		}
	else if (nelq) {
		msgcatch(3063L);
		msgcatch(7068L);
		msgzap1(3035L);
		ekkqslv_(&rc, dspace, &I1, startbasis ? &I0 : &I1);
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
			fprintf(fort, "\tcall ekkqslv(rc, dspace, 1, %d)\n",
				startbasis ? 0 : 1);
#endif
		}
	else {
		ekksslv_(&rc, dspace, &method, &pricing);
		switch(rc) {
			case 100:
				relaxed_infeas = 1;
			case 110:
				break;
			default:
				chk("ekksslv");
			}
#ifndef noDEBUG
		if (fort)
			fprintf(fort, "\tcall ekksslv(rc, dspace, %ld, %ld)\n",
				method, pricing);
#endif
		}

	aok = 0;
	if (nint && !relaxed_infeas) {
		msgcatch1(105L);	/* good solution found */
		if (mipdisplay && outlev < 3) {
			msgcatch(87L);	/* "good" node */
			msgcatch(102L);	/* infeasible node */
			}
		if (binsimp && nbv) {
			nk = 3 - binsimp;
			ekkmpre_(&rc, dspace, &nk);
			if (rc != 113)
				chk("ekkmpre");
#ifndef noDEBUG
			if (fort)
				fprintf(fort,
				"\tcall ekkmpre(rc, dspace, %ld)\n", nk);
#endif
			/* check for exceeded iteration limit */
			irget();
			if (icontrol[3] >= icontrol[4]) {
				amplout(1);
				goto done;
				}
			}
#ifdef bunit_and_munit_in_TMPDIR
		{
		char *td;
		fint k;
		/* Need to open bunit and munit as direct files		*/
		/* ==> need to know record length (RECL, which the	*/
		/* OSL manual does not provide).  This code is		*/
		/* commented out, pending info from IBM about RECL.	*/
		if (mfile | bfile && (td = getenv("TMPDIR")) && *td) {
			pid = getpid();
			k = strlen(td);
			tmpdir = (char *)Malloc(k+16);
			memcpy(tmpdir, td, k);
			tmpstub = tmpdir + k;
			if (tmpstub[-1] != '/')
				*tmpstub++ = '/';
			if (mfile)
				opn(munit[1]);
			if (bfile)
				opn(bunit[1]);
			}
		}
#endif
		need_rm = 1;
		msgcatch1(3011L); /* catch OSL memory-overwrite bug */
		msgcatch1(3L);	  /* catch OSL infeasibility message */
		ekkmslv_(&rc, dspace, &I1, munit + mfile, bunit + bfile);
		if (rc != 100 && rc != 113)
			chk("ekkmslv");
		msgzap1(3L);
#ifndef noDEBUG
		if (fort)
			fprintf(fort,
			"\tcall ekkmslv(rc, dspace, 1, %ld, %ld)\n",
				munit[mfile], bunit[bfile]);
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
					method, pricing);
#endif
			}
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

	amplout(aok);
 done:
	show_times();

	rec_uninit();	/* NDP shutdown */
	return 0;
	}
