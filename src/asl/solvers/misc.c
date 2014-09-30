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

/* include vararg/stdarg stuff first to avoid trouble with C++ */
#include "stddef.h"
#include "stdarg.h"
#include "asl.h"

#ifdef __cplusplus
extern "C" {
#endif

 extern void mpec_adjust_ASL(ASL*);
 extern void obj_adj_ASL(ASL*);


real edagread_one = 1.;
char *progname;
ASL *cur_ASL;
ASLhead ASLhead_ASL = {&ASLhead_ASL, &ASLhead_ASL};

 static char anyedag[] = "fg_read (or one of its variants)";
 static char psedag[] = "pfg_read, pfgh_read, or jacpdim";

 ASL *
set_cur_ASL(ASL *a)
{
	ASL *rv = cur_ASL;
	cur_ASL = a;
	return rv;
	}

  ASL *
get_cur_ASL(void)
{ return cur_ASL; }

 void
exit_ASL(EdRead *R, int n)
{
	Jmp_buf *J;
	if ((J = R->asl->i.err_jmp_) && n > 0)
		longjmp(J->jb, n);
	exit(n);
	}

 void
scream(EdRead *R, int n, const char *fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
	vfprintf(Stderr, fmt, ap);
	va_end(ap);
	exit_ASL(R, n);
	}

 static real
notread(const char *what, const char *pred)
{
	fprintf(Stderr, "\n*** %s called before %s.\n", what, pred);
	exit(1);
	/* not reached */
	return 0;
	}

 static real
obj0val(ASL *a, int nobj, real *X, fint *nerror)
{
	Not_Used(a);
	Not_Used(nobj);
	Not_Used(X);
	Not_Used(nerror);
	return notread("objval", anyedag);
	}

 static void
obj0grd(ASL *a, int nobj, real *X, real *G, fint *nerror)
{
	Not_Used(a);
	Not_Used(nobj);
	Not_Used(X);
	Not_Used(G);
	Not_Used(nerror);
	notread("objgrd", anyedag);
	}

 static void
con0val(ASL *a, real *X, real *R, fint *nerror)
{
	Not_Used(a);
	Not_Used(X);
	Not_Used(R);
	Not_Used(nerror);
	notread("conval", anyedag);
	}

 static void
jac0val(ASL *a, real *X, real *J, fint *nerror)
{
	Not_Used(a);
	Not_Used(X);
	Not_Used(J);
	Not_Used(nerror);
	notread("jacval", anyedag);
	}

 static real
con0ival(ASL *a, int i, real *X, fint *nerror)
{
	Not_Used(a);
	Not_Used(i);
	Not_Used(X);
	Not_Used(nerror);
	notread("conival", anyedag);
	return 0.;
	}

 static real
conivalmap(ASL *a, int i, real *X, fint *nerror)
{
	int *cm;
	if ((cm = a->i.cmap))
		i = cm[i];
	return a->p.Conival_nomap(a, i, X, nerror);
	}

 static void
congrdmap(ASL *a, int i, real *X, real *G, fint *nerror)
{
	int *cm;
	if ((cm = a->i.cmap))
		i = cm[i];
	a->p.Congrd_nomap(a, i, X, G, nerror);
	}

 static int
lcon0val(ASL *a, int i, real *X, fint *nerror)
{
	Not_Used(a);
	Not_Used(i);
	Not_Used(X);
	Not_Used(nerror);
	notread("conival", anyedag);
	return 0;
	}

 static void
con0grd(ASL *a, int i, real *X, real *G, fint *nerror)
{
	Not_Used(a);
	Not_Used(i);
	Not_Used(X);
	Not_Used(G);
	Not_Used(nerror);
	notread("congrd", anyedag);
	}

 static void
hv0comp(ASL *a, real *hv, real *p, int nobj, real *ow, real *y)
{
	Not_Used(a);
	Not_Used(hv);
	Not_Used(p);
	Not_Used(nobj);
	Not_Used(ow);
	Not_Used(y);
	notread("hvcomp", "pfgh_read or fgh_read");
	}

 static void
hv0compd(ASL *a, real *hv, real *p, int co)
{
	Not_Used(a);
	Not_Used(hv);
	Not_Used(p);
	Not_Used(co);
	notread("hvcompd", "pfgh_read or fgh_read");
	}

 static varno_t
hv0comps(ASL *a, real *hv, real *p, int co, varno_t nz, varno_t *z)
{
	Not_Used(a);
	Not_Used(hv);
	Not_Used(p);
	Not_Used(co);
	Not_Used(nz);
	Not_Used(z);
	notread("hvcomps", "pfgh_read or fgh_read");
	return 0;
	}

 static void
hv0init(ASL *a, int n, int no, real *ow, real *y)
{
	Not_Used(a);
	Not_Used(n);
	Not_Used(no);
	Not_Used(ow);
	Not_Used(y);
	notread("hvinit", "pfgh_read");
	}

 static void
hes0set(ASL *a, int flags, int obj, int nobj, int con, int ncon)
{
	Not_Used(a);
	Not_Used(flags);
	Not_Used(obj);
	Not_Used(nobj);
	Not_Used(con);
	Not_Used(ncon);
	notread("duthes, fullhes, or sputhes", "pfgh_read or jacpdim");
	}

 static void
x0known(ASL *a, real *x, fint *nerror)
{
	Not_Used(a);
	Not_Used(x);
	Not_Used(nerror);
	notread("xknown", psedag);
	}

 static void
dut0hes(ASL *a, real *H, int nobj, real *ow, real *y)
{
	Not_Used(a);
	Not_Used(H);
	Not_Used(nobj);
	Not_Used(ow);
	Not_Used(y);
	notread("duthes", "pfgh_read or jacpdim");
	}

 static void
ful0hes(ASL *a, real *H, fint LH, int nobj, real *ow, real *y)
{
	Not_Used(a);
	Not_Used(H);
	Not_Used(LH);
	Not_Used(nobj);
	Not_Used(ow);
	Not_Used(y);
	notread("fullhes", "pfgh_read or jacpdim");
	}

 static void
sut0hes(ASL *a, SputInfo **p, real *H, int nobj, real *ow, real *y)
{
	Not_Used(a);
	Not_Used(p);
	Not_Used(H);
	Not_Used(nobj);
	Not_Used(ow);
	Not_Used(y);
	notread("sputhes", "pfgh_read or jacpdim");
	}

 static fint
sut0set(ASL *a, SputInfo **p, int nobj, int have_ow, int have_y, int both)
{
	Not_Used(a);
	Not_Used(p);
	Not_Used(nobj);
	Not_Used(have_ow);
	Not_Used(have_y);
	Not_Used(both);
	notread("sputset", "pfgh_read or jacpdim");
	return 0;
	}

Edagpars edagpars_ASL = {
	{0,0},	/* h */
	1.,	/* hffactor */
	5,	/* FUNNEL_MIN */
	5,	/* maxfwd */
	1,	/* need_funcadd */
	100,	/* vrefGulp */
#ifdef ASL_OLD_DERIV_ERR_CHECK
	2,	/* want_derivs */
#else
	1,	/* want_derivs */
#endif
	12,	/* ihd_limit */
	-1,	/* solve_code */
	obj0val,
	obj0val,
	obj0grd,
	obj0grd,
	con0val,
	jac0val,
	con0ival,
	con0ival,
	con0grd,
	con0grd,
	hv0comp,
	hv0comp,
	hv0compd,
	hv0comps,
	hv0init,
	hv0init,
	hes0set,
	lcon0val,
	x0known,
	dut0hes,
	dut0hes,
	ful0hes,
	ful0hes,
	sut0hes,
	sut0hes,
	sut0set
	};

 int
edag_peek(EdRead *R)
{
	int c;
	R->Line++;
	R->lineinc = 0;
	R->rl_buf[0] = c = getc(R->nl);
	return c;
	}

 static void
eatcr(FILE *nl)
{
	int c;

	while((c = getc(nl)) == '\r');
	if (c >= 0 && c != '\n')
		ungetc(c, nl);
	}

 char *
read_line(EdRead *R)
{
	char *s, *se;
	int x;
	char *rv;
	FILE *nl = R->nl;

	s = R->rl_buf;
	se = s + sizeof(R->rl_buf) - 1;
	if (R->lineinc)
		R->Line++;
	else {
		s++;
		R->lineinc = 1;
		}
	rv = s;
	for(;;) {
		x = getc(nl);
		if (x < ' ') {
			if (x < 0) {
 eof:
				if (R->can_end)
					return 0;
				fprintf(Stderr,
				 "Premature end of file, line %ld of %s\n",
					R->Line, R->asl->i.filename_);
				exit_ASL(R,1);
				}
			if (x == '\n')
				break;
			if (x == '\r') {
				eatcr(nl);
				break;
				}
			}
		*s++ = x;
		if (s >= se) {
			for(;;) {
				x = getc(nl);
				if (x == '\r') {
					eatcr(nl);
					goto eol;
					}
				if (x == '\n')
					goto eol;
				if (x < 0)
					goto eof;
				}
			}
		}
 eol:
	*s = 0;
	return rv;
	}

 static void
memfailure(const char *who, const char *what, size_t len)
{
	if (progname)
		fprintf(Stderr, "%s: ", progname);
	fprintf(Stderr, "%s(%lu) failure: %s.\n", who, (long)len, what);
	exit(1);
	}

static char	ran_out[] =	"ran out of memory";

 void *
mymalloc(size_t len)
{
	void *rv;
	static char who[] = "malloc";
	rv = malloc(len);
	if (!rv) {
		/* Defend against stupid systems: malloc(0) */
		/* should return a nonzero value.  Routines in */
		/* amplsolver.a should never call malloc(0), but */
		/* solvers may do so. */
		if (!len)
			rv = malloc(sizeof(real));
		if (!rv)
			memfailure(who, ran_out, len);
		}
	return rv;
	}

 void *
myralloc(void *rv, size_t len)
{
	static char who[] = "realloc";
	rv = realloc(rv, len);
	if (!rv) {
		if (!len)
			rv = malloc(sizeof(real));
		if (!rv)
			memfailure(who, ran_out, len);
		}
	return rv;
	}

 void
what_prog(void)
{
	if (progname)
		fprintf(Stderr, "%s: ", progname);
	}

 void
badread(EdRead *R)
{
	what_prog();
	fprintf(Stderr, "error reading line %ld of %s:\n\t", R->Line, R->asl->i.filename_);
	}

 void
badline(EdRead *R)
{
	ASL *asl = R->asl;
	FILE *nl;
	char *s, *se;
	int x;

	fprintf(Stderr, "bad line %ld of %s", R->Line, filename);
	if (xscanf == ascanf) {
		if (!R->lineinc) {
			nl = R->nl;
			s = R->rl_buf;
			se = s + sizeof(R->rl_buf) - 1;
			while(s < se && (x = getc(nl)) >= ' ')
				*++s = x;
			*s = 0;
			}
		fprintf(Stderr, ": %s\n", R->rl_buf);
		}
	else
		fprintf(Stderr, "\n");
	exit_ASL(R,1);
	}

#undef asl

#define Mb_gulp 31
 typedef struct Mblock {
	struct Mblock *next;
	void *m[Mb_gulp];
	} Mblock;

 void **
M1record_ASL(Edaginfo *I, void *x)
{
	Mblock *mb;
	void **rv;

	if (I->Mbnext >= I->Mblast) {
		mb = (Mblock *)Malloc(sizeof(Mblock));
		mb->next = (Mblock*)I->Mb;
		I->Mb = (void*)mb;
		I->Mbnext = mb->m;
		I->Mblast = mb->m + Mb_gulp;
		}
	rv = I->Mbnext++;
	*rv = x;
	return rv;
	}

 void *
M1alloc_ASL(Edaginfo *I, size_t n)
{
	Mblock *mb;

	if (I->Mbnext >= I->Mblast) {
		mb = (Mblock *)Malloc(sizeof(Mblock));
		mb->next = (Mblock*)I->Mb;
		I->Mb = (void*)mb;
		I->Mbnext = mb->m;
		I->Mblast = mb->m + Mb_gulp;
		}
	return *I->Mbnext++ = Malloc(n);
	}

 void *
M1zapalloc_ASL(Edaginfo *I, size_t n)
{
	void *rv;

	memset(rv = M1alloc_ASL(I, n), 0, n);
	return rv;
	}

 void
M1free_ASL(Edaginfo *I, void **mnext, void **mlast)
{
	void **x, **x0;
	Mblock *Mb, *mb;
	void **Mblast;

	if (!(Mb = (Mblock *)I->Mb))
		return;
	x = (void **)I->Mbnext;
	Mblast = I->Mblast;
	I->Mbnext = mnext;
	I->Mblast = mlast;
	x0 = Mb->m;
	for(;;) {
		if (mlast == Mblast)
			x0 = mnext;
		while(x > x0)
			if (*--x)
				free(*x);
		if (mlast == Mblast) {
			I->Mb = (void*)Mb;
			return;
			}
		mb = Mb->next;
		free(Mb);
		if (!(Mb = mb))
			break;
		x0 = Mb->m;
		Mblast = x = x0 + Mb_gulp;
		}
	I->Mb = 0;
	}

 void
xknown_(real *x)
{
	ASL *asl;
	if (!(asl = cur_ASL))
		badasl_ASL(asl,0,"xknown");
	xknowne(x, (fint*)0);
	}

 void
xknowe_(real *x, fint *nerror)
{
	ASL *asl;
	if (!(asl = cur_ASL))
		badasl_ASL(asl,0,"xknown");
	xknowne(x, nerror);
	}


 void
xunkno_(void)
{
	ASL *asl;
	if (!(asl = cur_ASL))
		badasl_ASL(asl,0,"xunkno");
	asl->i.x_known = 0;
	}

 void
mnnzchk_ASL(ASL *asl, fint *M, fint *N, size_t NZ, const char *who1)
{
	int n;
	if (!asl)
		goto bad;
	n = asl->i.ASLtype;
	if (n < ASL_read_fg || n > ASL_read_pfgh)
		goto bad;
	ASL_CHECK(asl, n, who1);
	if (*M == n_con && *N == c_vars && NZ == nzjac)
		return;
	what_prog();
	fprintf(Stderr,
 "%s: got M = %ld, N = %ld, NZ = %ld\nexpected M = %d, N = %d, NZ = %d\n",
			who1, (long)*M, (long)*N, NZ, n_con, c_vars, nzjac);
	exit(1);
 bad:
	badasl_ASL(asl, ASL_read_fg, who1);
	}

 void
LUcopy_ASL(int nv, real *L, real *U, real *LU)
{
	real *LUe;
	for(LUe = LU + 2*nv; LU < LUe; LU += 2) {
		*L++ = LU[0];
		*U++ = LU[1];
		}
	}

 int
already_ASL(const char *who)
{
	fprintf(Stderr, "%s called after ASL_alloc().\n", who);
	return 1;
	}

 void
ASL_free(ASL **aslp)
{
	ASL *a;
	ASLhead *h;
	extern void at_end_ASL ANSI((Exitcall*));

	if (!(a = *aslp))
		return;
	if (a == cur_ASL)
		cur_ASL = 0;
	h = a->p.h.prev;
	(h->next = a->p.h.next)->prev = h;
	if (a->i.arprev)
		at_end_ASL(a->i.arprev);
	M1free(&a->i, (void**)0, (void**)0);
	free((void*)a);
	*aslp = 0;
	}

 void
badasl_ASL(ASL *a, int n, const char *who)
{
	if (!Stderr)
		Stderr_init_ASL();	/* set Stderr if necessary */
	if (a)
		fprintf(Stderr,
			"\n*** %s needs ASL_alloc(%d), not ASL_alloc(%d)\n",
			who, n, a->i.ASLtype);
	else if (n)
		fprintf(Stderr, "\n*** %s called before ASL_alloc(%d)\n",
			who, n);
	else
		fprintf(Stderr,
		 "\n*** %s called before ASL_alloc, jacdim, jac2dim, or jacpdim\n", who);
	exit(1);
	}

#define SKIP_NL2_DEFINES
#include "nlp.h"
#include "nlp2.h"
#include "asl_pfg.h"
#include "asl_pfgh.h"

 ASL *
ASL_alloc(int k)
{
	static int msize[5] = {
		sizeof(ASL_fg),
		sizeof(ASL_fg),
		sizeof(ASL_fgh),
		sizeof(ASL_pfg),
		sizeof(ASL_pfgh)
		};
	ASL *a;
	ASLhead *h;
	int n;

	if (!Stderr)
		Stderr_init_ASL();	/* set Stderr if necessary */
	Mach_ASL();
	if (k < 1 || k > 5)
		return 0;
	a = (ASL*) mymalloc(n = msize[k-1]);
	memcpy(a, &edagpars_ASL, sizeof(Edagpars));
	memset(&a->i, 0, n - sizeof(Edagpars));
	a->i.ASLtype = k;
	a->i.n_prob = 1;
	switch(k) {
	  case ASL_read_pfg:	((ASL_pfg*)a)->P.merge = 1; break;
	  case ASL_read_pfgh:	((ASL_pfgh*)a)->P.merge = 1;
	  }
	h = a->p.h.next = ASLhead_ASL.next;
	a->p.h.prev = h->prev;
	h->prev = ASLhead_ASL.next = &a->p.h;
	return cur_ASL = a;
	}

#define Egulp 400

 void *
mem_ASL(ASL *asl, unsigned int len)
{
	fint k;
	char *memNext;

	if (len >= 256)
		return M1alloc(len);
#ifdef Double_Align
	len = (len + (sizeof(real)-1)) & ~(sizeof(real)-1);
#else
	len = (len + (sizeof(int)-1)) & ~(sizeof(int)-1);
#endif
	ACQUIRE_DTOA_LOCK(MEM_LOCK);
	memNext = asl->i.memNext;
	if (memNext + len >= asl->i.memLast) {
		memNext = (char *)M1alloc(k = Egulp*Sizeof(expr) + len);
		asl->i.memLast = memNext + k;
		}
	asl->i.memNext = memNext + len;
	FREE_DTOA_LOCK(MEM_LOCK);
	return memNext;
	}

 EdRead *
EdReadInit_ASL(EdRead *R, ASL *asl, FILE *nl, void *S)
{
	R->asl = asl;
	R->nl = nl;
	R->S = S;
	R->Line = 10;
	R->lineinc = 1;
	R->can_end = 0;
	R->dadjfcn = asl->i.dadjfcn;
	R->iadjfcn = asl->i.iadjfcn;
	return R;
	}

 void
Suf_read_ASL(EdRead *R, int readall)
{
	ASL *asl = R->asl;
	SufDesc *D;
	char *s, sufname[128];
	const char *fmt;
	int *d, i, isreal, k, n, nx, nx1;
	real *r, t;

	if (xscanf(R, "%d %d %127s", &k, &n, sufname) != 3)
		badline(R);
	if (k < 0 || k > 7 || n <= 0)
		badline(R);
	isreal = k & ASL_Sufkind_real;
	k &= ASL_Sufkind_mask;
	nx = (&asl->i.n_var_)[k];
	if (k == 1)
		nx += n_lcon;
	if (n > nx)
		badline(R);
	if (readall) {
 new_D:
		D = (SufDesc*)M1zapalloc(sizeof(SufDesc) + strlen(sufname) + 1);
		D->next = asl->i.suffixes[k];
		asl->i.suffixes[k] = D;
		asl->i.nsuff[k]++;
		asl->i.nsuffixes++;
		strcpy(s = (char*)(D+1), sufname);
		D->sufname = s;
		D->kind = k;
		if (isreal)
			D->kind |= ASL_Sufkind_real;
		}
	else for(D = asl->i.suffixes[k]; ; D = D->next) {
		if (!D) {
			if (readall)
				goto new_D;
 skip:
			/* Skip this suffix table */
			fmt = (char*)(isreal ? "%d %lf" : "%d %d");
			do if (xscanf(R,fmt,&k,&t) != 2)
					badline(R);
				while(--n);
			return;
			}
		if (k == (D->kind & ASL_Sufkind_mask)
		 && !strcmp(sufname,D->sufname))
			break;
		}
	if ((D->kind & ASL_Sufkind_outonly) == ASL_Sufkind_outonly)
		goto skip;
	nx1 = nx + D->nextra + asl->i.nsufext[k];
	if (D->kind & ASL_Sufkind_real) {
		D->u.i = 0;
		if (!(r = D->u.r))
			D->u.r = r = (real*)mem(nx1*sizeof(real));
		if (n < nx)
			memset(r,0,nx*sizeof(real));
		if (nx < nx1)
			memset(r+nx, 0, (nx1-nx)*sizeof(real));
		if (isreal)
			do  {
				if (xscanf(R,"%d %lf",&i,&t) != 2
				 || i < 0 || i >= nx)
					badline(R);
				r[i] = t;
				}
				while(--n);
		else
			do  {
				if (xscanf(R,"%d %d",&i,&k) != 2
				 || i < 0 || i >= nx)
					badline(R);
				r[i] = k;
				}
				while(--n);
		}
	else {
		D->u.r = 0;
		if (!(d = D->u.i))
			D->u.i = d = (int*)mem(nx1*sizeof(int));
		if (n < nx)
			memset(d,0,nx*sizeof(int));
		if (nx < nx1)
			memset(d+nx, 0, (nx1-nx)*sizeof(int));
		if (isreal)
			do {
				if (xscanf(R,"%d %lf",&i,&t) != 2
				 || i < 0 || i >= nx)
					badline(R);
				d[i] = (int)(t + 0.5);
				} while(--n);
		else
			do {
				if (xscanf(R,"%d %d",&i,&k) != 2
				 || i < 0 || i >= nx)
					badline(R);
				d[i] = k;
				}
				while(--n);
		}
	D->kind |= ASL_Sufkind_input;
	}

 real
f_OPNUM_ASL(expr_n *e)
{
#ifdef _WIN32	/* Work around a Microsoft linker bug... */
		/* Without the following test, f_OPNUM gets confused */
		/* with f_OPVARVAL.  Both get mapped to the same address */
		/* in the r_ops_ASL array defined in fg_read.c. */
	if (!e) {
		printf("f_OPNUM(e) has e = 0\n");
		return 0.;
		}
#endif
	return e->v;
	}

 void
No_derivs_ASL(const char *who)
{
	fprintf(Stderr, "\nBUG: %s called with want_derivs == 0.\n", who);
	exit(1);
	}

#define ndcc asl->i.ndcc_
#define nzlb asl->i.nzlb_

 void
flagsave_ASL(ASL *asl, int flags)
{
	real t;
	ssize_t nc, nv, nz;

	t = nZc;
	if (t >= 2147483648. /* 2^31 */) {
		if (sizeof(size_t) <= 4) {
			fprintf(Stderr, "\n*** Problem too large for 32-bit "
					"addressing (%.g Jacobian nonzeros).\n", t);
			exit(1);
			}
		else if (!(flags & (ASL_allow_Z | ASL_use_Z))) {
			fprintf(Stderr, "\n*** Problem too large (%.g Jacobian nonzeros)\n", t);
			exit(1);
			}
		else if (sizeof(((cgrad*)0)->goff) <= 4)
			fprintf(Stderr, "\n*** Problem too large (%.g Jacobian nonzeros) for "
				"jacval().\nRecompile ASL with \"#define ASL_big_goff\" "
				"added to arith.h.\n", t);
		flags |= ASL_use_Z;
		}
	asl->i.rflags = flags;
	if (flags & ASL_cc_simplify && n_cc) {
		if (ndcc < 0)
			/* supply overestimates */
			ndcc = nzlb = n_cc;
		asl->i.nsufext[ASL_Sufkind_var] += 3*ndcc + n_cc + nzlb;
		asl->i.nsufext[ASL_Sufkind_con] += 2*ndcc + nzlb;
		/* use nsufext[ASL_Sufkind_prob] for # of extra Jacobian nonzeros */
		asl->i.nsufext[ASL_Sufkind_prob] += 5*ndcc + n_cc + 2*nzlb;
		}
	nv = n_var + asl->i.nsufext[ASL_Sufkind_var];
	nc = n_con + asl->i.nsufext[ASL_Sufkind_con];
	nz = nZc + asl->i.nsufext[ASL_Sufkind_prob];
	if (!LUv) {
		LUv = (real*)M1alloc(2*sizeof(real)*nv);
		if (flags & ASL_sep_U_arrays)
			Uvx = LUv + nv;
		}
	if (!LUrhs) {
		LUrhs = (real*)M1alloc(2*sizeof(real)*nc);
		if (flags & ASL_sep_U_arrays)
			Urhsx = LUrhs + nc;
		}
	if (flags & ASL_sep_U_arrays) {
		if (!Uvx)
			Uvx = (real*)M1alloc(nv*sizeof(real));
		if (!Urhsx)
			Urhsx = (real*)M1alloc(nc*sizeof(real));
		}
	if (flags & ASL_want_A_vals && !A_vals)
		A_vals = (real*)M1alloc(nz*sizeof(real));
	if (A_vals) {
		if (!A_rownos)
			A_rownos = (int *)M1alloc(nz*sizeof(int));
		}
	else if (nc)
		asl->i.Cgrad0 = asl->i.Cgrad_ = (cgrad **)M1zapalloc(nc*sizeof(cgrad *));
	}

 static void
zerograd_chk(ASL *asl)
{
	int j, n, nv, nx, *z, **zg;
	ograd *og, **ogp, **ogpe;

	nx =  asl->i.nsufext[ASL_Sufkind_var];
	if (!(nv = asl->i.nlvog)) {
		nv = n_var;
		if (nv > asl->i.n_var0)
			nx -= nv - asl->i.n_var0;
		}
	zerograds = 0;
	ogp = Ograd;
	ogpe = ogp + (j = n_obj);
	while(ogp < ogpe) {
		og = *ogp++;
		n = 0;
		while(og) {
			j += og->varno - n;
			n = og->varno + 1;
			if (n >= nv)
				break;
			og = og->next;
			}
		if (n < nv)
			j += nv - n;
		}
	if (j == n_obj)
		return;
	j += n_obj * nx;
	zerograds = zg = (int **)mem(n_obj*sizeof(int*)+j*sizeof(int));
	z = (int*)(zg + n_obj);
	ogp = Ograd;
	while(ogp < ogpe) {
		*zg++ = z;
		og = *ogp++;
		n = 0;
		while(og) {
			while(n < og->varno)
				*z++ = n++;
			og = og->next;
			if (++n >= nv)
				break;
			}
		while(n < nv)
			*z++ = n++;
		*z++ = -1;
		z += nx;
		}
	}

 void
adjust_zerograds_ASL(ASL *asl, int nnv)
{
	int i, j, k, n, *z, **zg, **zge;

	if (!(zg = zerograds)) {
		zerograd_chk(asl);
		return;
		}
	n = n_var;
	for(zge = zg + n_obj; zg < zge; ++zg) {
		z = *zg;
		for(i = 0; (k = z[i]) >= 0; ++i) {
			if (k >= n)
				break;
			}
		for(j = n, k = nnv; k > 0; --k)
			z[i++] = j++;
		z[i] = -1;
		}
	}


#ifndef NO_BOUNDSFILE_OPTION /*{*/
 extern void bswap_ASL(void *, size_t);

 static void
bad_bounds(ASL *asl, const char *fmt, ...)
{
	Jmp_buf *J;
	va_list ap;
	va_start(ap, fmt);
	if (progname)
		fprintf(Stderr, "\n%s: ", progname);
	else
		fprintf(Stderr, "\n");
	vfprintf(Stderr, fmt, ap);
	fprintf(Stderr, "\n");
	va_end(ap);
	if ((J = asl->i.err_jmp_))
		longjmp(J->jb, 1);
	exit(1);
	}
#endif /*}*/

 int
prob_adj_ASL(ASL *asl)
{
	cgrad *cg, **pcg, **pcge;
	int flags, k;
#ifndef NO_BOUNDSFILE_OPTION /*{*/
	FILE *f;
	char *bf, buf[4096], *s, *s1, *se;
	int a, i, n1, nr, nv, swap;
	real *L, *U, *x;
	size_t inc, m, n;
#endif /*}*/

	if (n_obj)
		adjust_zerograds_ASL(asl, 0);
	flags = asl->i.rflags;
	asl->i.Cgrad0 = asl->i.Cgrad_;
	if (flags & (ASL_obj_replace_eq | ASL_obj_replace_ineq))
		obj_adj_ASL(asl);
	if (!A_vals) {
		if (flags & ASL_cc_simplify && n_cc)
			mpec_adjust_ASL(asl);
		if (flags & ASL_rowwise_jac) {
			pcg = Cgrad;
			pcge = pcg + n_con;
			k = 0;
			while(pcg < pcge)
				for(cg = *pcg++; cg; cg = cg->next)
					cg->goff = k++;
			}
		}
	if (n_obj)
		zerograd_chk(asl);
#ifndef NO_BOUNDSFILE_OPTION /*{*/
	if ((bf = asl->i.boundsfile)) {
		if (!(f = fopen(bf, "rb")))
			bad_bounds(asl, "Cannot open boundsfile \"%s\".", bf);
		m = sizeof(buf);
		if ((n = fread(buf, 1, m, f)) < 24
		 || strncmp(buf, "Bounds, x; arith ", 17)) {
 badmagic:
			bad_bounds(asl, "Bad magic in boundsfile \"%s\".", bf);
			}
		a = (int)strtol(s = buf+17, &se, 10);
		if (se <= s || a < 0 || a > 2 || *se >= ' ')
			goto badmagic;
		if (a == 0 && *se == '\r')
			++se;
		if (*se++ != '\n')
			goto badmagic;
		nv = n_var;
		L = LUv;
		if ((U = Uvx))
			inc = 1;
		else {
			U = L + 1;
			inc = 2;
			}
		if (!(x = X0))
			X0 = x = (real*)M1alloc(nv*sizeof(real));
		swap = 0;
		if (a) {
			if (a != Arith_Kind_ASL)
				swap = 1;
			s = buf + 20;
			n1 = *(int*)s;
			if (swap)
				bswap_ASL(&n1, sizeof(n1));
			if (n1 != nv) {
 bad_n1:
				bad_bounds(asl, "Expected %d bounds triples in boundsfile"
					" \"%s\"; got %d.", nv, bf, n1);
				}
			s += sizeof(int);
			se = buf + n;
			for(i = 1; i <= nv; ++i) {
				if (se-s < 3*sizeof(real)) {
					for(s1 = buf; s < se; ++s, ++s1)
						*s1 = *s;
					m = s1 - buf;
					m  += n = fread(s1, 1, sizeof(buf) - m, f);
					se = buf + m;
					if (m < 3*sizeof(real)) {
 too_few:
						bad_bounds(asl, "%d too few bound triples "
							"in boundsfile \"%s\".",
							nv - i + 1, bf);
						}
					s = buf;
					}
				*L = *(real*)s;
				s += sizeof(real);
				*x = *(real*)s;
				s += sizeof(real);
				*U = *(real*)s;
				s += sizeof(real);
				if (swap) {
					bswap_ASL(L, sizeof(real));
					bswap_ASL(x, sizeof(real));
					bswap_ASL(U, sizeof(real));
					}
				if (*L > *U || *x < *L || *x > *U) {
 bad_triple:
					bad_bounds(asl, "bad bound triple %d in bounds file "
						"\"%s\":\n\tL = %.g\n\tx = %.g\n\tU = %.g",
						i, bf, *L, *x, *U);
					}
				L += inc;
				U += inc;
				++x;
				}
			}
		else {
			n1 = (int)strtol(se, &s, 10);
			if (n1 != nv)
				goto bad_n1;
			se = buf + n;
			for(i = 1; i <= nv; ++i) {
				nr = 0;
 tryagain:
				while(s < se && *s <= ' ')
					++s;
				s1 = s;
				for(k = 0; k < 3 && s1 < se; ++k) {
					while(s1 < se && *s1 > ' ')
						++s1;
					while(s1 < se && *s1 <= ' ')
						++s1;
					}
				if (k < 3) {
					if (nr)
						goto too_few;
					for(s1 = buf; s < se; ++s, ++s1)
						*s1 = *s;
					m = sizeof(buf) - (s1 - buf);
					n = fread(s1, 1, m, f);
					se = s1 + n;
					s = buf;
					nr = 1;
					goto tryagain;
					}
				*L = strtod(s, &s1);
				if (s1 <= s || *s1 > ' ') {
 badnumber:
					while(s1 < se && *s1 > ' ')
						++s1;
					bad_bounds(asl, "Bound triple %d: bad number \"%.*s\""
						" in boundsfile \"%s\".", i,
						(int)(s1-s), s, bf);
					}
				*x = strtod(s = s1, &s1);
				if (s1 <= s || *s1 > ' ')
					goto badnumber;
				*U = strtod(s = s1, &s1);
				if (s1 <= s || *s1 > ' ')
					goto badnumber;
				if (*L > *U || *x < *L || *x > *U)
					goto bad_triple;
				L += inc;
				U += inc;
				++x;
				s = s1;
				}
			}
		fclose(f);
		}
#endif /*}*/
	asl->i.err_jmp_ = 0;
	return 0;
	}

 void
suf_declare_ASL(ASL *asl, SufDecl *sd, int n)
{
	SufDesc *d, *dnext[4];
	SufDecl *sde;
	int i, j;

	if (!asl)
		badasl_ASL(asl, 0, "suf_declare");
	asl->i.nsuffixes = 0;
	if (n > 0) {
		asl->i.nsuffixes = n;
		d = (SufDesc*)M1alloc(n*sizeof(SufDesc));
		memset(asl->i.nsuff, 0, 4*sizeof(int));
		for(i = 0; i < n; i++)
			asl->i.nsuff[sd[i].kind & ASL_Sufkind_mask]++;
		for(i = 0; i < 4; i++)
			if ((j = asl->i.nsuff[i]))
				asl->i.suffixes[i] = d += j;
		memset(dnext, 0, 4*sizeof(SufDesc*));
		for(sde = sd + n; sd < sde; sd++) {
			d = --asl->i.suffixes[i = sd->kind & ASL_Sufkind_mask];
			d->next = dnext[i];
			dnext[i] = d;
			d->sufname = sd->name;
			d->table = sd->table;
			d->kind = sd->kind & ~ASL_Sufkind_input;
			d->nextra = sd->nextra;
			d->u.i = 0;
			d->u.r = 0;
			}
		}
	}

 SufDesc *
suf_get_ASL(ASL *asl, const char *name, int kind)
{
	SufDesc *d, *de;
	int ifread;

	if (!asl)
		badasl_ASL(asl, 0, "suf_get");
	ifread = kind & ASL_Sufkind_input;
	d = asl->i.suffixes[kind &= ASL_Sufkind_mask];
	de = d + asl->i.nsuff[kind];
	for(;; d++) {
		if (d >= de) {
			fprintf(Stderr, "suf_get(\"%s\") fails!\n", name);
			exit(1);
			}
		if (!strcmp(name, d->sufname))
			break;
		}
	if (ifread && !(d->kind & ASL_Sufkind_input))
		d = 0;
	return d;
	}

 SufDesc *
suf_iput_ASL(ASL *asl, const char *name, int kind, int *I)
{
	SufDesc *d = suf_get_ASL(asl, name, kind);
	d->u.i = I;
	d->kind &= ~ASL_Sufkind_real;
	d->kind |= ASL_Sufkind_output;
	return d;
	}

 SufDesc *
suf_rput_ASL(ASL *asl, const char *name, int kind, real *R)
{
	SufDesc *d = suf_get_ASL(asl, name, kind);
	d->u.r = R;
	d->kind |= ASL_Sufkind_output | ASL_Sufkind_real;
	return d;
	}

 int *
get_vcmap_ASL(ASL *asl, int k)
{
	cgrad **cgp;
	int i, m, n, *x;

	if ((x = (&asl->i.vmap)[k &= 1]))
		return x;
	m = 0;
	if (k == ASL_Sufkind_con && Cgrad)
		m = asl->i.n_con0 + asl->i.nsufext[ASL_Sufkind_con];
	n = (&asl->i.n_var0)[k] + asl->i.nsufext[k];
	cgp = (cgrad**)M1alloc(m * sizeof(cgrad*) + n*sizeof(int));
	x = (&asl->i.vmap)[k] = (int*)(cgp + m);
	for(i = 0; i < n; ++i)
		x[i] = i;
	asl->p.Conival = conivalmap;
	asl->p.Congrd = congrdmap;
	if (m)
		memcpy(asl->i.Cgrad0 = cgp, Cgrad, m*sizeof(cgrad*));
	return x;
	}

 int *
get_vminv_ASL(ASL *asl)
{
	int i, j, n, *vm, *x;

	if ((x = asl->i.vminv))
		return x;
	if (!(vm = asl->i.vmap))
		vm = get_vcmap_ASL(asl, ASL_Sufkind_var);
	n = asl->i.n_var0 + asl->i.nsufext[ASL_Sufkind_var];
	x = (int*)M1alloc(n*sizeof(int));
	for(i = 0; i < n; ++i)
		x[i] = -1;
	n = n_var;
	for(i = 0; i < n; ++i) {
		if ((j = vm[i]) >= 0)
			x[j] = i;
		}
	for(i = 0, j = n; i < n; ++i)
		if (x[i] < 0)
			x[i] = j++;
	return asl->i.vminv = x;
	}

 int
ka_read_ASL(ASL *asl, EdRead *R, int mode, int **kap, size_t **kapZ)
{
	int flags, *kai;
	size_t i, k, *ka, t;
	int j;
	unsigned Long u;

	k = asl->i.n_var0;
	if (!xscanf(R,"%d",&j) || j != k - 1)
		return 1;
	if ((i = k) < n_var)
		i = n_var;
	flags = asl->i.rflags;
	if (flags & ASL_use_Z) {
		*kap = kai = A_colstarts = 0;
		if (!(ka = A_colstartsZ))
			A_colstartsZ = ka = (size_t*)M1alloc((i+1)*Sizeof(size_t));
		*kapZ = ka + 1;
		}
	else {
		*kapZ = ka = A_colstartsZ = 0;
		if (!(kai = A_colstarts))
			A_colstarts = kai = (int*)M1alloc((i+1)*Sizeof(int));
		*kap = kai + 1;
		}
	if (sizeof(int) == sizeof(size_t)) {
		if (!ka)
			ka = (size_t*)kai;
		*ka++ = 0;
		*ka++ = 0;	/* sic */
		if (mode == 'K') {
			t = 0;
			while(--k > 0) {
				if (!xscanf(R, "%d", &u))
					return 1;
				*ka++ = t += u;
				}
			}
		else {
			while(--k > 0) {
				if (!xscanf(R, "%d", &u))
					return 1;
				*ka++ = u;
				}
			}
		}
	else if (flags & ASL_use_Z) {
		*ka++ = 0;
		*ka++ = 0;	/* sic */
		if (mode == 'K') {
			t = 0;
			while(--k > 0) {
				if (!xscanf(R, "%d", &u))
					return 1;
				*ka++ = t += u;
				}
			}
		else {
			while(--k > 0) {
				if (!xscanf(R, "%d", &u))
					return 1;
				*ka++ = u;
				}
			}
		}
	else {
		*kai++ = 0;
		*kai++ = 0;	/* sic */
		if (mode == 'K') {
			t = 0;
			while(--k > 0) {
				if (!xscanf(R, "%d", &u))
					return 1;
				*kai++ = (int)(t += u);
				}
			}
		else {
			while(--k > 0) {
				if (!xscanf(R, "%d", &u))
					return 1;
				*kai++ = (int)u;
				}
			}
		}
	return 0;
	}

 void
goff_comp_ASL(ASL *asl)
{
	cgrad *cg, **cgx, **cgxe;
	size_t *ka;


	/* The following "if" test and either its "then" or its "else" block */
	/* should be optimized away. */

	if (sizeof(size_t) == sizeof(int)) {
		cgx = Cgrad;
		cgxe = cgx + asl->i.n_con0;
		if (!(ka = A_colstartsZ))
			ka = (size_t*)A_colstarts;
		++ka;
		while(cgx < cgxe)
			for(cg = *cgx++; cg; cg = cg->next)
				cg->goff = ka[cg->varno]++;
		}
	else {
		int *kai;

		cgx = Cgrad;
		cgxe = cgx + asl->i.n_con0;
		if ((kai = A_colstarts)) {
			++kai;
			while(cgx < cgxe)
				for(cg = *cgx++; cg; cg = cg->next)
					cg->goff = kai[cg->varno]++;
			}
		else {
			ka = A_colstartsZ + 1;
			while(cgx < cgxe)
				for(cg = *cgx++; cg; cg = cg->next)
					cg->goff = ka[cg->varno]++;
			}
		}
	}

 void
colstart_inc_ASL(ASL *asl)
{
	size_t *ka, *kae;

	ka = A_colstartsZ;
	if (sizeof(size_t) == sizeof(int)) {
		if (!ka)
			ka = (size_t*)A_colstarts;
		kae = ka + asl->i.n_var0;
		while(ka <= kae)
			++*ka++;
		}
	else if (ka) {
		kae = ka + asl->i.n_var0;
		while(ka <= kae)
			++*ka++;
		}
	else {
		int *kai, *kaie;

		kai = A_colstarts;
		kaie = kai + asl->i.n_var0;
		while(kai <= kaie)
			++*kai++;
		}
	}

#ifdef __cplusplus
}
#endif
