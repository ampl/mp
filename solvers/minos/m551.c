/****************************************************************
Copyright (C) 1997-1998 Lucent Technologies
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

#include "signal.h"
#include "ctype.h"
#include "getstub.h"

#define N_FNAMES 100
static char *fnames[N_FNAMES];

#ifdef __cplusplus
extern "C" {
#else
#endif

 static fint
m1user(fint *phase, fint *ierror, fint *majits, fint *minits, fint *ninf, real *sinf,
	char *label, real *step, real *obj, real *prinf, real *duinf, fint *nssave,
	fint *nfcon, real *penpar, fint *nswap)
{
	return 0;
	}

 typedef fint (*m1user_t)(fint *phase, fint *ierror, fint *majits, fint *minits,
			  fint *ninf, real *sinf, char *label, real *step, real *obj,
			  real *prinf, real *duinf, fint *nssave, fint *nfcon,
			  real *penpar, fint *nswap);

 typedef void (*opkey_t)(char*, char*, fint*, fint*, fint*, ftnlen, ftnlen);

 extern void m1init_(void);
 extern void m1time_(fint *clock, fint *prtopt);
 extern void m2core_(fint *mode, fint *mincor);
 extern void m3dflt_(fint *mode);
 extern void m3file_(fint *ncalls, fint *ioptns, opkey_t opkey,
	fint *iprint, fint *isumm, fint *inform);
 extern void m3inpt_(real *objadd, real *z, fint *nwcore);
 extern void m3key_(char *buffer, char *key, fint *lprnt, fint *lsumm,
	fint *inform, ftnlen buffer_len, ftnlen key_len);
 extern void mifile_(fint *mode);
 extern int minoss_(char *start, fint *m, fint *n, fint *nb, fint *ne,
			fint *nname, fint *nncon, fint *nnobj, fint *nnjac,
			fint *iobj, real *objadd, char *names, real *a,
			fint *ha, fint *ka, real *bl, real *bu, fint *
			name1, fint *name2, fint *hs, real *xn, real *pi,
			real *rc, fint *inform, fint *mincor, fint *ns,
			fint *ninf, real *sinf, real *obj, real *z,
			fint *nwcore, ftnlen start_len, ftnlen names_len);
/*!! in misolv_(), hs and ha are now int* rather than shortint*; */
/*!! args ifuser, m1user, majitn, minitn added */
 extern void misolv_(fint *mimode, char *start, fint *mxx, fint *nxx,
	fint *nbxx, fint *nexx, fint *nkax, fint *nnamex, fint *iobjxx,
	real *objadd, real *a, fint *ha, fint *ka, real *bl, real *bu,
	fint *name1, fint *name2, fint *hs, real *xn, real *pi, real *rc,
	fint *inform, fint *ns, real *z, fint *nwcore, fint *ifuser,
	m1user_t m1user, fint *majitn, fint *minitn, ftnlen start_len);
 extern void mistart_(fint*, fint*, fint*);

#ifdef NO_F2C
 char **xargv;
 extern void MAIN__(void);

 int
main(int argc, char **argv)
{
	xargv = argv;
	MAIN__();
	return 0;
	}
#else /* !NO_F2C */
extern char **xargv;
#endif

 static int functimes, timing, wantfuncs;
 static int outlev = 1;
 static real meminc, memused, memused0, memused1, tsolve;
 static char xxxvers[] = "AMPL/MINOS 5.51\0\nAMPL/MINOS 5.51 Driver Version 20090923\n";
#define asl cur_ASL

 extern struct { fint jverif[4], lverif[2]; } m8veri_;

 extern struct {
	fint alone, ampl, gams, mint, page1, page2;
	} m1env_;

 extern struct {
	real eps, eps0, eps1, eps2, eps3, eps4, eps5, plinfy;
	} m1eps_;

 extern struct {
	fint iread, iprint, isumm;
	} m1file_;

 extern struct {
	real tlast[5], tsum[5];
	fint numt[5], ltime;
	} m1tim_;

 extern struct {
	fint mrows, mcols, melms;
	} m2len_;

 extern struct {
	fint minlu, maxlu, lena, nbelem, ip, iq, lenc, lenr, locc, locr, iploc,
		iqloc, lua, indc, indr;
	} m2lu1_;

 extern struct {
	real parmlu[30];
	fint luparm[30];
	} m2lu4_;

 extern struct {
	fint maxw, maxz;
	} m2mapz_;

 extern struct {
	real aijtol, bstruc[2];
	fint mlst, mer;
	real aijmin, aijmax;
	fint na0, line, ier[20];
	} m3mps3_;

 extern struct {
	real sclobj, scltol;
	fint lscale;
	} m3scal_;

 extern struct {
	fint kchk, kinv, ksav, klog, ksumm, i1freq, i2freq, msoln;
	} m5freq_;

 extern struct {
	fint maxr, maxs, mbs, nn, nn0, nr, nx;
	} m5len_;

 extern struct {
	real sinf, wtobj;
	fint minimz, ninf, iobj, jobj, kobj;
	} m5lobj_;

 extern struct {
	fint prnt0, prnt1, summ0, summ1, newhed;
	} m5log4_;

 extern struct {
	fint itn, itnlim, nphs, kmodlu, kmodpi;
	} m5lp1_;

 extern struct {
	real fobj, fobj2;
	fint nnobj, nnobj0;
	} m7len_;

 extern struct {
	fint nfcon[4], nfobj[4], nprob, nstat1, nstat2;
	} m8func_;

 extern struct {
	real penpar, rowtol;
	fint ncom, nden, nlag, nmajor, nminor;
	} m8al1_;

 extern struct {
	fint njac, nncon, nncon0, nnjac;
	} m8len_;

#define objsen m5lobj_.minimz

 static char *usage_msg[] = {
 "  where  stub  is from  `ampl -obstub`  or  `ampl -ogstub`.  Assignments",
 "  have the form  u=filename  or  spec=value  in which the integers u are",
 "  file unit numbers that appear in spec assignments (don't use 5 or 6",
 "  for u) and  spec  is a MINOS SPECS file keyword phrase, with keywords",
 "  in the phrase separated by _ (underscore).  Use",
 "	outlev=0 for no options echoed on stdout,",
 "	outlev=1 for neither log nor summary file on stdout (default),",
 "	outlev=2 to see summary output on stdout,",
 "	outlev=3 to see detailed (log file) output on stdout, and",
 "	outlev=4 to get log file plus solution on stdout.",
 "  For outlev <= 2, use 7=logfilename to get the log file.  Assignments",
#ifdef Student_Edition
 "  may also be given in $minos_options .",
#else
 "  may also be given in $minos_options .  No stub ==> read SPECS file on",
 "  stdin, followed (unless the SPECS file says otherwise) by an MPS file.",
#endif
		0 };

#ifdef Student_Edition
#ifndef NO_MINOS1
#define NO_MINOS1
#endif
#endif

#ifndef NO_MINOS1
 /* stuff derived from mi00main.f in the MINOS source */
 extern struct {
	fint iback, idump, iload, imps, inewb, insrt, ioldb, ipnch, iprob,
		iscr, isoln, ispecs, ireprt;
	} m2file_;

 extern struct {
	fint ne, nka, la, lha, lka;
	} m2mapa_;

 extern struct {
	fint m, n, nb, nscl;
	} m3len_;

 extern struct {
	fint lascal, lbl, lbu, lbbl, lbbu, lhrtyp, lhs, lkb;
	} m3loc_;

 extern struct {
	fint lname1, lname2, lkeynm, nname;
	} m3mps1_;

 extern struct {
	fint lpi, lpi2, lw, lw2, lx, lx2, ly, ly2, lgsub, lgsub2, lgrd, lgrd2,
		lr, lrg, lrg2, lxn;
	} m5loc_;

 extern struct {
	fint idebug, ierr, lprint;
	} m5log1_;

 static void
minos1_(real *z, fint *nz)
{
	fint i, inform, j, k, k1, lrc, majitn, mimode, mincor, minitn;
	fint ifuser, ncalls, ns, nwcore, nz1;
	real objadd, *z0, *z1;
	static fint I0 = 0, I1 = 1, IM1 = -1, I2 = 2, I3 = 3;

	m2file_.ispecs = 5;
	m1file_.isumm = 6;
	m1file_.iprint = 0;

	mistart_(&m1file_.iprint, &m1file_.isumm, &m2file_.ispecs);
	z0 = z;
	nz1 = *nz;

	for(i = 1;; i++) {
		ncalls = i;
		m1tim_.ltime = 1;
		m1time_(&I0, &I0);
		inform = 0;
		m3file_(&ncalls, &m2file_.ispecs, m3key_, &m1file_.iprint,
			&m1file_.isumm, &inform);
		if (inform >= 2)
			return;
		m3dflt_(&I2);
		mifile_(&I2);
		if ((j = m8len_.nncon) < 10)
			j = 10;
		m8len_.njac = j * m8len_.nnjac;
		m2core_(&I1, &mincor);
		if (meminc > 0.) {
			k = (1 << 17) * meminc; /* 17 because sizeof(real) == 8 == 2^3 */
			k1 = mincor + k;
			if (mincor >= k1) {
				fprintf(Stderr,
				 "meminc = %g megabytes is too big a memory increment.\n",
				 meminc);
				exit(1);
				}
			mincor = k1;
			}
		if (mincor < m2mapz_.maxz)
			mincor = m2mapz_.maxz;
		if (mincor > nz1) {
			if (z != z0)
				free(z);
			z = Malloc(mincor*sizeof(real));
			nz1 = mincor;
			}
		/* m3inpt_ ignores nwcore!  Must supply it via m2mapz_.maxz */
		m2mapz_.maxz = nwcore = nz1;
		ifuser = 0;
		m3dflt_(&I3);
		m1time_(&I1, &I0);
		m3inpt_(&objadd, z, &nwcore);
		m1time_(&IM1, &I0);
		if (m5log1_.ierr)
			return;
		mimode = 1;
		m3mps1_.nname = m3len_.nb;
		lrc = m5loc_.lpi + m3len_.m;

		z1 = z - 1;
		misolv_(&mimode, "Cold", &m3len_.m, &m3len_.n, &m3len_.nb,
			&m2mapa_.ne, &m2mapa_.nka, &m3mps1_.nname, &m5lobj_.iobj,
			&objadd, &z1[m2mapa_.la], (fint*)&z1[m2mapa_.lha],
			(fint*)&z1[m2mapa_.lka], &z1[m3loc_.lbl], &z1[m3loc_.lbu],
			(fint*)&z1[m3mps1_.lname1], (fint*)&z1[m3mps1_.lname2],
			(fint*)&z1[m3loc_.lhs], &z1[m5loc_.lxn], &z1[m5loc_.lpi],
			&z1[lrc], &inform, &ns, &z1[1], &nwcore, &ifuser, m1user,
			&majitn, &minitn, (ftnlen)4);

		m1time_(&I0, &I2);
		}
	}
#endif

 static fint
mkey(char *phrase, ftnlen len)
{
	char key[16];
	static fint lprint, lsumm = 6;
	fint inform = 0;
	m3key_(phrase, key, &lprint, &lsumm, &inform, len, (ftnlen)sizeof(key));
	return (int)inform;
	}

 static fint
nkey(fint *np, char *fname, ftnlen L)
{
	fint n = *np;
	char *s;

	if (n < 0 || n >= N_FNAMES) {
		printf("Unit number > %d\n", N_FNAMES-1);
		return 1;
		}
	fnames[n] = s = (char*)M1alloc(L+1);
	memcpy(s, fname, L);
	s[L] = 0;
	return 0;
	}

 static char *
set_outlev(Option_Info *oi, keyword *kw, char *s)
{
	char *rv = I_val(oi,kw,s);
	if (!*(int *)kw->info)
		oi->option_echo = 0;
	return rv;
	}

 static keyword keywds[] = {
	KW("ftimes",  I_val, &functimes, "report function eval. times"),
	KW("maxfwd",  IA_val, voffset_of(ASL,p.maxfwd_), "max vars in fwd AD of common exprs (default 5)"),
	KW("meminc",  D_val, &meminc, "extra megabytes of memory to give MINOS initially"),
	KW("objno",   L_val, &m8func_.nprob, "objective number: 0 = none, 1 = first (default)"),
	KW("outlev",  set_outlev, &outlev, "output level; 1 = default"),
	KW("timing",  I_val, &timing,  "report I/O and solution times: 1 = stdout, 2 = stderr, 3 = both"),
	KW("version", Ver_val, 0, "report version"),
	KW("wantsol", WS_val, 0, WS_desc_ASL+5)
	};

 static keyword options[] = {
	KW("f", IK1_val, &wantfuncs, "list available user-defined functions"),
	KW("t", IK1_val, &functimes, "time function evaluations")
	};

 static Option_Info Oinfo = {
	"minos", "MINOS 5.51", "minos_options", keywds, nkeywds, 1, xxxvers,
	usage_msg, mkey, nkey, options, sizeof(options)/sizeof(keyword), 20090923
	};

 static SufDecl
suftab[] = {
	{ "sstatus", 0, ASL_Sufkind_var, 1 },
	{ "sstatus", 0, ASL_Sufkind_con, 1 }
	};

 static int
envopt(char **argv)
{
	int rv;

	rv = getopts(argv, &Oinfo);
	if (fnames[7])
		m1file_.iprint = 7;
	if (fnames[0] && !freopen(fnames[0], "w", Stderr))
		printf("Can't redirect Stderr to %s\n", fnames[0]);
	fflush(stdout);
	return rv;
	}

#undef asl

 void
envopt_(void)
{
	Oinfo.n_keywds = 0;	/* stand-alone usage */
	envopt(xargv);
	}


 static fint
objmunge(fint M, fint *mp, fint N, fint NZ, fint *nzp, fint *hs,
		real *lb, real *ub, real *A, fint *ha, fint *ja,
		real *objadj)
{
	fint i, j, k, na0, ne, nlclim, nlvlim, nz, rv;
	cgrad *cg, **cgx, *ncg;
	ograd *og, **ogp;
	char *h, *he;
	real aijtol, ninf, pinf, t, t1;
	fint *ha1, *ha2, *hae, *ja1, *kadj, *kadj1, si;
	real *a1, *a2, *lbe;
	ASL *asl = cur_ASL;

#define objno m8func_.nprob
	rv = ne = na0 = 0;
	i = N*sizeof(fint);
	kadj = (fint *)Malloc(i);
	memset((char *)kadj, 0, i);
	aijtol = m3mps3_.aijtol;
	nz = NZ;
	nlvlim = nlvc;
	nlclim = nlc;

	/* omit tiny components of A */

	if (aijtol > 0) {
		a1 = a2 = A;
		ha1 = ha2 = ha;
		ja1 = ja;
		for(i = 0; i < N; i++) {
			j = *++ja1;
			hae = ha + j - 1;
			while(ha1 < hae) {
				t = *a1++;
				if ((si = *ha1++) > nlclim
				 && i >= nlvlim
				 && (t < 0 ? -t : t) < aijtol) {
					na0++;
					continue;
					}
				*a2++ = t;
				*ha2++ = si;
				}
			*ja1 = j - na0;
			}
		nz -= na0;
		}

	/* find objective, count gradient components */

	*objadj = 0;
	if (objno >= 0 && objno < n_obj) {
		obj_no = objno;
		if (!objsen)
			objsen = objtype[objno] ? -1 : 1;
		if (nl_obj(objno))
			rv = nlvo;
		else
			*objadj = objconst(objno);
		ogp = &Ograd[objno];
		for(og = *ogp; og; og = og->next)
			if (og->varno < rv)
				goto keep;
			else {
				if ((t = og->coef) < 0)
					t = -t;
				if (t >= aijtol) {
				 keep:
					kadj[og->varno] = 1;
					ne++;
					}
				else
					na0++;
				}
		}
	h = havex0;
	he = h + N;

	pinf = m1eps_.plinfy;
	ninf = -pinf;

	while(h < he) {
		if ((t = *ub) > pinf)
			t = *ub = pinf;
		if ((t1 = *lb) < ninf)
			t1 = *lb = ninf;
		lb++;
		ub++;
		if (*h++ || t1 <= ninf && t >= pinf)
			si = -1;
		else
			si = t <= 0.;
		*hs++ = si;
		}

	/* adjust for MINOS's surprising notion of bounds */

	lbe = lb + M;
	while(lb < lbe) {
		if ((t = -*lb) > pinf)
			t = pinf;
		if ((t1 = -*ub) < ninf)
			t1 = ninf;
		*lb++ = t1;
		*ub++ = t;
		*hs++ = t1 < t ? -1 : t <= 0.;
		}

	/* insert objective */

	si = M;
	if (ne) {
		m5lobj_.iobj = *mp = ++si;
		*lb = ninf;
		*ub = pinf;
		*hs = 0;
		a1 = A + nz;
		ha1 = ha + nz;
		k = nz;
		nz += ne;
		a2 = A + nz;
		ha2 = ha + nz;
		ja1 = ja + N;
		kadj1 = kadj + N;
		for(;;) {
			*ja1-- += ne;
			if (*--kadj1) {
				*kadj1 = a2 - A;
				*--a2 = 0;
				*--ha2 = si;
				if (!--ne)
					break;
				}
			for(j = *ja1; k >= j; k--) {
				*--a2 = *--a1;
				*--ha2 = *--ha1;
				}
			}
		while(og = *ogp)
			if ((i = kadj[og->varno]) && og->varno >= rv) {
				A[i - 1] = og->coef;
				*ogp = og->next;
				}
			else
				ogp = &og->next;
		}
	else
		*mp = M;
	*nzp = nz;
	free((char *)kadj);

	/* adjust for computing Jacobian */

	if (nlvlim) {
		if (!Cgrad) {
			Cgrad = (cgrad **)M1alloc(nlclim*sizeof(cgrad *));
			memset((char *)Cgrad, 0, nlclim*sizeof(cgrad *));
			}
		hae = ha + (ja[nlvlim] - 1);
		for(ha1 = ha, k = 0; ha1 < hae;)
			if (*ha1++ <= nlclim)
				k++;
		ncg = (cgrad *)M1alloc(k*sizeof(cgrad));
		i = k = 0;
		ja1 = ja + 1;
		for(a1 = A, ha1 = ha; ha1 < hae; a1++)
			if ((j = *ha1++) <= nlclim) {
				cgx = Cgrad + j - 1;
				cg = ncg++;
				cg->next = *cgx;
				*cgx = cg;
				cg->goff = k++;
				j = ha1 - ha;
				while(j >= *ja1) {
					ja1++;
					i++;
					}
				cg->varno = i;
				cg->coef = *a1;
				*a1 = 0;
				}
		}
	m3mps3_.na0 = na0;
	c_vars = nlvlim;
	return o_vars = rv;
	}

 static void
time_out(real tot, FILE *f)
{
	if (timing) {
		fprintf(f," MINOS times:\n read: %10.2f\n solve: %9.2f",
			m1tim_.tsum[0], tsolve);
		if (m1tim_.tsum[1] < tsolve)
			fprintf(f, "\texcluding minos setup: %.2f",
				m1tim_.tsum[1]);
		fprintf(f, "\n write: %9.2f\n total: %9.2f\n",
			m1tim_.tsum[2], tot);
		}
	if (m1tim_.ltime > 1 && m1tim_.numt[3] | m1tim_.numt[4]) {
		if (timing)
			fprintf(f, "\n");
		if (m1tim_.numt[3])
			fprintf(f,
	" constraints: %9.2f sec. for %ld evaluations, %ld Jacobians\n",
				m1tim_.tsum[3], (long)m1tim_.numt[3],
				(long)m8func_.nfcon[1]);
		if (m1tim_.numt[4])
			fprintf(f,
	" objective: %11.2f sec. for %ld evaluations, %ld gradients\n",
				m1tim_.tsum[4], (long)m1tim_.numt[4],
				(long)m8func_.nfobj[1]);
		}
	}

 static void
negate(fint M, real *x)
{
	real *xe;
	for(xe = x + M; x < xe; x++)
		*x = -*x;
	}

 static int
vtrans[] = { 0, 3, 2, 0, 1, 0, 0 },
ctrans[] = { 0, 3, 2, 1, 0, 0, 0 };

 static void
hs1_adjust(fint N, int *s, fint *ss, int *trans)
{
	fint i;

	for(i = 0; i < N; i++)
		ss[i] = trans[s[i]];
	}

 static char *
hs_adjust(ASL *asl, fint N, fint M, fint m1, fint *vss, int *vs, SufDesc *vsd, SufDesc *csd, real *A, fint *ha, fint *ja)
{
	fint *hae, i1, i2, is, *jae, nerror;
	int nlin;
	real t, *x, *xe;

	if (!(vsd->kind & ASL_Sufkind_input)
	 || !(csd->kind & ASL_Sufkind_input))
		return "Cold";
	hs1_adjust(N, vs, vss, vtrans);
	hs1_adjust(M, vs+N, vss+N, ctrans);
	if (m1 > M)
		vss[M+N] = 3;
	/* Why can't minos compute the initial slacks? */
	/* Then we could eliminate the following mess... */
	if (nlc) {
		nerror = 0;
		conval(X0, x = X0+N, &nerror);
		if (nerror)
			memset(X0+N, 0, nlc*sizeof(real));
		else
			for(xe = x + nlc; x < xe; x++)
				*x = -*x;
		}
	if ((nlin = n_con - nlc) > 0) {
		x = X0 + N + nlc;
		i1 = nlc + 1;
		i2 = n_con;
		memset(x, 0, nlin*sizeof(real));
		hae = ha;
		xe = X0;
		for(jae = ja + N; ja < jae; ja++) {
			hae += ja[1] - ja[0];
			for(t = *xe++; ha < hae; A++)
				if ((is = *ha++) >= i1 && is <= i2)
					x[is - i1] -= t**A;
			}
		}
	return "Warm";
	}

 static void
send_status(fint N, fint *ss, int *s, int *trans, real *L, real *U)
{
	fint i;

	for(i = 0; i < N; i++)
		if (((s[i] = trans[ss[i]]) + 1 & ~1) == 4 && L[i] == U[i])
			s[i] = 5;
	}

#ifndef Student_Edition
 static int
file_assignment(char *s)
{
	int c;

	if ((c = *s++) < '0' || c > '9')
		return 0;
	while((c = *s++) >= '0' && c <= '9');
	return c == '=';
	}
#endif

static Jmp_buf Jb;

 void
MAIN__(void)
{
	char names[40], *stub;
	real *A, *lb, *pi, *rc, *ub, *z, z1;
	static fint I1 = 1, I2 = 2, I3 = 3, I4 = 4;
	FILE *nl;
	fint INFORM, M, MXROW, MXCOL, N, N1, NB, NO, NS, NZ;
	fint i, k, k1, m1, mincor, ninf, nint, nresize;
	fint nwcore, nz;
	fint *ha, *hs, *ja, name1[1], name2[1];
	char buf[32], msg[400], *msg1, *start;
	real obj, objadj, sinf, t, tout, tstart;
	ASL *asl;
	int *varstat;
	SufDesc *csd, *vsd;
	static int	sctrans[] = { 4, 3, 2, 1 },
			svtrans[] = { 3, 4, 2, 1 };
	typedef struct { char *msg; int code; int wantobj; } Sol_info;
	static Sol_info solinfo[] = {
	 {/* 0 */ "optimal solution found", 000, 1},
	 {/* 1 */ "infeasible problem (or bad starting guess)", 200, 0},
	 {/* 2 */ "unbounded (or badly scaled) problem", 300, 0},
	 {/* 3 */ "too many iterations", 400, 1},
	 {/* 4 */ "the objective has not changed for the last %ld iterations", 500, 1},
	 {/* 5 */ "the superbasics limit (%ld) is too small", 520, 1},
	 {/* 6 */ "error evaluating nonlinear expressions", 521, 1},
	 {/* 7 */ "incorrect gradients from funobj", 530, 0},
	 {/* 8 */ "incorrect gradients from funcon", 531, 0},
	 {/* 9 */ "the current point cannot be improved", 501, 1},
	 {/* 10*/ "numerical error: the general constraints\ncannot be satisfied accurately", 201, 1},
	 {/* 11*/ "cannot find superbasic to replace basic variable", 532, 1},
	 {/* 12*/ "basis factorization requested twice in a row", 533, 1},
	 {/* 13*/ "optimal solution found?  Optimality\ntests satisfied, but reduced gradient is large", 100, 1},
	 {/* 14*/ "not enough storage for the basis factors.\nTry rerunning with workspace_(total)=%ld in $minos_options ", 522, 0},
	 {/* 15*/ "error in basis package", 534, 0},
	 {/* 16*/ "singular basis after several factorization attempts", 510, 1},
	 {/* 17*/ "input basis had wrong dimensions", 535, 0},
	 {/* 18*/ "unexpected return code (%ld)", 536, 0},
	 {/* 19*/ "solution aborted", 540, 1},
	 {/* 20*/ "too many major iterations", 401, 1},
	 {/* 21*/ "infeasible problem", 200, 0},
	 {/* 22*/ "cannot allocate enough memory to solve the problem", 550, 0},
	 {/* 23*/ "more than 2^31 double-precision words of memory needed", 551, 0}
	 };

	tstart = xectim_();
	asl = ASL_alloc(ASL_read_fg);
	asl->i.congrd_mode = 2;	/* sparse Jacobians */
	stub = getstub(&xargv, &Oinfo);
	if (wantfuncs) {
		show_funcs();
		exit(0);
		}
#ifdef NO_MINOS1
	if (!stub)
		usage_ASL(&Oinfo, 1);
#else
	if (!stub || file_assignment(stub)) {
		minos1_(&z1, &I1);
		return;
		}
#endif

	nl = jacdim(stub, &M, &N, &NO, &NZ, &MXROW, &MXCOL, (fint)0);

	if (N <= 0) {
		fprintf(Stderr, "%s has no variables\n", filename);
		exit(4);
		}

	suf_declare(suftab, sizeof(suftab)/sizeof(SufDecl));
	N1 = N + 1;
	m1 = M + 1;
	NB = N + m1;
	nz = NZ + N;
	Fortran = 1;

	LUv = lb = (real *)M1alloc((4*NB + m1 + nz)*sizeof(real)
					+ (N1 + NB + nz)*sizeof(fint) + N);
	LUrhs = lb + N;
	Uvx = ub = LUrhs + m1;
	Urhsx = ub + N;
	X0 = Urhsx + m1;
	A = A_vals = X0 + NB;
	memset((char *)(X0 + N), 0, m1*sizeof(real));
	rc = A + nz;
	pi = pi0 = rc + NB;
	ja = (fint *)(pi + m1);
	ha = ja + N1;
	hs = ha + nz;
	varstat = (int*)M1alloc(NB*sizeof(int));
	vsd = suf_iput("sstatus", ASL_Sufkind_var, varstat);
	csd = suf_iput("sstatus", ASL_Sufkind_con, varstat + N);
	havex0 = (char *)(hs + NB);
	A_rownos = (int *)ha;
	if (sizeof(int) == sizeof(fint))
		A_colstarts = (int *)ja;
	else
		A_colstarts = (int *)M1alloc(N1*sizeof(int));
	Fortran = 1;
	m1init_();
	m1env_.ampl = 1; /* prevent "Solution not printed" */
	asl->i.nlvog = nlvo;
	fg_read(nl,0);
	if (sizeof(int) != sizeof(fint)) {
		fint *ja1 = ja;
		fint *ja1e = ja + N1;
		int *ja2 = A_colstarts;
		while(ja1 < ja1e)
			*ja1++ = *ja2++;
		}

	m3dflt_(&I1);
	objsen = 0;
	m7len_.nnobj = objmunge(M, &m1, N, NZ, &nz, hs, lb, ub, A, ha, ja, &objadj);
	m8al1_.nden = 2;	/* sparse Jacobian! */
	m8len_.nncon = n_conjac[1] = nlc;
	m8len_.nnjac = nlvc;
	m8veri_.lverif[0] = m5freq_.msoln = -11111;
	m5lp1_.itnlim = 99999999;
	objno = n_obj > 0;
	if (functimes)
		timing = 1;
	m2len_.mrows = m1;
	m2len_.mcols = N;
	if (envopt(xargv))
		exit(2);
	m8len_.njac = nzjac;
	m3dflt_(&I2);
	m1tim_.ltime = functimes + 1;
	if (objno > n_obj) {
		printf("objno = %ld must be <= %d\n", (long)objno, n_obj);
		exit(2);
		}
	--objno;
	if ((obj_no = (int)m8func_.nprob)
	 && (obj_no < -1 || obj_no >= n_obj)) {
		printf("Bad problem number %d: not in [0,%d]\n",obj_no,n_obj);
		exit(2);
		}

	if (nint = nlogv + niv + nlvbi + nlvci + nlvoi) {
		printf("ignoring integrality of %ld variables\n", (long)nint);
		need_nl = 0;
		}

	if (m5freq_.msoln == -11111)
		m5freq_.msoln = 0;
	if (outlev > 1) {
		if (outlev == 2)
			m1file_.iprint = m1file_.isumm = 6;
		else {
			m1file_.iprint = 6;
			if (outlev > 3)
				m5freq_.msoln = 2;
			}
		}
	if (m8veri_.lverif[0] == -11111
	 && (m1file_.iprint > 0 ||  m1file_.isumm > 0))
		m8veri_.lverif[0] = 0;
	if ((m1file_.isumm == 6 || m1file_.iprint == 6)
	 && need_nl) {
		printf("\n");
		need_nl = 0;
		}
	mifile_(&I1);
	start = hs_adjust(asl, N, M, m1, hs, varstat, vsd, csd, A, ha, ja);
	NB = N + m1;
	i = m8len_.njac;	/* used not to be necessary */
	m2core_(&I1, &mincor);
	memused0 = memused = mincor * (real)sizeof(real);
	if (meminc > 0.) {
		k = (1 << 17) * meminc;
		k1 = mincor + k;
		if (mincor >= k1) {
			fprintf(Stderr,
				"meminc = %g megabytes is too big a memory increment.\n",
				meminc);
			i = 23;
			t = xectim_();
			goto no_mem;
			}
		mincor = k1;
		memused = mincor * (real)sizeof(real);
		}
	memused1 = memused;
	m8len_.njac = i;	/* scrogged by m2core */
	if (mincor < m2mapz_.maxz)
		mincor = m2mapz_.maxz;
	else
		m2mapz_.maxz = mincor;
	nwcore = mincor;
	m2core_(&I4, &i);
	if (i > mincor) {
		fprintf(Stderr,
			"MINOS first wanted %ld words, now it wants %ld\n",
			(long)mincor, (long)i);
		exit(2);
		}

	/* open files mentioned on command line or in $minos_options */
	mifile_(&I2);

	/* show parameters if print_level >= 1 and outlev > 1 */
	if (outlev > 1 && m5log4_.prnt1)
		m3dflt_(&I3);

	/* The following nonsense used not to be necessary. */
	if (outlev == 2)
		m1file_.iprint = m5log4_.prnt0 = m5log4_.prnt1 = 0;

	t = xectim_();
	m1tim_.tsum[0] = t - tstart;
	err_jmp1 = &Jb;
	if (setjmp(Jb.jb))
 longjumped:
		i = 19;
	else {
		if (setjmp(fpe_jmpbuf)) {
			report_where(asl);
			printf("\nFloating point error.\n");
			fflush(stdout);
			need_nl = 0;
			goto longjumped;
			}
		signal(SIGFPE, fpecatch);
		obj = Infinity;
		if (objsen < 0)
			negate(M,pi);
		if (m1 > M)
			pi[M] = 0;
		for(nresize = 0;; ++nresize) {
			z = (real *)malloc(mincor*sizeof(real));
			if (!z) {
				i = 22;
				goto no_mem;
				}
			memset(z, 0, mincor*sizeof(real));
			INFORM = NS = 0;
			memset(names, ' ', 40);
			minoss_(start, &m1, &N, &NB, &nz, &I1,
				&m8len_.nncon, &m7len_.nnobj, &m8len_.nnjac,
				&m5lobj_.iobj, &objadj, names,
				A, ha, ja, lb, ub, name1, name2,
				hs, X0, pi, rc,
				&INFORM, &mincor, &NS, &ninf, &sinf, &obj,
				z, &nwcore, (ftnlen)4, (ftnlen)8);
			if (INFORM != 42)
				break;
			free(z);
			if (nwcore == 0x7fffffff) {
				i = 23;
				goto no_mem;
				}
			mincor += mincor >> 1;
			if (mincor < 0)
				mincor = 0x7fffffff;
			memused = mincor * (real)sizeof(real);
			m2mapz_.maxz = nwcore = mincor;	
			}
		M1record(z);	/* to be freed by ASL_free */
		if ((i = INFORM) < 0 || i > 13)
			i = INFORM >= 20 && INFORM <= 22 ? INFORM + (14-20)
				: INFORM == 30 ? 17 : 18;
		}
 no_mem:
	tout = xectim_();
	tsolve = tout - t;
	if (need_nl && timing & 1) {
		printf("\n");
		need_nl = 0;
		}
	msg1 = msg + Sprintf(msg, "%s: ", Oinfo.bsname);
	switch(i) {
		case 1:
			if (!nlc)
				i = 21;
			goto have_i;
		case 3:
			if (m5lp1_.itn < m5lp1_.itnlim)
				i = 20;
			goto have_i;
		case 4:
			nint = 2*NB;
			if (nint < 200)
				nint = 200;
			else if (nint > 1000)
				nint = 1000;
			msg1 += Sprintf(msg1, solinfo[4].msg, (long)nint);
			break;
		case 5:
			msg1 += Sprintf(msg1, solinfo[5].msg,
					(long)m5len_.maxs);
			break;
		case 14:	/* imitate computation of minlen in m2bfac */
			INFORM = m2lu1_.nbelem*5/4;
			if (INFORM <= m2lu1_.lena)
				INFORM = m2lu1_.lena;
			else if (INFORM < m2lu4_.luparm[12])
				INFORM = m2lu4_.luparm[12];
			INFORM = m2mapz_.maxz + 3*(INFORM - m2lu1_.lena);
			/* no break */
		default:
 have_i:
			msg1 += Sprintf(msg1, solinfo[i].msg, (long)INFORM);
		}
	msg1 += Sprintf(msg1, ".\n%ld iterations", (long)m5lp1_.itn);
	if (solinfo[i].wantobj && obj != Infinity) {
		g_fmtop(buf, m5lobj_.iobj ? obj : objadj);
		msg1 += Sprintf(msg1, ", objective %s", buf);
		}
	solve_result_num = solinfo[i].code;
	if (m8func_.nfcon[0] + m8func_.nfobj[0] > 0 && !functimes) {
		msg1 += Sprintf(msg1, "\nNonlin evals: ");
		if (m8func_.nfobj[0] > 0)
			msg1 += Sprintf(msg1, "obj = %ld, grad = %ld%s",
				(long)m8func_.nfobj[0],
				(long)m8func_.nfobj[1],
				m8func_.nfcon[0] > 0 ? ", " : "");
		if (m8func_.nfcon[0] > 0)
			msg1 += Sprintf(msg1, "constrs = %ld, Jac = %ld",
				(long)m8func_.nfcon[0],
				(long)m8func_.nfcon[1]);
		msg1 += Sprintf(msg1, ".");
		}
	if (memused > memused1)
		Sprintf(msg1, "\nAdding  meminc=%.3g  to $minos_options might save time.",
			1e-6*(memused-memused0));
	send_status(N, hs, varstat, svtrans, LUv, Uvx);
	send_status(M, hs+N, varstat+N, sctrans, LUrhs, Urhsx);
	write_sol(msg, X0, pi, &Oinfo);
	if (timing | functimes) {
		t = xectim_();
		m1tim_.tsum[2] = t - tout;
		t -= tstart;
		if (!timing || timing & 1)
			time_out(t, stdout);
		if (timing & 2) {
			fflush(stdout);
			time_out(t, Stderr);
			}
		}
	ASL_free(&asl);	/* for Purify */
	}

 void
gfname_(fint *i0, char *fname, ftnlen fname_len)
{
	int i = *i0, L;
	char *s;

	if (i >= N_FNAMES  || i < 0) {
		fprintf(Stderr, "gfname called with i = %d\n", i);
		exit(1);
		}
	if (s = fnames[i])
		L = Sprintf(fname, "%s", s);
	else
		L = Sprintf(fname, "fort.%d", i);
	while (L < fname_len)
		fname[L++] = ' ';
	}

#undef scream
 static void
scream(char *fmt, int j, int k)
{
	fprintf(Stderr, fmt, j, k);
	exit(1);
	}
 void
funcon_(fint *MODE, fint *M, fint *N, fint *NJAC,
		real *X, real *F, real *G,
		fint *NSTATE, fint *NPROB, real *Z, fint *NWCORE)
{
	ASL *asl = cur_ASL;

	if (*NSTATE) {
		if (*NSTATE != 1)
			return;
		if (*N != c_vars) {
			scream("funcon expected N = %d but got %d\n",
				c_vars, (int)*N);
			 /* suppress warning about unused vars: */
			Not_Used(NJAC);
			Not_Used(Z);
			Not_Used(NWCORE);
			Not_Used(NPROB);
			}
		}
	want_deriv = (int)*MODE & 2;
	xknown(X);
	conval(X, F, 0);
	if (want_deriv)
		jacval(X, G, 0);
	xunknown();
	}

 void
funobj_(fint *MODE, fint *N,
		real *X, real *F, real *G,
		fint *NSTATE, fint *NPROB, real *Z, fint *NWCORE)
{
	int i;
	ASL *asl = cur_ASL;

	if (*NSTATE) {
		if (*NSTATE != 1)
			return;
		if (*N != o_vars) {
			scream("funobj expected N = %d but got %d\n",
				o_vars, (int)*N);
			printf("", Z, NWCORE); /* use unused vars */
			}
		}
	want_deriv = (int)*MODE & 2;
	i = (int)*NPROB;
	if (i < 0 || i >= n_obj) {
		*F = 0.;
		if (want_deriv)
			memset(G, 0, *N*sizeof(real));
		}
	else {
		xknown(X);
		*F = objval(i,X,0);
		if (want_deriv)
			objgrd(i,X,G,0);
		xunknown();
		}
	}

 void
matmod_(fint *ncycle, fint *nprob, fint *finish,
	fint *m, fint *n, fint *nb, fint *ne, fint *nka,
	fint *ns, fint *nscl, fint *nname, real *a, fint *ha,
	fint *ka, real *bl, real *bu, real *ascale,
	fint *hs, fint *name1, fint *name2, real *x,
	real *pi, real *rc, real *z, fint *nwcore)
{
	*finish = 1;
	}

#ifdef __cplusplus
}
#endif
