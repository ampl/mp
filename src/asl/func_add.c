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

#include "asl.h"

#ifdef __cplusplus
extern "C" {
#endif
extern ASLhead ASLhead_ASL;

const char *i_option_ASL;
unsigned long randseed_ASL;
static int n_added;

 void
addrandinit_ASL(AmplExports *ae, RandSeedSetter rss, void *v)
{
	char *s, *se;
	unsigned long x;

	/* Something more elaborate will be needed if we add a "reset" facility to ASL. */

	if (!randseed_ASL) {
		randseed_ASL = 1;
		if ((s = getenv("randseed"))
		 && (x = (unsigned long)strtol(s,&se,10))
		 && !se)
			randseed_ASL = x;
		}
	(*rss)(v, randseed_ASL);
	}

 func_info *
func_lookup(ASL *asl, register const char *s, int add)
{
	register unsigned x = 0;
	func_info *fi, **finext;
	Const char *s0 = s;

	while(*s)
		x = 31*x + *s++;
	finext = &fhash[x % NFHASH];
	for(fi = *finext; fi; fi = fi->next)
		if (!strcmp(s0, fi->name)) {
			if (add) {
				fprintf(Stderr,
				"addfunc: duplicate function %s\n", s0);
				fi = 0;
				}
			return fi;
			}
	if (add) {
		fi = (func_info *)mem_ASL(asl, sizeof(func_info));
		fi->next = *finext;
		*finext = fi;
		fi->name = s0;
		}
	return fi;
	}

 void
addfunc_ASL(const char *fname, ufunc *f, int ftype, int nargs, void *funcinfo, AmplExports *ae)
{
	register func_info *fi;
	ASL *asl = (ASL*)ae->asl;
	if (ftype && ftype != 1) {
#ifndef COMPLAIN_AT_BAD_FTYPE
		if (ftype < 0 || ftype > 6)
#endif
		{
		fprintf(Stderr, "function %s: ftype = %d; expected 0 or 1\n",
			fname, ftype);
		exit(1);
		}
#ifndef COMPLAIN_AT_BAD_FTYPE
		return;
#endif
		}
	if ((fi = func_lookup(asl, fname, 1))) {
		n_added++;
		fi->funcp = f;
		fi->ftype = ftype;
		fi->nargs = nargs;
		fi->funcinfo = funcinfo;
		if (!funcsfirst)
			funcsfirst = fi;
		else
			funcslast->fnext = fi;
		funcslast = fi;
		fi->fnext = 0;
		}
	}

enum { NEFB = 5, NEFB0 = 2 };

 static Exitcall a_e_info[NEFB0];
 static Exitcall *a_e_next = a_e_info;
 static Exitcall *a_e_last = a_e_info + NEFB0;
 static Exitcall *a_e_prev;

 typedef struct
ExitCallInfo { Exitcall *cur, **curp, *last, **lastp; } ExitCallInfo;

 static void
AtReset1(AmplExports *ae, Exitfunc *ef, void *v, ExitCallInfo *eci)
{
	Exitcall *ec;
	ASL *asl = (ASL*)ae->asl;
	if (eci) {
		eci->cur = asl->i.arprev;
		eci->curp = &asl->i.arprev;
		eci->last = asl->i.arlast;
		eci->lastp = &asl->i.arlast;
		}
	if (asl->i.arnext >= asl->i.arlast) {
		asl->i.arnext = (Exitcall*)M1alloc(NEFB*sizeof(Exitcall));
		asl->i.arlast = asl->i.arnext + NEFB;
		}
	asl->i.arnext->prev = asl->i.arprev;
	asl->i.arprev = ec = asl->i.arnext++;
	ec->ef = ef;
	ec->v = v;
	}

 static void
AtReset(AmplExports *ae, Exitfunc *ef, void *v)
{ AtReset1(ae, ef, v, 0); }


 void
at_end_ASL(Exitcall *ec)
{
	while(ec) {
		(*ec->ef)(ec->v);
		ec = ec->prev;
		}
	}

 void
at_exit_ASL(VOID)
{
	Exitcall *ec;
	ASLhead *h, *h0;

	h0 = &ASLhead_ASL;
	h = ASLhead_ASL.next;
	h0->next = h0->prev = h0;
	for(; h != h0; h = h->next)
		if ((ec = ((ASL*)h)->i.arprev))
			at_end_ASL(ec);
	if ((ec = a_e_prev)) {
		a_e_prev = 0;
		at_end_ASL(ec);
		}
	}

 static void
AtExit1(AmplExports *ae, Exitfunc *ef, void *v, ExitCallInfo *eci)
{
	Exitcall *ec;
	Not_Used(ae);
#ifndef NO_ONEXIT
	if (!a_e_prev)
		atexit(at_exit_ASL); /* in case mainexit() is bypassed */
#endif
	if (eci) {
		eci->cur = a_e_prev;
		eci->curp = &a_e_prev;
		eci->last = a_e_last;
		eci->lastp = &a_e_last;
		}
	if (a_e_next >= a_e_last) {
		a_e_next = (Exitcall*)mymalloc(NEFB*sizeof(Exitcall));
		a_e_last = a_e_next + NEFB;
		}
	a_e_next->prev = a_e_prev;
	a_e_prev = ec = a_e_next++;
	ec->ef = ef;
	ec->v = v;
	}

 static void
AtExit(AmplExports *ae, Exitfunc *ef, void *v)
{ AtExit1(ae, ef, v, 0); }

 static Char *
Tempmem(TMInfo *T, size_t L)
{
	TMInfo *T1 = (TMInfo *)mymalloc(L + sizeof(TMInfo));
	T1->u.prev = T->u.prev;
	T->u.prev = T1;
	return (Char*)(T1+1);
	}

 static void
No_table_handler(
	int (*Dbread)(AmplExports*, TableInfo*),
	int (*Dbwrite)(AmplExports*, TableInfo*),
	char *hname,
	int flags,
	void *vinfo)
{}

 static cryptblock*
No_crypto(char *key, size_t scrbytes)
{ return 0; }

typedef void Funcadd ANSI((AmplExports*));

#ifndef CLOSE_AT_RESET
static Funcadd *Fa0[4], **Fa = Fa0;
static int nFa = 0, nFamax = 4;
#endif

#ifdef SYMANTEC
#define No_popen_or_pclose
typedef char *(*Tempnamtype)(const char*, const char*);
#define Tempnam_cast (Tempnamtype)
#endif

#ifdef WATCOM
#define tempnam _tempnam
#endif

#ifdef NO_tempnam

/* If the system does not provide a true tempnam function */
/* the AMPL/solver interface library will not do so either. */

 static char *
tempnam(const char *dir, const char *pfx)
{ return 0; }
#endif /* NO_tempnam */

#ifdef _WIN32
#define popen _popen
#define pclose _pclose
#else
#ifdef MSDOS
#undef  No_popen_or_pclose
#define No_popen_or_pclose
#endif
#endif

#ifdef No_popen_or_pclose
#undef popen
#define popen no_popen
#undef pclose
#define pclose no_pclose

 static int
no_pclose(FILE*f) { return 1; }

 static FILE*
no_popen(const char*cmd, const char*type) { return 0; }
#endif

 static AmplExports AE;

#ifdef clearerr
 static void
myclearerr(FILE *f)
{ clearerr(f); }
#undef clearerr
#define clearerr myclearerr
#endif /*clearerr*/

#ifdef feof
 static int
myfeof(FILE *f)
{ return feof(f); }
#undef feof
#define feof myfeof
#endif /*feof*/

#ifdef ferror
 static int
myferror(FILE *f)
{ return ferror(f); }
#undef ferror
#define ferror myferror
#endif /*ferror*/

#ifdef _fileno
#undef fileno
#define fileno _fileno
#endif

#ifdef fileno
 static int
myfileno(FILE *f)
{ return fileno(f); }
#undef fileno
#define fileno myfileno
#endif /* fileno */

#ifndef Tempnam_cast
#define Tempnam_cast /*nothing*/
#endif

#ifdef __linux__
#define USE_MKSTEMP
#endif
#ifdef USE_MKSTEMP
 /* Shut up warnings about tempnam and tmpnam. */
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

 static int
isdir(const char *s)
{
	struct stat sbuf;
	if (stat(s, &sbuf))
		return 0;
	return S_ISDIR(sbuf.st_mode);
	}

 static char *
my_tempnam(const char *dir, const char *pfx, char *s)
{
	const char *c;
	int i;
	size_t Ld, Lp;

	if ((c = getenv("TMPDIR")) && isdir(c))
		dir = c;
	else if (!dir || !isdir(dir))
		dir = "/tmp";
	if (!pfx)
		pfx = "";
	Ld = strlen(dir);
	Lp = strlen(pfx);
	if (!s)
		s = (char*)Malloc(Ld + Lp + 8);
	strcpy(s, dir);
	if (s[Ld-1] != '/')
		s[Ld++] = '/';
	strcpy(s+Ld, pfx);
	strcpy(s + Ld + Lp, "XXXXXX");
	if ((i = mkstemp(s)))
		close(i);
	else {
		free(s);
		s = 0;
		}
	return s;
	}

 static char *
Tempnam(const char *dir, const char *pfx)
{ return my_tempnam(dir,pfx,0); }

 static char *
Tmpnam(char *s)
{
	static char *s0;
	if (s)
		return my_tempnam(0,"Temp_",s);
	if (s0)
		free(s0);
	return s0 = my_tempnam(0,"Temp_",0);
	}
#undef tempnam
#define tempnam Tempnam
#undef tmpnam
#define tmpnam Tmpnam
#endif /* USE_MKSTEMP */

 void (*breakfunc_ASL) ANSI((int,void*)), *breakarg_ASL;

 void
func_add(ASL *asl)
{
	AmplExports *ae;

	if (need_funcadd) {
		if (!i_option_ASL
		 && !(i_option_ASL = getenv("ampl_funclibs")))
		      i_option_ASL = getenv("AMPLFUNC");
		if (!AE.PrintF) {
			AE.StdIn = stdin;
			AE.StdOut = stdout;
			AE.StdErr = Stderr;
			AE.ASLdate = ASLdate_ASL;
			AE.Addfunc = addfunc_ASL;
			AE.PrintF = printf;
			AE.FprintF = fprintf;
			AE.SprintF = sprintf;
			AE.SnprintF = snprintf;
			AE.VfprintF = vfprintf;
			AE.VsprintF = vsprintf;
			AE.VsnprintF = vsnprintf;
			AE.Strtod = strtod;
			AE.AtExit = AtExit;
			AE.AtReset = AtReset;
			AE.Tempmem = Tempmem;
			AE.Add_table_handler = No_table_handler;
			AE.Crypto = No_crypto;
			AE.Qsortv = qsortv;
			AE.Clearerr = clearerr;
			AE.Fclose = fclose;
			AE.Fdopen = fdopen;
			AE.Feof = feof;
			AE.Ferror = ferror;
			AE.Fflush = fflush;
			AE.Fgetc = fgetc;
			AE.Fgets = fgets;
			AE.Fileno = fileno;
			AE.Fopen = fopen;
			AE.Fputc = fputc;
			AE.Fputs = fputs;
			AE.Fread = fread;
			AE.Freopen = freopen;
			AE.Fscanf = fscanf;
			AE.Fseek = fseek;
			AE.Ftell = ftell;
			AE.Fwrite = fwrite;
			AE.Pclose = pclose;
			AE.Perror = perror;
			AE.Popen = popen;
			AE.Puts = puts;
			AE.Rewind = rewind;
			AE.Scanf = scanf;
			AE.Setbuf = setbuf;
			AE.Setvbuf = setvbuf;
			AE.Sscanf = sscanf;
			AE.Tempnam = Tempnam_cast tempnam;
			AE.Tmpfile = tmpfile;
			AE.Tmpnam = tmpnam;
			AE.Ungetc = ungetc;
			AE.Getenv = getenv_ASL;
			AE.Breakfunc = breakfunc_ASL;
			AE.Breakarg = breakarg_ASL;
			AE.Addrandinit = addrandinit_ASL;
			}
		if (AE.asl)
			memcpy(ae = (AmplExports*)M1alloc(sizeof(AmplExports)),
				&AE, sizeof(AmplExports));
		else
			ae = &AE;
		asl->i.ae = ae;
		ae->asl = (Char*)asl;
		auxinfo_ASL(ae);
#ifndef CLOSE_AT_RESET
		if (nFa > 0) {
			/* not the first nl_reader call */
			int i = 0;
			while(i < nFa)
				(*Fa[i++])(ae);
			}
		else
#endif
			funcadd(ae);
		need_funcadd = 0;
		}
	}

 void
show_funcs_ASL(ASL *asl)
{
	func_info *fi;
	int nargs;
	const char *atleast;

	func_add(asl);
	fprintf(Stderr, "Available nonstandard functions:%s\n",
		(fi = funcsfirst) ? "" : " none");
	for(; fi; fi = fi->fnext) {
		if ((nargs = fi->nargs) >= 0)
			atleast = "";
		else {
			nargs = -(1 + nargs);
			atleast = "at least ";
			}
		fprintf(Stderr, "\t%s(%s%d %sarg%s)\n", fi->name,
			atleast, nargs, fi->ftype ? "" : "real ",
			nargs == 1 ? "" : "s");
		}
	fflush(Stderr);
	}

 void
note_libuse_ASL(void)
{ ++n_added; }

 int
aflibname_ASL(AmplExports *ae, const char *fullname, const char *name, int nlen,
	Funcadd *fa, int save_fa, void (*dl_close)(void*), void *h)
{
	Exitcall *ec;
	ExitCallInfo eci;

	af_libnamesave_ASL(ae, fullname, name, nlen);
	n_added = 0;
	if (save_fa)
		AtExit1( ae, dl_close, h, &eci);
	else
		AtReset1(ae, dl_close, h, &eci);
	(*fa)(ae);
	if (!n_added) {
		for(ec = *eci.curp; ec != eci.cur; ec = ec->prev)
			(*ec->ef)(ec->v);
		*eci.curp = ec;
		*eci.lastp = eci.last;
		/* A small storage leak is possible if a new block of */
		/* Exitcalls was allocated, but since the present !n_added */
		/* case is unlikely, this leak should be of little concern. */
		}
#ifndef CLOSE_AT_RESET
	else if (save_fa) {
		if (++nFa >= nFamax) {
			Funcadd **Fa1;
			nFamax <<= 1;
			Fa1 = (Funcadd**)Malloc(nFamax * sizeof(Funcadd*));
			memcpy(Fa1, Fa, nFa*sizeof(Funcadd*));
			if (Fa != Fa0)
				free(Fa);
			Fa = Fa1;
			}
		Fa[nFa-1] = fa;
		}
#endif /* !CLOSE_AT_RESET */
	return n_added;
	}
#ifdef __cplusplus
}
#endif
/* Last relevant change to asl.h: 19991013. */
