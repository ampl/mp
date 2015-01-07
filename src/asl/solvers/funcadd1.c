/****************************************************************
Copyright (C) 1998, 1999, 2000 Lucent Technologies
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

#ifdef NO_FUNCADD
#include "funcadd.h"

const char *ix_details_ASL[] = {0};

 void
funcadd(AmplExports *ae)
{ ae = ae; /* shut up non-use warning */ }

#else

#ifdef _WIN32
#undef WIN32
#define WIN32
#endif

#ifdef WIN32
#include "windows.h"
#undef void
#endif

#define _POSIX_SOURCE	/* for HP-UX */

#include "stdlib.h"	/* for free() */
#include "string.h"
#include "funcadd.h"
#include "arith.h"	/* for X64_bit_pointers */
#include <sys/types.h>
#include <sys/stat.h>
#ifndef S_IFREG /*{*/
#ifdef __S_IFREG
#define S_IFREG __S_IFREG
#define S_IFDIR __S_IFDIR
#elif defined(_S_IFREG)
#define S_IFREG _S_IFREG
#define S_IFDIR _S_IFDIR
#endif
#endif /*}*/
#ifdef X64_bit_pointers
static char Bits[] = "64", BitsAlt[] = "32";
#else
static char Bits[] = "32", BitsAlt[] = "64";
#endif

#ifdef Old_APPLE	/* formerly __APPLE__, for earlier versions of Mac OS X */
#define FUNCADD "_funcadd_ASL"
#endif
#ifndef FUNCADD
#define FUNCADD "funcadd_ASL"
#endif

#ifdef __cplusplus
extern "C" {
extern int libload_ASL(AmplExports *ae, const char *s, int ns, int warn);
#endif

typedef void Funcadd ANSI((AmplExports*));

extern void *mymalloc_ASL ANSI((size_t));
#undef mymalloc
#define mymalloc(x) mymalloc_ASL((size_t)(x))

const char *ix_details_ASL[] = {
	"? {show -i options}",
	"- {do not import functions: do not access amplfunc.dll}",
	"dir {look for amplfunc.dll in directory dir}",
	"file {import functions from file rather than amplfunc.dll}",
	"",
	"When the x of -ix is suitably quoted, multiple files may appear on",
	"separate lines or may appear on the same line if each is enclosed",
	"by single or double quotes.",
	"",
	"If no -i option appears but $ampl_funclibs is set, assume",
	"-i $ampl_funclibs.  Otherwise, if $AMPLFUNC is set, assume",
	"-i $AMPLFUNC.  Otherwise look for amplfunc.dll in the",
	"directory that is current when execution begins.",
	"",
	"-ix and -i x are treated alike.",
	0 };
#define afdll afdll_ASL
extern int aflibname_ASL ANSI((AmplExports*, const char*, const char*, int, Funcadd*, int, void(*)(void*), void*));
extern const char *i_option_ASL;

#ifdef __cplusplus
	}
#endif

static int first = 1;

 static int
file_kind(const char *name) /* 1 == regular file, 2 ==> directory; else 0 */
{
	struct stat sb;

	if (stat(name,&sb))
		return 0;
	if (sb.st_mode & S_IFDIR)
		return 2;
	if (sb.st_mode & S_IFREG)
		return 1;
	return 0;
	}

#ifdef WIN32

#define SLASH '\\'
char afdll[] = "\\amplfunc.dll";
typedef HINSTANCE shl_t;
#define dlopen(x,y) LoadLibrary(x)
#define find_dlsym(a,b,c) (a = (Funcadd*)GetProcAddress(b,c))
#define dlclose(x) FreeLibrary((HMODULE)x)
#define NO_DLERROR

 static int
Abspath(const char *s)
{
	int c = *s;
	if ((c >= 'a' && c <= 'z' || c >= 'A' && c <= 'Z')
	 && s[1] == ':'
	 && (c = s[2]) == '\\' || c == '/')
		return 1;
	return 0;
	}

#else /* !WIN32 */

#define SLASH '/'

char afdll[] = "/amplfunc.dll";

#define Abspath(s) (*(s) == '/')

#include "unistd.h"	/* for getcwd */
#define GetCurrentDirectory(a,b) getcwd(b,(int)(a))

#ifdef __hpux
#include "dl.h"
#define dlopen(x,y) shl_load(x, BIND_IMMEDIATE, 0)
#define find_dlsym(a,b,c) !shl_findsym(&b, c, TYPE_PROCEDURE, &a)
#define dlclose(x) shl_unload((shl_t)x)
#define NO_DLERROR
#else
#ifdef Old_APPLE
#include <mach-o/dyld.h>
typedef struct {
	NSObjectFileImage ofi;
	NSModule m;
	char *name;
	} NS_pair;
typedef NS_pair *shl_t;

 static void*
find_sym_addr(NS_pair *p, const char *name)
{
	NSSymbol nss;

	if (nss = NSLookupSymbolInModule(p->m, name))
		return NSAddressOfSymbol(nss);
	return 0;
	}

#define find_dlsym(a,b,c) (a = find_sym_addr(b,c))

 static void
dlclose(NS_pair *p)
{
	if (NSUnLinkModule(p->m, NSUNLINKMODULE_OPTION_NONE))
		NSDestroyObjectFileImage(p->ofi);
	free(p);
	}
#define NO_DLERROR
#else
#include "dlfcn.h"
typedef void *shl_t;
#define find_dlsym(a,b,c) (a = (Funcadd*)dlsym(b,c))
#ifdef sun
#ifndef RTLD_NOW
#define RTLD_NOW RTLD_LAZY
#endif
#endif /* sun */
#endif /* Old_APPLE */
#endif /* __hpux */
#endif /* WIN32 */

#ifdef __cplusplus
extern "C" {
#endif

#ifdef WIN32
 static int
wrong_bits(AmplExports *ae, char *name)
{
	FILE *f;
	IMAGE_DOS_HEADER dh;
	int rc;
	union { WORD w[2]; DWORD dw; } u;
	struct {
		IMAGE_FILE_HEADER ifh;
		IMAGE_OPTIONAL_HEADER ioh;
		} h;
#ifdef X64_bit_pointers
#define Bits_MAGIC 0x20b
#else
#define Bits_MAGIC 0x10b
#endif
	if (!(f = fopen(name, "rb")))
		return 1;
	rc = 0;
	 if (fread(&dh, sizeof(IMAGE_DOS_HEADER), 1, f) != 1
	 || dh.e_magic != IMAGE_DOS_SIGNATURE
	 || fseek(f, dh.e_lfanew, SEEK_SET)
	 || fread(&u, sizeof(u), 1, f) != 1
	 || u.dw != IMAGE_NT_SIGNATURE
	 || fread(&h, sizeof(h), 1, f) != 1
	 || h.ioh.Magic != Bits_MAGIC)
		rc = 1;
	fclose(f);
	return rc;
	}
#undef Bits_MAGIC
#endif

 static shl_t
dl_open(AmplExports *ae, char *name, int *warned, int *pns)
{
	FILE *f;
	char *d, *d0, *dz, *s;
	const char *cs;
	int ns;
	shl_t h;
#ifdef Old_APPLE
	NS_pair p;
#endif
	d = d0 = dz = 0;
	for(s = name; *s; ++s)
		switch(*s) {
		 case '.':
			d = s;
			break;
		 case '/':
#ifdef WIN32
		 case '\\':
#endif
			d = 0;
		 }
	ns = s - name;
	if (d
	 && d - name > 3
	 && d[-3] == '_') {
		if (d[-2] == BitsAlt[0]
		 && d[-1] == BitsAlt[1]) {
			d[-2] = Bits[0];
			d[-1] = Bits[1];
			dz = d;
			d = 0;
			}
		else if (d[-2] == Bits[0]
		 && d[-1] == Bits[1]) {
			dz = d;
			d = 0;
			}
		}
 tryagain:
#ifdef Old_APPLE
	NSObjectFileImageReturnCode irc;
	irc = NSCreateObjectFileImageFromFile(name,&p.ofi);
	h = 0;
	if (irc == NSObjectFileImageSuccess) {
		p.m = NSLinkModule(p.ofi, name,
			  NSLINKMODULE_OPTION_BINDNOW
			| NSLINKMODULE_OPTION_PRIVATE
			| NSLINKMODULE_OPTION_RETURN_ON_ERROR);
		if (!p.m)
			fprintf(stderr, "NSLinkModule(\"%s\") failed.\n", name);
		else {
			h = (NS_pair*)mymalloc(sizeof(NS_pair) + strlen(name) + 1);
			strcpy(p.name = (char*)(h+1), name);
			memcpy(h, &p, sizeof(NS_pair));
			}
		}
	else if (irc != NSObjectFileImageAccess)
		fprintf(stderr,
			"return %d from NSCreateObjectFileImageFromFile(\"%s\")\n",
			irc, name);
#else
#ifdef WIN32 /*{*/ /* make sure name is for the right number of bits */
	if (wrong_bits(ae, name))
		h = 0;
	else
#endif /*}*/
	h = dlopen(name, RTLD_NOW);
#endif
	if (!h) {
		if (d) {
			do s[3] = s[0]; while(--s >= d);
			d[0] = '_';
			d[1] = Bits[0];
			d[2] = Bits[1];
			d0 = d;
			d = 0;
			ns += 3;
			goto tryagain;
			}
		if (dz) {
			for(d = dz-3; (*d = *dz); ++d, ++dz);
			d = dz = 0;
			goto tryagain;
			}
		if (d0)
			for(s = d0; (s[0] = s[3]); ++s);
		if (!warned && (f = fopen(name,"rb"))) {
			fclose(f);
			if (file_kind(name) == 1) {
				*warned = 1;
#ifdef NO_DLERROR
				fprintf(Stderr, "Cannot load library \"%s\".\n", name);
#else
				fprintf(Stderr, "Cannot load library \"%s\"", name);
				cs = dlerror();
				fprintf(Stderr, cs ? ":\n%s\n" : ".\n", cs);
#endif
				}
			}
		}
	*pns = ns;
	return h;
	}

 static void
dl_close(void *h)
{
#ifdef CLOSE_AT_RESET
	first = 1;
#endif
	if (h)
		dlclose(h);
	}

 int
libload_ASL(AmplExports *ae, const char *s, int ns, int warn)
{
	Funcadd *fa;
	char buf0[2048], *buf;
	int ns1, rc, rcnf, warned;
	shl_t h;
	size_t n, nx;

	nx = 0;
	buf = buf0;
	if (!Abspath(s)) {
		if (!GetCurrentDirectory(sizeof(buf0),buf0))
			return 2;
		nx = strlen(buf0);
		}
	n = ns + sizeof(afdll) + nx + 3; /* +3 for inserting _32 or _64 */
	if (n > sizeof(buf0)) {
		buf = (char*)mymalloc(n);
		if (nx)
			memcpy(buf, buf0, nx);
		}
	if (nx)
		buf[nx++] = SLASH;
	strncpy(buf+nx, s, ns);
	buf[nx+ns] = 0;
	rc = warned = 0;
	rcnf = warn >> 1;
	warn &= 1;
	if ((h = dl_open(ae, buf, &warned, &ns1))) {
 found:
		if (find_dlsym(fa, h, FUNCADD)
		 || find_dlsym(fa, h, "funcadd")) {
#ifdef CLOSE_AT_RESET
			aflibname_ASL(ae,buf,buf+nx,ns1-nx,fa,0,dl_close,h);
				/* -DCLOSE_AT_RESET is for use in shared */
				/* libraries, such as MATLAB mex functions, */
				/* that may be loaded and unloaded several */
				/* times during execution of the program. */
#else
			aflibname_ASL(ae,buf,buf+nx,ns1-nx,fa,1,dl_close,h);
#endif
			}
		else {
			fprintf(stderr, "Could not find funcadd in %s\n", buf);
			dl_close(h);
			rc = 3;
			}
		}
	else if (warn) {
		if (!warned) {
			strcpy(buf+nx+ns, afdll);
			if ((h = dl_open(ae, buf, &warned, &ns1)))
				goto found;
			}
		if (warned)
			rc = 2;
		else
			goto notfound;
		}
	else {
 notfound:
		rc = rcnf;
		if (warn) {
			buf[nx+ns] = 0;
			if (file_kind(buf) == 2) {
				buf[nx+ns] = SLASH;
				fprintf(Stderr, "Cannot find library \"%s\".\n", buf);
				}
			else
				fprintf(Stderr, "Cannot find library \"%.*s\".\n", ns, s);
			}
		}
	if (buf != buf0)
		free(buf);
	return rc;
	}

 static int
libloop(AmplExports *ae, const char *s)
{
	const char *s1, *s2;
	int c, ns, rc;

	for(rc = 0;; s = s1) {
		while(*s <= ' ')
			if (!*s++)
				goto ret;
		if (*s == '"' || *s == '\'') {
			c = *s++;
			for(s1 = s; *s1 != c; ++s1)
				if (!*s1)
					goto ret;
			if (s1 == s)
				goto ret;
			s2 = s1++;
			}
		else {
			for(s1 = s; *++s1 >= ' '; );
			for(s2 = s1; s2[-1] == ' '; --s2);
			}
		ns = s2 - s;
		if (libload_ASL(ae, s, ns, 1))
			++rc;
		}
 ret:
	return rc;
	}

int n_badlibs_ASL;

 void
funcadd(AmplExports *ae)
{
	const char *s;
	int nb = 0;

	if (first) {
		first = 0;
		if ((s = i_option_ASL)) {
			if (!*s || (*s == '-' && !s[1]))
				return;
			nb += libloop(ae, s);
			}
		else
			nb = libload_ASL(ae, afdll+1, (int)sizeof(afdll)-2, 0);
		}
	n_badlibs_ASL = nb;
	}

#ifdef __cplusplus
}
#endif

#endif /* NO_FUNCADD */
