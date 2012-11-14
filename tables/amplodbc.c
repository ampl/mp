/****************************************************************
Copyright (C) 2000-2001 Lucent Technologies
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

 static char Version[] = "\n@(#) AMPL ODBC driver, version 20121108.\n";

#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

#include <sql.h>
#include <sqlext.h>
#include <stdlib.h>
#include <string.h>
#include "arith.h"	/* for LONG_LONG_POINTERS */
#include "funcadd.h"

#ifndef SQL_NO_DATA
#define SQL_NO_DATA SQL_NO_DATA_FOUND	/* VC++ 4 */
#endif

#define UNUSED(x) ((void)(1 ? 0 : (x)))

static unsigned char tolc[256];	/* for converting to lower case */

#ifndef DEFAULT_MAXLEN
#define DEFAULT_MAXLEN 480
#endif

#define CC (const char*)
#define UC (UCHAR*)
#define TM(len) (*ae->Tempmem)(TI->TMI,len)

 typedef struct
UnknownType {
	struct UnknownType *next;
	SDWORD type;
	int mytype;
	} UnknownType;

 typedef struct
HInfo {
	AmplExports *AE;
	TableInfo *TI;
	HENV	env;
	HDBC	hc;
	HSTMT	hs;
	int	verbose;
	unsigned int maxlen;
	char	*s, *se; /* for use in Adjust_ampl_odbc */
	int	totbadtimes;
	int	firstbadtimes;
	int	firstbadcol;
	int	nsmix;
	int	ntimes;
	int	oldquotes;
	int	sqldb;
	int	thiscol;
	int	wrmode;
	long	firstbadrow;
	int	*coltypes;
	int	*colperm;
	TIMESTAMP_STRUCT *ts;
	/* for scrunch and ascrunch... */
	char	*Missing;
	real	*dd;
	UnknownType *ut;
	SQLCHAR quote; /* A quote character for identifiers. */
	char **quoted_colnames; /* Quoted column names. */
	} HInfo;

 enum { /* for wrmode */ wr_drop=0, wr_append=1 };

#ifndef NO_Adjust_ampl_odbc
 static int Adjust_ampl_odbc(HInfo*, char*, TIMESTAMP_STRUCT****, int**, int**, int*, int**);
#endif

 static void *
Malloc(AmplExports *ae, size_t len)
{
	void *rv = malloc(len);
	if (!rv) {
		fprintf(Stderr, "malloc(%lu) failure!\n", (unsigned long)len);
		exit(1);
		}
	return rv;
	}

#define offset_of(t,x) ((size_t)&((t*)0)->x)

 static int
prc(HInfo *h, char *who, int i)
{
	typedef struct HandleStuff {
		const char *desc;
		int hkind;
		size_t hoff;
		} HandleStuff;
	AmplExports *ae;
	SDWORD native_errno;
	SQLSMALLINT errmsglen;
	SWORD emlen;
	UCHAR *errmsg, errmsg0[SQL_MAX_MESSAGE_LENGTH], sqlstate[64];
	int rv;
#if 0
	HandleStuff *hs, *hse;
	SQLINTEGER k, ne;
	SQLPOINTER sptr;
	SQLSMALLINT elen, htype, j;
	static HandleStuff HS[3] = {
		{ "SQL_HANDLE_STMT", SQL_HANDLE_STMT, offset_of(HInfo, hs) },
		{ "SQL_HANDLE_DBC",  SQL_HANDLE_DBC,  offset_of(HInfo, hc) },
		{ "SQL_HANDLE_ENV",  SQL_HANDLE_ENV,  offset_of(HInfo, env) }};
#endif

	if (i == SQL_SUCCESS)
		return 0;
	rv = 1;
	if (i == SQL_SUCCESS_WITH_INFO) {
		if (!h->verbose)
			return 0;
		rv = 0;
		}
	errmsg = errmsg0;
	errmsglen = sizeof(errmsg0);
	ae = h->AE;
	printf("%s returned %d\n", who, i);
#if 0
	elen = 0;
	k = 0;
	for(hs = HS, hse = HS + 3; hs < hse; ++hs) {
		htype = hs->hkind;
		sptr = (SQLPOINTER)((char*)h + hs->hoff);
		i = SQLGetDiagField(htype, sptr, 0, SQL_DIAG_NUMBER,
			&k, 0, &elen);
		if (i == SQL_SUCCESS)
			break;
		printf("SQLGetDiagField(%s) returned %d\n", hs->desc, i);
		}
	if (hs >= hse)
		goto use_SQLError;
	if (k <= 0) {
		printf("Surprising number of diagnostic records %d reported by "
			"SQLGetDiagField(...,SQL_DIAG_NUMBER,...)\n", (int)k);
		goto use_SQLError;
		}
	for(j = 1; j <= k; ++j) {
		while((i = SQLGetDiagRec(SQL_HANDLE_STMT, h->hs, j, sqlstate, &ne,
				errmsg, errmsglen, &elen)) != SQL_SUCCESS) {
			if (i != SQL_SUCCESS_WITH_INFO) {
				printf("Surprise return %d from SQLGetDiagRec\n", i);
				goto use_SQLError;
				}
			if (elen < errmsglen)
				break;
			do errmsglen <<= 1; while (elen >= errmsglen);
			if (errmsg != errmsg0)
				free(errmsg);
			errmsg = UC Malloc(ae, errmsglen);
			printf("Trying errmsglen = %d\n", (int)errmsglen);
			}
		if (j > 1)
			printf("Diagnostic record %d:\n", (int)j);
		printf("SQL state %.5s: %s\n", sqlstate, errmsg);
		}
	goto done;
 use_SQLError:
#endif
	i = SQLError(h->env, h->hc, h->hs, sqlstate, &native_errno,
	    errmsg, sizeof(errmsg0), &emlen);
	if (i != SQL_SUCCESS && i != SQL_SUCCESS_WITH_INFO)
		printf("SQLError returned %d\n", i);
	else {
		printf("sqlstate = \"%s\"\n", sqlstate);
		printf("errmsg = \"%s\"\n", errmsg);
		printf("native_errno = %d\n", native_errno);
		}
/* done: */
	if (errmsg != errmsg0)
		free(errmsg);
	return rv;
	}

 static int
prcnr(HInfo *h, char *who, int i)	/* variant for SQLGetNumResultCols */
{
	AmplExports *ae;
	SWORD emlen;
	SDWORD native_errno;
	UCHAR errmsg[SQL_MAX_MESSAGE_LENGTH], sqlstate[64];
	int j, rv;

	if (i == SQL_SUCCESS)
		return 0;
	rv = 1;
	if (i == SQL_SUCCESS_WITH_INFO) {
		if (!h->verbose)
			return 0;
		rv = 0;
		}
	ae = h->AE;
	j = SQLError(h->env, h->hc, h->hs, sqlstate, &native_errno,
			errmsg, sizeof(errmsg), &emlen);
	if (j == SQL_SUCCESS && !strcmp((char*)sqlstate,"S0002") && !h->verbose)
		return 1;
	printf("%s returned %d\n", who, i);
	if (j != SQL_SUCCESS)
		printf("SQLError returned %d\n", j);
	else {
		printf("sqlstate = \"%s\"\n", sqlstate);
		printf("errmsg = \"%s\"\n", errmsg);
		printf("native_errno = %d\n", native_errno);
		}
	return rv;
	}

 static void
freestmt(HSTMT *hs)
{
	if (*hs != SQL_NULL_HSTMT) {
		SQLFreeStmt(*hs, SQL_DROP);
		*hs = SQL_NULL_HSTMT;
		}
	}

 static int
match(const char *dext0, UCHAR *ext, UCHAR *eend)
{
	UCHAR *dext = UC dext0;
	int c, c1;
	for(;;) {
		c = tolc[*dext++];
		if (ext == eend) {
			if (!c)
				break;
			return 0;
			}
		c1 = tolc[*ext++];
		if (c != c1)
			return 0;
		if (c == 0)
			break;
		}
	return 1;
	}

 static UCHAR *
firstext(UCHAR *s, int *n_ext, int *elen, UCHAR **sp)
{
	UCHAR *rv;
	int len, n;

	len = n = 0;
	rv = 0;

	for(;;) {
		if (!*s)
			goto done;
		if (match("FileExtns=",s,s+10)) {
			s += 10;
			break;
			}
		while(*s++);
		}
	if (sp)
		*sp = s;
	for(;;) {
		while(*s != '.') {
			if (!*s++)
				goto done;
			}
		if (!rv)
			rv = s;
		n++;
		while(*++s != ',') {
			if (!*s)
				goto done;
			len++;
			}
		}
 done:
	*elen = len;
	*n_ext = n;
	return rv;
	}

 typedef struct
DRV_desc {
	struct DRV_desc *next;
	UCHAR *driver;
	UCHAR *dsn;
	UCHAR *ext;
	char *ntype;	/* SQL numeric  type name */
#if 0
	char *stype;	/* SQL string   type name */
#endif
	char *ttype;	/* SQL datetime type name */
	unsigned int tprec;
	} DRV_desc;

static UCHAR no_dsn[] = " ";
static DRV_desc *ds0, **ds1p;
static int nds1p;

 static void
mfree(void *v) { free(v); }

 static DRV_desc*
get_ds0(HInfo *h)
{
	/* Return pointer to first DRV_desc -- generate list first time. */

	AmplExports *ae = h->AE;
	DRV_desc *ds, **ds1, **dsp, *dsx;
	SWORD attr_len, desc_len;
	UCHAR attr[1024], desc[128], *drv1, *s, *s1, *sx;
	UWORD dir;
	int elen, elt, i, n, nd, nt, verbose;
	static int beenhere;

	if (!(ds = ds0)) {
		if (beenhere++)
			return ds;
		dsp = &ds0;
		if ((verbose = h->verbose))
			printf("ODBC drivers...\n");
		elt = nd = nt = 0;
		desc_len = sizeof(desc);
		attr_len = sizeof(attr);
		for(dir = SQL_FETCH_FIRST;; dir = SQL_FETCH_NEXT) {
			memset(desc, 0, desc_len);
			memset(attr, 0, attr_len);
			i = SQLDrivers(h->env, dir, desc, sizeof(desc), &desc_len,
						attr, sizeof(attr)-2, &attr_len);
			if (i != SQL_SUCCESS && i != SQL_SUCCESS_WITH_INFO)
				break;
			if (desc_len > sizeof(desc)) {
				printf("Surprise in get_ds0: desc_len = %d\n", desc_len);
				desc_len = sizeof(desc);
				}
			if (attr_len > sizeof(attr) - 2) {
				printf("Surprise in get_ds0: attr_len = %d\n", attr_len);
				attr_len = sizeof(attr);
				}
			if (verbose)
				printf("\t\"%s\"", desc);
			if (!(s = firstext(attr, &n, &elen, &s1))) {
				if (verbose)
					printf(" -- no associated extensions.\n");
				continue;
				}
			if (verbose)
				printf(" for %s\n", s1);
			elt += elen + desc_len + 1;
			nd++;
			nt += n;
			}
		if (!(nds1p = nd))
			return ds;
		dsx = (DRV_desc*)Malloc(ae, nt*(sizeof(DRV_desc) + 1)
					+ nd*sizeof(DRV_desc*) + elt);
		at_exit(mfree, dsx);
		ds1 = ds1p = (DRV_desc**)(dsx + nt);
		sx = UC (ds1 + nd);
		for(dir = SQL_FETCH_FIRST;; dir = SQL_FETCH_NEXT) {
			memset(desc, 0, desc_len);
			memset(attr, 0, attr_len);
			i = SQLDrivers(h->env, dir, desc, sizeof(desc), &desc_len,
						attr, sizeof(attr)-2, &attr_len);
			if (i != SQL_SUCCESS && i != SQL_SUCCESS_WITH_INFO)
				break;
			if (!(s = firstext(attr, &n, &elen, 0)))
				continue;
			s1 = desc;
			drv1 = sx;
			while((*sx++ = *s1++));
			*ds1++ = dsx;
			for(;;) {
				ds = dsx++;
				ds->driver = drv1;
				ds->dsn = 0;
				ds->ntype = 0;
				*dsp = ds;
				dsp = &ds->next;
				ds->ext = sx;
				while((*sx = *++s) && *s != ',')
					sx++;
				*sx++ = 0;
				if (!--n)
					break;
				while(*++s != '.');
				}
			}
		*dsp = 0;
		ds = ds0;
		}
	return ds;
	}

 static int
dr_cmp(const void *a, const void *b, void *v)
{
	UNUSED(v);
	return strcmp((char*)(*(DRV_desc**)a)->driver, (char*)(*(DRV_desc**)b)->driver);
	}

 static DRV_desc*
desc_search(UCHAR *desc)
{
	/* binary search for desc */

	DRV_desc *ds, **ds1, **dsp;
	int i, n, n1;

	n = nds1p;
	dsp = ds1p;
	while(n > 0) {
		ds1 = dsp + (n1 = n >> 1);
		ds = *ds1;
		if (!(i = strcmp((char*)desc, (char*)ds->driver)))
			return ds;
		if (i < 0)
			n = n1;
		else {
			n -= n1 + 1;
			dsp = ds1 + 1;
			}
		}
	return 0;
	}

 static void
gen_dsn(HInfo *h)
{
	AmplExports *ae;
	DRV_desc *ds;
	SWORD desc_len, dsn_len;
	UCHAR desc[256], *driver, dsn[SQL_MAX_DSN_LENGTH+1], *s, *s0, *s1;
	UWORD dir;
	int i, tlen;
	static UCHAR yes_dsn[] = "<<yes>>";

	ae = h->AE;
	if (nds1p > 1)
		qsortv(ds1p, nds1p, sizeof(DRV_desc*), dr_cmp, 0);
	tlen = 0;
	for(dir = SQL_FETCH_FIRST;; dir = SQL_FETCH_NEXT) {
		i = SQLDataSources(h->env, dir, dsn, sizeof(dsn), &dsn_len,
					desc, sizeof(desc), &desc_len);
		if (i != SQL_SUCCESS && i != SQL_SUCCESS_WITH_INFO)
			break;
		if ((ds = desc_search(desc)) && !ds->dsn) {
			tlen += dsn_len + 1;
			ds->dsn = yes_dsn;
			}
		}
	if (!tlen)
		return;
	s = Malloc(ae, tlen);
	at_exit(mfree, s);
	for(dir = SQL_FETCH_FIRST;; dir = SQL_FETCH_NEXT) {
		i = SQLDataSources(h->env, dir, dsn, sizeof(dsn), &dsn_len,
					desc, sizeof(desc), &desc_len);
		if (i != SQL_SUCCESS && i != SQL_SUCCESS_WITH_INFO)
			break;
		if ((ds = desc_search(desc)) && ds->dsn == yes_dsn) {
			ds->dsn = s0 = s;
			for(s1 = dsn; (*s++ = *s1); s1++);
			driver = ds->driver;
			while((ds = ds->next) && ds->driver == driver)
				ds->dsn = s0;
			}
		}
	}

 static char*
get_dsn(HInfo *h, DRV_desc *ds)
{
	static int beenhere;

	if (ds->dsn == no_dsn)
		return 0;
	if (!beenhere++)
		gen_dsn(h);
	return (char*)ds->dsn;
	}

 static DRV_desc *
ext_to_drv(UCHAR *ext, UCHAR *eend, HInfo *h)
{
	/* Return ds and driver corresponding to extension ext. */

	DRV_desc *ds;

	for(ds = get_ds0(h); ds; ds = ds->next)
		if (match(CC ds->ext, ext, eend))
			break;
	return ds;
	}

 static char *
copy(char *a, char *b)
{
	while((*a = *b++))
		a++;
	return a;
	}

 static int
dsncompar(const void *a, const void *b, void *v)
{
	UNUSED(v);
	return strcmp((char*)(*(DRV_desc**)a)->ext,
		(char*)(*(DRV_desc**)b)->ext);
	}

 static void
bad_ending(char *fmt, char *dsname, HInfo *h)
{
	AmplExports *ae = h->AE;
	DRV_desc *ds = ds0, *ds1, **z, **z1;
	TableInfo *TI = h->TI;
	char *plural, *s;
	int n;
	size_t L;

	if (!ds) {
		TI->Errmsg = "ODBC is not available.";
		return;
		}
	get_dsn(h,ds); /* force gen_dsn() if necessary */
	for(n = 0; ds; ds = ds->next)
		if (ds->dsn)
			n++;
	z = z1 = (DRV_desc**)TM((n+1)*sizeof(DRV_desc*));
	for(ds = ds0; ds; ds = ds->next)
		if (ds->dsn)
			*z1++ = ds;
	*z1 = 0;
	qsortv(z, n, sizeof(DRV_desc*), dsncompar, 0);
	z1 = z;
	n = 1;
	ds = *z;
	while((ds1 = *++z1))
		if (strcmp((char*)ds->ext, (char*)ds1->ext)) {
			z[n++] = ds1;
			ds = ds1;
			}
	z[n] = 0;
	L = strlen(fmt) + strlen(dsname) + 16; /*elbow room*/
	for(z1 = z; (ds = *z1++);)
		L += strlen(CC ds->ext) + 3;
	if (h->verbose) {
		L += 210;
		for(ds = ds0; ds; ds = ds->next) {
			L += strlen(CC ds->driver) + strlen(CC ds->ext) + 16;
			if (ds->dsn)
				L += strlen(CC ds->dsn) + 16;
			while((ds1 = ds->next) && ds->driver == ds1->driver) {
				ds = ds1;
				L += strlen(CC ds->ext);
				}
			}
		}
	TI->Errmsg = s = (char*)TM(L);
	s += sprintf(s, fmt, dsname);
	for(z1 = z; (ds = *z1++);)
		s += sprintf(s, "\n\t.%s", ds->ext);
	if (!h->verbose)
		return;
	s += sprintf(s, "\n\n%s\n%s\n%s\n",
	 "Possible DRIVER= and sample DSN= strings and associated extensions",
	 "(.ext values) for use with \"DRIVER=...;DBQ=somthing.ext\"",
	 "or with \"DSN=...;DBQ=something.ext\":");
	for(ds = ds0; ds; ds = ds->next) {
		s += sprintf(s, "\n\tDRIVER=%s;\n", ds->driver);
		if (ds->dsn)
			s += sprintf(s, "or\n\tDSN=%s;\n", ds->dsn);
		plural = "";
		if ((ds1 = ds->next) && ds1->driver == ds->driver)
			plural = "s";
		s += sprintf(s, "with extension%s\n\t.%s\n", plural, ds->ext);
		while((ds1 = ds->next) && ds1->driver == ds->driver) {
			ds = ds1;
			s += sprintf(s, "\t.%s\n", ds->ext);
			}
		}
	}

 static int
ODBC_check(AmplExports *AE, TableInfo *TI, HInfo *h)
{
	char *s;
	int i = TI->nstrings;

	if (i < 1 || (strcmp(s = TI->strings[0], "ODBC") && strcmp(s, "odbc")))
		return DB_Refuse;
	if (i < 2 || i > 5) {
		TI->Errmsg =
		  "AMPL ODBC handler: expected 2-8 strings before \":[...]\":\n\
  'ODBC', connection_spec [ext_name] [optional_strings]\n\n"
	"For more details, use the AMPL command\n\n\tprint _handler_desc['odbc'];";
		return DB_Error;
		}
	memset(h, 0, sizeof(HInfo));
	h->AE = AE;
	h->TI = TI;
	h->Missing = TI->Missing;
	h->maxlen = DEFAULT_MAXLEN;

	if (!tolc[1]) {
		for(i = 0; i < 256; i++)
			tolc[i] = i;
		for(i = 'A'; i <= 'Z'; i++)
			tolc[i] = i + 'a' - 'A';
		}

	return 0;
	}

typedef int Compare(const void*, const void*, void*);

 static int
compar(const void *a, const void *b, void *v)
{
	return strcmp(((char**)v)[*(int*)a], ((char**)v)[*(int*)b]);
	}

 static void
strsort(AmplExports *ae, int n, void *b, Compare *cmp, int *z)
{
	int i;
	for(i = 0; i < n; i++)
		z[i] = i;
	if (n > 1)
		qsortv(z, n, sizeof(int), cmp, b);
	}

 static int *
ct_alloc(HInfo *h, int **parray)
{
	AmplExports *ae = h->AE;
	TableInfo *TI = h->TI;
	int k, n, *rv, *z, *ze;

	n = TI->arity + TI->ncols;
	k = (rv = h->coltypes) ? 0 : n;
	if (parray)
		k += n;
	z = (int*)TM(k*sizeof(int));
	if (!rv) {
		h->coltypes = rv = z;
		for(ze = z + n; z < ze; ++z)
			*z = 0;
		}
	if (parray) {
		*parray = h->colperm = z;
		strsort(ae, n, TI->colnames, compar, z);
		}
	return rv;
	}

 static int
strcmpv(const void *a, const void *b, void *v)
{
	int i, j;
	const char *s, *t;
	UNUSED(v);
	s = (const char*)a;
	t = (const char*)b;
	for(;;) {
		i = *s++;
		j = *t++;
		if (i != j)
			return i - j;
		if (!i)
			break;
		}
	return 0;
	}

 static int
get_times(HInfo *h, char *names)
{
	AmplExports *ae = h->AE;
	TableInfo *TI = h->TI;
	char **cn, *s, *s1, **times;
	int cc, *ct, i, j, k, n, nc, m, *p, q;

	for(cc =n = 0, s = names;;) {
 next_name:
		while(*s <= ' ')
			if (!*s++)
				goto have_n;
		if (cc) {
			cc = 0;
			if (*s == ',') {
				s++;
				goto next_name;
				}
			}
		n++;
		while((q = *s) > ' ') {
			s++;
			switch(q) {
			  case '"':
			  case '\'':
				for(;;s++) {
					if (*s == q) {
						if (*++s != q)
							break;
						}
					else if (!*s) {
						TI->Errmsg = (char*)
						 "Ill-formed quoted string in \"times=...\"";
						return 1;
						}
					}
				break;
			  case ',':
				goto next_name;
			  }
			}
		cc = 1;
		}
 have_n:
	if (!n)
		return 0;
	ct = ct_alloc(h, &p);
	times = (char**)TM(n*sizeof(char*) + (s-names));
	s1 = (char*)(times + n);
	for(cc = n = 0, s = names;;) {
 next_name1:
		while(*s <= ' ')
			if (!*s++)
				goto done;
		if (cc) {
			cc = 0;
			if (*s == ',') {
				s++;
				goto next_name1;
				}
			}
		times[n++] = s1;
		while((q = *s) > ' ') {
			s++;
			switch(q) {
			  case '"':
			  case '\'':
				for(;;s++) {
					if (*s == q) {
						if (*++s != q)
							break;
						*s1++ = q;
						}
					else
						*s1++ = *s;
					}
				break;
			  case ',':
				goto end_name;
			  default:
				*s1++ = q;
			  }
			}
		cc = 1;
 end_name:
		*s1++ = 0;
		}
 done:
	if (n > 1)
		qsortv(times, n, sizeof(char*), strcmpv, 0);
	cn = TI->colnames;
	nc = TI->arity + TI->ncols;
	for(i = j = m = 0; i < n; i++) {
		s = times[i];
		while((k = strcmp(cn[j],s)) < 0) {
			if (++j >= nc) {
 notfound:
				sprintf(TI->Errmsg = (char*)TM(strlen(s)+32),
					"times column \"%s\" not found.", s);
				return 1;
				}
			}
		if (k > 0)
			goto notfound;
		if (ct[j] == 0) {
			ct[j] = 2;
			++m;
			}
		}
	h->ntimes += m;
	return 0;
	}

 static int
not_found(AmplExports *ae, TableInfo *TI, char *dsn, FILE **fp)
{
	FILE *f, **fp1;
	char *mode, *tname;

	if ((fp1 = fp))
		mode = "r";
	else {
		fp1 = &f;
		mode = "rb";
		}
	if ((*fp1 = fopen(dsn, mode))) {
		if (!fp)
			fclose(f);
		return 0;
		}
	tname = TI->nstrings >= 3 ? TI->strings[2] : TI->tname;
	sprintf(TI->Errmsg = (char*)TM(strlen(tname) + strlen(dsn) + 40),
		"Cannot open file %s for table %s.", dsn, tname);
	return 1;
	}

 enum { Bsize = 1024 };
 typedef union {
		double align;
		char clbuf[Bsize];
		long L;
		} uBuf;

 static DRV_desc *
ds_alloc(HInfo *h, char *dsname)
{
	AmplExports *ae = h->AE;
	DRV_desc *ds;
	TableInfo *TI = h->TI;
	size_t len;

	len = strlen(dsname) + 1;
	ds = (DRV_desc*)TM(sizeof(DRV_desc)+len);
	memset(ds, 0, sizeof(DRV_desc));
	memcpy(ds->driver = UC (ds+1), dsname, len);
	ds->dsn = no_dsn;
	return ds;
	}

 static DRV_desc *
dsname_ds(HInfo *h, char *dsname)
{
	/* return ds for dsname */

	DRV_desc *ds;
	unsigned char *a, *b;

	for(ds = get_ds0(h);; ds = ds->next) {
		if (!ds) {
			ds = ds_alloc(h, dsname);
			break;
			}
		a = ds->driver;
		b = (unsigned char*)dsname;
		while(tolc[*a] == tolc[*b++])
			if (!*a++)
				goto ret;
		}
 ret:
	return ds;
	}

 static int
parse_dsn_file(HInfo *h, char **dsnp)
{
	AmplExports *ae;
	FILE *f;
	TableInfo *TI;
	char buf[Bsize], *dsn, *s, *s0, *s1;
	int rv;
	long L;

	rv = 0;
	if (not_found(ae = h->AE, TI = h->TI, dsn = *dsnp, &f))
		goto ret;

	/* Search for DBQ= and FIL= lines */

	fseek(f, 0L, SEEK_END);
	L = ftell(f);
	fseek(f, 0L, SEEK_SET);
	if (fgets(buf, sizeof(buf), f) && !strncmp(buf,"[ODBC]",6)) {
		s = s0 = (char*)TM(L);
		while(fgets(buf, sizeof(buf), f)) {
			for(s1 = buf; (*s = *s1) && *s1 != '\n'; s1++)
				s++;
			*s++ = ';';
			}
		if (s > s0) {
			rv = 1;
			*dsnp = s0;
			}
		else
			sprintf(TI->Errmsg = (char*)TM(strlen(dsn) + 40),
				"File %s only had \"[ODBC]\".", dsn);
		}
	else
		sprintf(TI->Errmsg = (char*)TM(strlen(dsn) + 40),
			"File %s did not start with \"[ODBC]\".", dsn);
	fclose(f);
 ret:
	return rv;
	}

 static DRV_desc *
conn_ds(HInfo *h, UCHAR *cs)
{
	/* search connection string for DRIVER=... or DSN=..., return corresp. ds */

	DRV_desc *ds;
	UCHAR c, *e;

	for(;;cs++) {
		while(*cs != 'D')
			if (!*cs++)
				return 0;
		if (cs[1] == 'S' && cs[2] == 'N' && cs[3] == '=') {
			cs += 4;
			break;
			}
		if (match("river=", cs+1, cs+7)) {
			cs += 7;
			break;
			}
		}
	for(e = cs; *e && *e != ';'; e++);
	c = *e;
	*e = 0;
	ds = dsname_ds(h, (char*)cs);
	*e = c;
	return ds;
	}

 static char*
get_ext(char *s)
{
	char *t;

 top:
	while(*s) {
		if (*s++ == '.')
			for(t = s;;)
				switch(*t++) {
				  case 0:
					return s;
				  case '/':
				  case ':':
				  case '\\':
					s = t;
					goto top;
				  case '.':
					s = t;
				  }
		}
	return 0;
	}

 static char*
fully_qualify(char *dsname, char *buf, size_t len)
{
	DWORD n;
	int c;
	size_t Ldsn;

	c = *dsname;
	if (c == '/' || c == '\\')
		return dsname;
	if (((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z')) && dsname[1] == ':')
		return dsname;
	Ldsn = strlen(dsname);
	if (!getcwd(buf, len))
          return dsname;
	n = strlen(buf);
	if (n + Ldsn + 1 >= len)
		return dsname;
	buf[n++] = '\\';
	strcpy(buf+n, dsname);
	return buf;
	}

 static struct {
	HWND hw;
	int score;
	} winfo;

#ifdef _WIN32

 static char *wprocnames[] = {
	"start menu",
	"program manager",
	"ampl_lic",
	"ampl",
	"command prompt",
	"sw:",
	0};

 static int
match1(char *s, char *t)
{
	while(*t) {
		if (*s++ != *t++)
			return 0;
		}
	return *s <= ' ';
	}

 static BOOL CALLBACK
mywproc(HWND hw, LPARAM all)
{
	char buf[2048], *s, *se;
	int i, len;

	len = GetWindowText(hw, buf, sizeof(buf));
	if (len <= 0)
		goto done;
	se = buf + 15;
	if (len < 15)
		se = buf + len;
	for(s = buf; s < se; ++s)
		*s = tolc[(int)*s];
	for(i = 0; (s = wprocnames[i++]); ) {
		if (match1(buf,s)) {
			if (winfo.score > i) {
				winfo.score = i;
				winfo.hw = hw;
				if (i == 1)
					return FALSE;
				break;
				}
			}
		}
 done:
	return TRUE;
	}

 static void
hw_get(AmplExports *ae)
{
	const char *s;
	/* work around MS bugs with stdlib (missing routines) */
	int c;
#ifdef LONG_LONG_POINTERS
	unsigned long long x;
#else
	unsigned long x;
#endif

	winfo.score = 9999;
	if ((s = getenv("sw_HWND"))) {
		for(x = 0; (c = *s); ++s) {
			if (c >= '0' && c <= '9')
				c -= '0';
			else if (c >= 'a' && c <= 'f')
				c += 10 - 'a';
			else if (c >= 'A' && c <= 'F')
				c += 10 - 'A';
			else
				continue;
			x <<= 4;
			x |= c;
			}
		winfo.hw = (HWND)x;
		/*DEBUG*/ /*printf("Got $sw_HWND = #%x\n", wi.hw);*/
		}
	else {
		winfo.hw = 0;
		EnumWindows(mywproc, 0);
		}
	}

#endif /* _WIN32 */

 static void
colname_adjust(HInfo *h, TableInfo *TI)
{
	/* adjust names in TI->colnames of the form */
	/* Time:xyz, Mixed:xyz, or Strcol:xyz to xyz . */

	UCHAR *s, *s1;
	char **cn;
	int *ct, i, k, nc, nt;

	cn = TI->colnames;
	nc = TI->arity + TI->ncols;
	h->coltypes = h->colperm = ct = 0;
	for(i = nt = 0; i < nc; ++i) {
		s = UC cn[i];
		if (match("Time:", s, s1 = s+5)) {
			k = 2;
			++nt;
			goto cn_adj;
			}
		if (match("Strcol:", s, s1 = s+7)) {
			k = 1;
			goto cn_adj;
			}
		if (match("Mixed:", s, s1 = s+6)) {
			k = 3;
 cn_adj:
			if (!ct)
				ct = ct_alloc(h,0);
			ct[i] = k;
			cn[i] = (char*)s1;
			}
		}
	h->ntimes = nt;
	}

/* Quotes a table or column name. */
static char *quote_name(HInfo *h, char *name)
{
	char quote = h->quote;
	size_t quoted_name_length = 0;
	int num_quotes = 0;
	char *quoted_name = 0;
	char *s = 0;
	if (quote == ' ')
		return name;
	for (s = name; *s; ++s) {
		if (*s == quote)
			++num_quotes;
	}
	quoted_name_length = strlen(name) + num_quotes + 2;
	quoted_name = (*h->AE->Tempmem)(h->TI->TMI, quoted_name_length + 1);
	quoted_name[0] = quote;
	for (s = quoted_name + 1; *name; ++name, ++s) {
		*s = *name;
		if (*name == quote)
			*++s = *name;
	}
	quoted_name[quoted_name_length - 1] = quote;
	quoted_name[quoted_name_length] = '\0';
	return quoted_name;
}

 static char*
Connect(HInfo *h, DRV_desc **dsp, int *rc, char **sqlp)
{
	AmplExports *ae = h->AE;
	DRV_desc *ds;
	SWORD cs_len;
	TableInfo *TI = h->TI;
	UCHAR cs[512];
	char dsbuf[1024];
	char *dreq, *driver, *dsname, *dsname0, *e1, *ee, *s, *s1, *s2, **strs, *tname;
	int i, j, nstr, verbose, wantretry;
	real t;
	unsigned int ui;
	int dsn = 0;

	*dsp = 0;
	h->env = SQL_NULL_HENV;
	h->hc = SQL_NULL_HDBC;
	h->hs = SQL_NULL_HSTMT;
	if (SQLAllocEnv(&h->env) != SQL_SUCCESS) {
		TI->Errmsg = "SQLAllocEnv failed!";
 eret:
		*rc = DB_Error;
		return 0;
		}
	if (prc(h, "SQLAllocConnect", SQLAllocConnect(h->env, &h->hc))) {
 unexpected:
		TI->Errmsg = "Unexpected ODBC failure.";
		goto eret;
		}
	strs = TI->strings;
	dsname = dsname0 = strs[1];
	tname = 0;
	verbose = 0;
	h->wrmode = wr_drop;
	h->nsmix = 1;
	h->oldquotes = 0;
	if ((nstr = TI->nstrings) < 3) {
		s = 0; /* shut up erroneous warning */
		nstr = -1;
		goto know_verbose;
		}
	s = strs[2];
	if (strchr(s, '=') || match("verbose", UC s, UC s+7)) {
		strs += 2;
		nstr -= 2;
		}
	else {
		tname = s;
		strs += 3;
		nstr -= 3;
		if (nstr <= 0)
			goto know_verbose;
		}
	if (nstr > 1)
		qsortv(strs, nstr, sizeof(char*), strcmpv, 0);

	s = *strs++;
	--nstr;
	if (sqlp && match("SQL=", UC s, UC s + 4)) {
		*sqlp = s + 4;
		if (nstr <= 0)
			goto know_verbose;
		s = *strs++;
		--nstr;
		}
	if (match("maxlen=", UC s, UC s + 7)) {
		s1 = s2 = "";
		if ((ui = t = strtod(s+7,&s2)) > 0 && !*s2 && t == ui)
			h->maxlen = ui;
		else {
			sprintf(TI->Errmsg = TM(strlen(s) + 72),
				"Inappropriate value in \"%s\":\n\t%s",
				s, "expected a positive decimal integer.");
			goto eret;
			}
		if (nstr <= 0)
			goto know_verbose;
		s = *strs++;
		--nstr;
		}
	if (match("nsmix=", UC s, UC s + 6)) {
		s1 = s + 6;
		j = 2;
		switch(*s1) {
		 case '!':
			if (!s1[1]) {
				j = 1;
				h->oldquotes = 1;
				}
			break;
		 case '*':
			if (!s1[1])
				j = 1;
			break;
		 case '-':
			if (s1[1])
				break;
		 case 0:
			j = 0;
			break;
		 default:
			sprintf(TI->Errmsg = TM(strlen(s) + 72),
				"Unexpected nsmix value in \"%s\"", s);
			goto eret;
		 }
		h->nsmix = j;
		if (nstr <= 0)
			goto know_verbose;
		s = *strs++;
		--nstr;
		}
	if (match("time=", UC s, UC s + 5)) {
		if (get_times(h, s+5))
			goto eret;
		if (nstr <= 0)
			goto know_verbose;
		s = *strs++;
		--nstr;
		}
	if (match("write=", UC s, UC s + 6)) {
		s += 6;
		if (match("append", UC s, UC s+6))
			h->wrmode = wr_append;
		else if (match("drop", UC s, UC s+4))
			h->wrmode = wr_drop;
		else {
			sprintf(TI->Errmsg = TM(strlen(s) + 100),
				"Inappropriate value in \"write=%s\":\n"
				"expected \"write=append\" or \"write=drop\".",
				s);
			goto eret;
			}
		if (nstr <= 0)
			goto know_verbose;
		s = *strs++;
		--nstr;
		}
	if (match("verbose=", UC s, UC s + 8)) {
		verbose = (int)strtol(CC s + 8, &s2, 10);
		s = *strs++;
		--nstr;
		if (verbose & 1)
			goto verbose1;
		goto know_verbose;
		}
	if (!match("verbose", UC s, UC s + 7)) {
 badopt:
		sprintf(TI->Errmsg = TM(strlen(s) + 140),
			"Unexpected option \"%s\".\n"
			"For a summary of possible options, issue the AMPL command\n\n"
			"\tprint _handler_desc['odbc'];\n", s);
		goto eret;
		}
	verbose = 1;
 verbose1:
	printf("%s", Version+6);
 know_verbose:
	if (nstr > 0)
		goto badopt;
	if (!tname)
		tname = TI->tname;
	wantretry = 1;
	h->sqldb = verbose & 2;
	h->verbose = verbose &= 1;
#ifdef _WIN32
	if (!winfo.score)
		hw_get(ae);
#endif
	if ((dsn = match("DSN=", UC dsname, UC dsname+4))
	 || match("DRIVER=",UC dsname, UC dsname+7)) {
		SQLUSMALLINT completion = 0;
 try_dsname:
		completion = SQL_DRIVER_COMPLETE;
#ifndef _WIN32
		/* Don't try to show prompt on systems other than Windows because it
		 * will fail anyway because of invalid hwnd. */
		if (!dsn)
			completion = SQL_DRIVER_NOPROMPT;
#endif
		i = SQLDriverConnect(h->hc, winfo.hw, UC dsname, SQL_NTS, cs,
				sizeof(cs), &cs_len, completion);
		if ((i == SQL_SUCCESS
		 ||  i == SQL_SUCCESS_WITH_INFO) && !(*dsp = conn_ds(h, cs))) {
			sprintf(TI->Errmsg = (char*)TM(strlen(CC cs) + 64),
			 "Bug? conn_ds did not find DRIVER= or DSN= in\n\t\"%s\".",
				cs);
			goto eret;
			}
		}
	else {
		if (!(s = get_ext(dsname))) {
			wantretry = 0;
			i = SQLConnect(h->hc, UC dsname, SQL_NTS, 0, 0, 0, 0);
			if (i == SQL_SUCCESS || i == SQL_SUCCESS_WITH_INFO) {
				*dsp = ds_alloc(h, dsname);
				cs[0] = 0;
				goto connected;
				}
			if (verbose)
				prc(h, "SQLConnect", i);
			sprintf(TI->Errmsg = (char*)TM(strlen(dsname)+48),
				"Could not connect with data source \"%s\".", dsname);
			goto eret;
			}
		if (match("dsn", UC s, UC (ee = s+strlen(s)))) {
			wantretry = 0;
			if (parse_dsn_file(h, &dsname))
				goto try_dsname;
			goto eret;
			}
		if (!(*dsp = ext_to_drv(UC s, UC ee, h))) {
			bad_ending("\"%s\" does not end in any of", dsname, h);
			goto eret;
			}
		e1 = s;
		dsname = fully_qualify(dsname, dsbuf, sizeof(dsbuf));

		/* Kludge around bad behavior of Access when *.mdb does not exist. */

		if (!sqlp && match("mdb", UC e1, UC ee) && not_found(ae, TI, dsname, 0))
			goto eret;
		ds = *dsp;

		/* Kludge around bad behavior of MS Excel when creating new */
		/* tables:  must use DSN=... rather than DRIVER=...  Prior  */
		/* to 20010816, we used "if (driver = get_dsn(h,ds))" on    */
		/* the theory that for consistency, it would be best to use */
		/* get_dsn whenever possible, but this can lead the wrong   */
		/* DBQ= assignment by SQLDriverConnect.			    */

		if (!sqlp && (driver = get_dsn(h,ds)))
			dreq = "DSN=";
		else {
			driver = (char*)ds->driver;
			dreq = "DRIVER=";
			}
 retry:
		s = s1 = TM(strlen(driver) + strlen(dsname) + 16);
		s = copy(s, dreq);
		s = copy(s, driver);
		*s++ = ';';
		s = copy(s, "DBQ=");
		s = copy(s, dsname);
		if (verbose) {
			printf("Calling SQLDriverConnect(\"%s\")\n", s1);
			fflush(stdout);
			}
		i = SQLDriverConnect(h->hc, winfo.hw, UC s1, SQL_NTS, cs, sizeof(cs),
			&cs_len, SQL_DRIVER_COMPLETE);
		if (i != SQL_SUCCESS && i != SQL_SUCCESS_WITH_INFO) {
			if (verbose)
				printf("Failed to connect with DRIVER=\"%s\"\n", driver);
			if (wantretry)
			    while((ds = ds->next))
				if (match(CC ds->ext, UC e1, UC ee)) {
					*dsp = ds;
					driver = (char*)ds->driver;
					goto retry;
					}
			}
		}
	if (i != SQL_SUCCESS && i != SQL_SUCCESS_WITH_INFO) {
		if (verbose)
			prc(h, "SQLDriverConnect", i);
		sprintf(TI->Errmsg = (char*)TM(strlen(dsname0)+32),
			"Could not connect to \"%s\".", dsname0);
		goto eret;
		}
 connected:
	{
		/* Quote table and column names. */
		SQLCHAR quote[2] = "\"";
		SQLSMALLINT length = 0;
		int num_cols = TI->arity + TI->ncols;
		int i = 0;
		prc(h, "SQLGetInfo", SQLGetInfo(h->hc, SQL_IDENTIFIER_QUOTE_CHAR,
				quote, sizeof(quote) / sizeof(*quote), &length));
		h->quote = quote[0];
		tname = quote_name(h, tname);
		h->quoted_colnames = TM(num_cols);
		for (i = 0; i < num_cols; ++i)
			h->quoted_colnames[i] = quote_name(h, TI->colnames[i]);
	}
	if (prc(h, "SQLAllocStmt", SQLAllocStmt(h->hc,&h->hs)))
		goto unexpected;
	if (verbose && cs[0])
		printf("Connection string: \"%s\"\n", cs);
	if (!*dsp) { /*DEBUG*/
		TI->Errmsg = "Botch in Connect(): *dsp = 0";
		goto eret;
		}
	return tname;
	}

 static void
cleanup(HInfo *h)
{
	freestmt(&h->hs);
	if (h->hc != SQL_NULL_HDBC) {
		SQLDisconnect(h->hc);
		SQLFreeConnect(h->hc);
		}
	if (h->env != SQL_NULL_HENV)
		SQLFreeEnv(h->env);
	}

 static char *
getname(HInfo *h, int *dbq)
{
	AmplExports *ae;
	TableInfo *TI;
	char *dsn, *s, *s0, *s1;
	size_t L;

	TI = h->TI;
	s = TI->strings[1];	/* DSN=...;... */
	*dbq = 0;
	if (match("DSN=", UC s, UC s+4))
		s0 = s + 4;
	else if (match("DRIVER=", UC s, UC s+7))
		s0 = s + 7;
	else {
		if (get_ext(s))
			*dbq = 1;
		return s;
		}

	/* if DBQ=filename, return dsn = filename; */
	/* else return dsname, where DSN=dsname. */

	ae = h->AE;
	for(s = s0; *s != ';' && *s; s++);
	if (*(s1 = s)) for(;;) {
		while(*++s != 'D')
			if (!*s)
				goto use_s0;
		if (s[1] != 'B' || s[2] != 'Q' || s[3] != '=')
			continue;
		s += 4;
		if (s[0] == 0 || s[0] == ';')
			break;
		*dbq = 1;
		s0 = s;
		while(*++s && *s != ';');
		s1 = s;
		break;
		}
 use_s0:
	L = s1 - s0;
	memcpy(dsn = (char*)TM(L + 1), s0, L);
	dsn[L] = 0;
	return dsn;
	}

 static char *
sql_type(HInfo *h, int odbctype, unsigned int *tprec)
{
	AmplExports *ae;
	SQLLEN len;
	HSTMT hs = h->hs;
	char nbuf[256], *rv;
	int i;
	unsigned long L;
	static char *slast, *snext;

	if (prc(h, "SQLGetTypeInfo", SQLGetTypeInfo(hs, (SQLSMALLINT)odbctype)))
		return 0;
	rv = 0;
	if (prc(h, "SQLBindCol_1",
			SQLBindCol(hs, (UWORD)1, SQL_C_CHAR, nbuf, sizeof(nbuf), &len))
	 || prc(h, "SQLBindCol_2",
			SQLBindCol(hs, (UWORD)3, SQL_C_LONG, &L, sizeof(L), &len))
	 || (i = SQLFetch(hs)) == SQL_NO_DATA
	 || prc(h, "SQLFetch in sql_types", i))
		goto done;
	if (slast - snext < len + 1) {
		i = (int)len + 1;
		if (i < 128)
			i = 128;
		snext = Malloc(h->AE, i);
		slast = snext + i;
		ae = h->AE;
		at_exit(mfree, snext);
		}
	rv = snext;
	for(i = 0; i < sizeof(nbuf) && nbuf[i]; i++)
		*snext++ = nbuf[i];
	*snext++ = 0;
	if (tprec)
		*tprec = L;
 done:
	prc(h, "SQLFreeStmt", SQLFreeStmt(hs, SQL_CLOSE));
	prc(h, "SQLFreeStmt", SQLFreeStmt(hs, SQL_UNBIND));
	return rv;
	}

 typedef struct
Sql_Type {
	char *sqlname;
	int sqltype;
	} Sql_Type;

 static int
get_types(DRV_desc *ds, HInfo *h)
{
	AmplExports *ae = h->AE;
	Sql_Type *st;
	static Sql_Type numtypes[] = {
		{ "SQL_DOUBLE", SQL_DOUBLE },
		{ "SQL_FLOAT", SQL_FLOAT },
		{ "SQL_REAL", SQL_REAL },
		{ "SQL_NUMERIC", SQL_NUMERIC },
		{0,0} };

	for(st = numtypes;; st++) {
		if (!st->sqlname) {
			h->TI->Errmsg = "Could not find SQL synonym for SQL_DOUBLE.";
			return 1;
			}
		if ((ds->ntype = sql_type(h, st->sqltype, 0)))
			break;
		if (h->verbose)
			printf("%s is not supported\n", st->sqlname);
		}
	if (!(ds->ttype = sql_type(h, SQL_TIMESTAMP, &ds->tprec))) {
		h->TI->Errmsg = "Could not find SQL synonym for SQL_TIMESTAMP.";
		return 1;
		}
	return 0;
	}

 static void
askusing(AmplExports *ae, TableInfo *TI, char *what, char *tname)
{
	char *s = TI->strings[1];
	sprintf(TI->Errmsg = (char*)TM(strlen(s) + strlen(tname) + 64),
		"%s %s failed.\n\tIs another application using \"%s\"?",
		what, tname, s);
	}

 static int
mustquote(char *s, int *slen)
{
	/* return 2 if s looks like a number, 1 if s has */
	/* white space at either end, else 0 */

	char *s0 = s, *s1;
	int dig, dot, e;

	if (*s <= ' ') {
		if (slen) {
			while(*s)
				s++;
			*slen = (int)(s - s0);
			}
		return 1;
		}
	switch(*s) {
	 case '-':
	 case '+':
		s++;
	 }
	s1 = s;
	dig = dot = e = 0;
	for(;;) switch(*s++) {
		case 0:
			goto break2;
		case '0': case '1': case '2': case '3': case '4':
		case '5': case '6': case '7': case '8': case '9':
			dig++;
			break;
		case '.':
			if (e || dot++)
				goto end_check;
			break;
		case 'e':
		case 'E':
			if (!dig || e++)
				goto end_check;
			dig = 0;
			switch(*s) {
			 case '-':
			 case '+':
				s++;
			 }
			break;
		default:
 end_check:
			if (s != s1+1)
				s1 = 0;
			while(*s)
				s++;
			if (slen)
				*slen = (int)(s - s0);
			if (s1 && (match("infinity", UC s1, UC s)
				   || match("nan", UC s1, UC s)))
				return 2;
			return s[-1] <= ' ';
		}
 break2:
	if (slen)
		*slen = (int)(s - s0);
	return dig > 0 ? 2 : 0;
	}

 static int
Write_odbc(AmplExports *ae, TableInfo *TI)
{
	DbCol *db, *db0;
	DRV_desc *ds;
	HSTMT hs;
	HInfo h;
	TIMESTAMP_STRUCT *ts, *ts0, *ts1, **tsp, ***tsq;
	UWORD u;
	char *Missing;
	char buf[32], **cn, *ct, *dt, *it, *s, **sb, **sp, **spe, *tname;
	double *rb;
	int deltry, i, i1, j, k, nc, nodrop, ntlen, nts, rc;
	int *slen, *sw;
	long ir, nr;
	size_t L, sblen, tnlen;
#ifdef NO_Adjust_ampl_odbc
#define p(x) x
#define pi(x) x
#else
	int *p_, *pi_;
#define p(x) p_[x]
#define pi(x) pi_[x]
#endif

	if ((i = ODBC_check(ae, TI, &h)))
		return i;
	rc = DB_Error;
	nc = TI->arity + TI->ncols;
	cn = TI->colnames;

	/* check validity of column names */

	colname_adjust(&h, TI);
	for(i = 0; i++ < nc; cn++) {
		s = *cn;
		if (!*s) {
			sprintf(TI->Errmsg = (char*)TM(strlen(*cn) + 60),
				"Column %d's name is the empty string.", i);
			return rc;
			}
		}

	if (!(tname = Connect(&h, &ds, &rc, 0))
	 || (!ds->ntype && get_types(ds, &h))) {
		cleanup(&h);
		return rc;
		}
	sw = 0;
	tsq = 0;
#ifndef NO_Adjust_ampl_odbc
	if ((nodrop = Adjust_ampl_odbc(&h, tname, &tsq, &p_, &pi_, &deltry, &sw))) {
		if (nodrop == 1) {
			rc = DB_Error;
			goto done;
			}
		}
	nc = TI->arity + TI->ncols;
#endif
	Missing = TI->Missing;
	hs = h.hs;
	ntlen = strlen(ds->ntype);
	tnlen = strlen(tname);
	L = tnlen + 30;	/*INSERT INTO name VALUES ()*/
	slen = (int*)TM(nc*sizeof(int));
	db0 = TI->cols;
	sblen = 0;
	cn = TI->colnames;
	nts = 0;
	for(i1 = 0; i1 < nc; i1++) {
		i = p(i1);
		db = db0 + i;
		k = 0;
		if (tsq && tsq[i])
			nts++;
		else if ((sp = db->sval)) {
			if (sw)
				k = sw[i];
			if (db->dval && k < 25)
				k = 25;
			for(spe = sp + TI->nrows; sp < spe; )
				if ((s = *sp++)) {
					if (s != Missing) {
						if (!h.nsmix)
							j = strlen(s);
						else if (mustquote(s,&j))
							j += 2;
						if (k < j)
							k = j;
						}
					}
			L += sprintf(buf, "VARCHAR(%d)", ++k);
			sblen += k;
			}
		else
			L += ntlen;
		slen[i] = k;
		L += strlen(*cn++) + 5;
		}
	if (nts)
		L += nts*strlen(ds->ttype);
	dt = (char*)TM(2*L+tnlen+16);
	ct = dt+tnlen+16;
	it = ct + L;
	rb = (double*)TM(nc*(sizeof(char*) + sizeof(double))
			+ nts*sizeof(TIMESTAMP_STRUCT) + sblen);
	ts = ts0 = (TIMESTAMP_STRUCT*)(rb + nc);
	sb = (char**)(ts + nts);
	cn = TI->colnames;
	if (nodrop) {
		j = sprintf(it = ct, "INSERT INTO %s (\"%s\"" /*)*/, tname, cn[0]);
		for(i = 1; i < nc; ++i)
			j += sprintf(it+j, ", \"%s\"", cn[i]);
		j += /*(*/ sprintf(it+j, ") VALUES (?"); /*)*/
		goto finish_insprep;
		}
	j = sprintf(ct, "CREATE TABLE %s ("/*)*/, tname);
	for(i1 = 0; i1 < nc; i1++) {
		i = p(i1);
		db = db0 + i;
		j += sprintf(ct+j, "%s%s ", i1 ? ", " : "", h.quoted_colnames[i]);
		if (tsq && tsq[i])
			j += sprintf(ct+j, "%s", ds->ttype);
		else if (db->sval)
			j += sprintf(ct+j, "VARCHAR(%d)", slen[i]);
		else
			j += sprintf(ct+j, "%s", ds->ntype);
		}
	ct[j++] = /*(*/')';
	ct[j] = 0;
	if (deltry) {
		sprintf(dt, "DELETE FROM %s", tname);
		if ((i = SQLExecDirect(hs, UC dt, SQL_NTS)) == SQL_SUCCESS
		  || i == SQL_SUCCESS_WITH_INFO) {
			if (h.verbose)
				printf("\"%s\" succeeded.\n", dt);
			goto prepare_insert;
			}
		if (h.verbose)
			printf("\"%s\" failed.\n", dt);
		goto droptry;
		}
	if ((i = SQLExecDirect(hs, UC ct, SQL_NTS)) != SQL_SUCCESS
	 && i != SQL_SUCCESS_WITH_INFO) {
		if (h.verbose || nodrop)
			prc(&h, ct, i);
		if (nodrop)
			goto done;
 droptry:
		sprintf(dt, "DROP TABLE %s", tname);
		if ((i = SQLExecDirect(hs, UC dt, SQL_NTS)) == SQL_SUCCESS
		  || i == SQL_SUCCESS_WITH_INFO) {
			if (h.verbose)
				printf("\"%s\" succeeded.\n", dt);
			if ((i = SQLExecDirect(hs, UC ct, SQL_NTS)) != SQL_SUCCESS
			 && i != SQL_SUCCESS_WITH_INFO) {
				if (h.verbose)
					prc(&h, ct, i);
				askusing(ae, TI, "CREATE TABLE", tname);
				goto done;
				}
			}
		else  {
			if (h.verbose)
				prc(&h, dt, i);
			askusing(ae, TI, "DROP TABLE", tname);
			goto done;
			}
		}
 prepare_insert:
	j = sprintf(it, "INSERT INTO %s VALUES (?"/*)*/, tname);
 finish_insprep:
	for(i = 1; i < nc; i++) {
		strcpy(it+j, ", ?");
		j += 3;
		}
	strcpy(it+j,/*(*/")");
	if (prc(&h, "Prepare(INSERT)", SQLPrepare(hs, UC it, SQL_NTS))) {
 failed:
		TI->Errmsg = "Unexpected ODBC failure";
		goto done;
		}
	s = (char*)(sb + nc);
	db = TI->cols;
	for(i = 0; i < nc; i++, db++) {
		u = pi(i) + 1;
		if (tsq && tsq[i])
			SQLBindParameter(hs, u, SQL_PARAM_INPUT, SQL_C_TIMESTAMP,
					SQL_TIMESTAMP, ds->tprec, 0, ts++, 0, NULL);
		else if (db->sval) {
			SQLBindParameter(hs, u, SQL_PARAM_INPUT, SQL_C_CHAR,
					SQL_VARCHAR, slen[i], 0, sb[i] = s, 0, NULL);
			s += slen[i];
			}
		else
			SQLBindParameter(hs, u, SQL_PARAM_INPUT, SQL_C_DOUBLE,
					SQL_DOUBLE, 0, 0, rb + i, 0, NULL);
		}
	nr = TI->nrows;
	for(ir = 0; ir < nr; ir++) {
		db = TI->cols;
		ts1 = ts0;
		for(i = 0; i < nc; i++, db++)
			if (tsq && (tsp = tsq[i]))
				*ts1++ = *tsp[ir];
			else if ((sp = db->sval)) {
				if ((s = sp[ir])) {
					if (s == Missing)
						*sb[i] = 0;
					else if (h.oldquotes && mustquote(s,0))
						sprintf(sb[i], "'%s'", s);
					else
						strcpy(sb[i], s);
					}
				else if (db->dval)
					sprintf(sb[i], "%.g", db->dval[ir]);
				else
					*sb[i] = 0;
				}
			else
				rb[i] = db->dval[ir];
		i = SQLExecute(hs);
		if (i != SQL_SUCCESS && i != SQL_SUCCESS_WITH_INFO) {
			if (!ir && !nodrop) {
				askusing(ae, TI, "INSERT INTO", tname);
				goto done;
				}
			prc(&h, "SQLExecute", i);
			goto failed;
			}
		}
	if (prc(&h, "COMMIT CREATE INSERT", SQLTransact(h.env, h.hc, SQL_COMMIT)))
		TI->Errmsg = "Write commit failure";
	else
		rc = DB_Done;
 done:
	cleanup(&h);
	return rc;
	}
#undef p
#undef pi

 typedef struct
DBColinfo {
	char *name;
	char *val;
	char *val0;
	SQLLEN type;
	SQLLEN prec;
	SQLLEN len;
	int mytype;
	int myoffset;
	} DBColinfo;

 static int
dbcompar(const void *a, const void *b, void *v)
{
	return strcmp(((DBColinfo*)v)[*(int*)a].name,
		      ((DBColinfo*)v)[*(int*)b].name);
	}

 static int
icompar(const void *a, const void *b, void *v)
{
	return ((int*)v)[*(int*)a] - ((int*)v)[*(int*)b];
	}

 static int
permute(int nf, int nt, int *p, char **namf, int *zf, DBColinfo *dbc, int *zt)
{
	int i, j, k, zi, zj;
	char *s;

	for(i = j = 0; i < nf; i++, j++) {
		s = namf[zi = zf[i]];
		for(;;j++) {
			if (j >= nt)
 notfound:			return zi+1;
			if (!(k = strcmp(dbc[zj = zt[j]].name, s)))
				break;
			if (k > 0)
				goto notfound;
			}
		p[zi] = zj;
		}
	return 0;
	}

 static int
needprec(HInfo *h, DBColinfo *dbc, int col_index)
{
	AmplExports *ae;
	HSTMT hs;
	SDWORD t;
	SQLLEN len;
	TableInfo *TI;
	UnknownType *ut, **utp;
	char sbuf[256];
	int i;

	dbc->prec = 0;
	switch(dbc->type) {
	 case SQL_NUMERIC:
	 case SQL_DECIMAL:
	 case SQL_INTEGER:
	 case SQL_SMALLINT:
	 case SQL_FLOAT:
	 case SQL_REAL:
	 case SQL_DOUBLE:
		return dbc->mytype = 0;
	 case SQL_CHAR:
	 case SQL_VARCHAR:
		return dbc->mytype = 1;
	 case SQL_DATE:
	 case SQL_TIME:
	 case SQL_TIMESTAMP:
		dbc->mytype = 2;
		return 0;
	 }
	dbc->mytype = 3;
	t = dbc->type;
	for(utp = &h->ut; (ut = *utp); utp = &ut->next)
		if (ut->type == t)
			goto found;
	if (prc(h, "SQLAllocStmt", SQLAllocStmt(h->hc,&hs)))
		return 0;
	sbuf[0] = 0;
	i = prc(h, "SQLGetTypeInfo in needprec", SQLGetTypeInfo(hs, (SQLSMALLINT)t))
	 || prc(h, "SQLBindCol_3",
			SQLBindCol(hs, (UWORD)5, SQL_C_CHAR, sbuf, sizeof(sbuf), &len));
	if (!i) {
		i = SQLFetch(hs);
		/* SQLFetch after SQLGetTypeInfo returns no data for some MySQL
		   types. */
		if (i != SQL_SUCCESS && i != SQL_NO_DATA)
			prc(h, "SQLFetch in needprec()", i);
	}
	prc(h, "SQLFreeStmt", SQLFreeStmt(hs, SQL_DROP));
	if (i == SQL_NO_DATA) {
		SQLSMALLINT length = 0;
		/* Fall back to checking column type name if there is no type info. */
		if (prc(hs, "SQLColAttributes(SQL_COLUMN_TYPE_NAME)",
				SQLColAttributes(h->hs, col_index, SQL_COLUMN_TYPE_NAME, sbuf,
						(SWORD)sizeof(sbuf), &length, 0))) {
			return 0;
		}
		if (strcmp(sbuf, "varchar") == 0)
			i = 1;
	} else if (!i) {
		switch(sbuf[0]) {
			case 0:
				i = 0;
				break;
			case '#':
				i = 2;
				break;
			default:
				i = 1;
			}
	} else
		return 0;
	ae = h->AE;
	TI = h->TI;
	ut = *utp = (UnknownType*)TM(sizeof(UnknownType));
	ut->next = 0;
	ut->type = t;
	ut->mytype = i;
 found:
	return (dbc->mytype = ut->mytype) == 1;
	}

 static char*
Mem(HInfo *h, size_t L)
{
	AmplExports *ae;
	TableInfo *TI;
	char *rv;
	size_t L1;

	if (h->se - h->s < L) {
		ae = h->AE;
		TI = h->TI;
		L1 = ((L+32) >> 12) + 1;
		L1 = (L1 << 12) - 32;
		h->se = (h->s = (char*)TM(L1)) + L1;
		}
	rv = h->s;
	h->s += L;
	return rv;
	}

 static char *
scrunch(HInfo *h, char *s, real *d, int mix)
{
	AmplExports *ae;
	char *t;

	for(t = s; *t; t++);
	while(t > s && *--t == ' ');
	t[1] = 0;
	if (t == s && *s == ' ')
		*s = 0;
	switch(*s) {
	  case '\'':
		if (t > s && *t == '\'') {
			*t = 0;
			s++;
			}
		break;
	  case 0:
		s = h->Missing;
		break;
	  default:
		if (mix && mustquote(s,0) == 2) {
			ae = h->AE;
			*d = strtod(s,0);
			s = 0;
			}
	  }
	return s;
	}

#ifdef SELECT_JUST_DESIRED_COLUMNS
 static char *
select_stmt(AmplExports *ae, TableInfo *TI, char *tname, size_t L)
{
	char **cn, **cne, *rv, *s;
	int nc;

	/* input L == strlen(tname) */

	nc = TI->arity + TI->ncols;
	cn = TI->colnames;
	cne = cn + nc;
	L += 4*nc + 20;
	do L += strlen(*cn++);
		while(cn < cne);
	rv = s = (char*)TM(L);
	cn = TI->colnames;
	s += sprintf(s, "SELECT \"%s\"", *cn);
	while(++cn < cne)
		s += sprintf(s, ", \"%s\"", *cn);
	sprintf(s, " FROM %s", tname);
	return rv;
	}
#endif /* SELECT_JUST_DESIRED_COLUMNS */

 static SWORD sqlc[3] = { SQL_C_DOUBLE, SQL_C_CHAR, SQL_C_TIMESTAMP };

 static TIMESTAMP_STRUCT Missing_time = { 1900,1,1,12,34,56,0 }, No_time;
	/* The Missing_time value is a kludge for Excel, which does not */
	/* correctly handle years before 1900.  For Access, we could use */
	/* { 0, 3, 3,  0,0,0, 0 }. */

 static int
Read_odbc(AmplExports *ae, TableInfo *TI)
{
	DBColinfo *dbc0, *dbc;
	DRV_desc *ds;
	DbCol *db;
	HInfo h;
	HSTMT hs;
	PTR ptr;
	SWORD len, ncols;
	TIMESTAMP_STRUCT *td, *ts;
	UWORD u;
	char **cd, *dsn, nbuf[512], *s, *sbuf, *tname;
	double *dd, t;
	int a, *ct, dbq, i, j, k, mix, nc, nf, nk[4], nt, *p, *z, *zt;
	size_t L;
	typedef struct NameType { const char *name; int type; } NameType;
	NameType *ntp;
	static NameType NT[] = {
		{ "SQL_CHAR", SQL_CHAR },
		{ "SQL_NUMERIC", SQL_NUMERIC },
		{ "SQL_DECIMAL", SQL_DECIMAL },
		{ "SQL_INTEGER", SQL_INTEGER },
		{ "SQL_SMALLINT", SQL_SMALLINT },
		{ "SQL_FLOAT", SQL_FLOAT },
		{ "SQL_REAL", SQL_REAL },
		{ "SQL_DOUBLE", SQL_DOUBLE },
		{ "SQL_DATE", SQL_DATE },
		{ "SQL_TIME", SQL_TIME },
		{ "SQL_TIMESTAMP", SQL_TIMESTAMP },
		{ "SQL_VARCHAR", SQL_VARCHAR },
		{ 0, 0}
		};
	static int first = 1;

	if ((i = ODBC_check(ae, TI, &h)))
		return i;
	if (!(dsn = getname(&h, &dbq)))
		return DB_Error;
	if (dbq && not_found(ae, TI, dsn, 0))
		return DB_Error;

	sbuf = 0;
	colname_adjust(&h, TI);
	if (!(tname = Connect(&h, &ds, &i, &sbuf))) {
		cleanup(&h);
		return i;
		}

	L = strlen(tname);
	if (!(s = sbuf)) {
#ifdef SELECT_JUST_DESIRED_COLUMNS
		/* This may be more efficient, but causes case to be ignored */
		/* with at least some ODBC drivers, threatening confusion. */
		sbuf = select_stmt(ae, TI, tname, L);
		if ((i = SQLPrepare(hs = h.hs, sbuf, SQL_NTS)) == SQL_SUCCESS
		 ||  i == SQL_SUCCESS_WITH_INFO)
			goto select_worked;
#endif
		sprintf(sbuf = (char*)TM(L + 32), "SELECT ALL * FROM %s", tname);
		}
	if (prc(&h, "SQLPrepare", SQLPrepare(hs = h.hs, UC sbuf, SQL_NTS))) {
		if (s) {
			sprintf(TI->Errmsg = (char*)TM(L + 32),
				"%s did not work.", tname);
			goto bailout;
			}
 badret:
		TI->Errmsg = "surprise ODBC failure";
 bailout:
		cleanup(&h);
		return DB_Error;
		}
 /*select_worked:*/
	if (prcnr(&h, "SQLNumResultCols", SQLNumResultCols(hs, &ncols))) {
		L += strlen(dsn);
		sprintf(TI->Errmsg = (char*)TM(L + 32),
			"Table %s does not appear in \"%s\".", tname, dsn);
		goto bailout;
		}
	dbc0 = dbc = (DBColinfo *)TM(ncols*sizeof(DBColinfo));
	a = TI->arity;	/* number of indexing columns */
	nc = TI->ncols;	/* number of data columns desired */
	nf = a + nc;	/* total number of columns of interest */
	i = nt = ncols;	/* columns in the database table */
	if (i < nf)
		i = nf;	/* we'll complain about a specific missing column */
	p = (int*)TM(3*sizeof(int)*i);
	z = p + i;
	zt = z + i;
	strsort(ae, nf, TI->colnames, compar, z);
	for(i = 1; i <= nt; i++, dbc++) {
		if (prc(&h, "SQLColAttributes(NAME)",
 				SQLColAttributes(hs, (UWORD)i, SQL_COLUMN_NAME, nbuf,
						(SWORD)sizeof(nbuf), &len, 0)))
 			goto badret;
		memcpy(dbc->name = (char*)TM(len+1), nbuf, len);
		dbc->name[len] = 0;
		}
	strsort(ae, nt, dbc0, dbcompar, zt);
	if ((i = permute(nf, nt, p, TI->colnames, z, dbc0, zt))) {
		L += strlen(s = TI->colnames[i-1]) + strlen(dsn) + 200 + 6*nt;
		for(i = nt; i > 9; i /= 10)
			L += nt;
		for(i = 0, dbc = dbc0; i < nt; i++, dbc++)
			L += strlen(dbc->name);
		i = sprintf(TI->Errmsg = (char*)TM(L),
			"Column \"%s\" does not appear in table %s of \"%s\",\n",
			s, tname, dsn);
		s = TI->Errmsg + i;
		s += sprintf(s, "which has %d column%s:\n\n\tCol\tName\n\n",
			nt, nt != 1 ? "s" : "");
		for(i = 0; i < nt; i++) {
			j = zt[i];
			s += sprintf(s, "\t%d\t\"%s\"\n", j+1, dbc0[j].name);
			}
		goto bailout;
		}
	memset(nk, 0, sizeof(nk));
	if (h.sqldb) {
		if (first) {
			first = 0;
			printf("\nPossible SQL column types:\n\n");
			for(ntp = NT; ntp->name; ++ntp)
				printf("\t%s\t= %d\n", ntp->name, ntp->type);
			}
		printf("\nTable %s:\n\n", TI->tname);
		}
	for(i = 0; i < nf; i++) {
		j = p[i];
		dbc = dbc0 + j;
		u = j + 1;
		if (prc(&h, "SQLColAttributes(TYPE)",
		 SQLColAttributes(hs, u, SQL_COLUMN_TYPE, 0, 0, 0, &dbc->type)))
			goto badret;
		if (needprec(&h, dbc, u)) {
			if (prc(&h, "SQLColAttributes(DISPLAY_SIZE)",
			 SQLColAttributes(hs, u, SQL_COLUMN_DISPLAY_SIZE,
				0,0,0, &dbc->prec)))
				goto badret;
			if (dbc->prec <= 0 || ++dbc->prec > h.maxlen)
				dbc->prec = h.maxlen;
			}
		if (h.sqldb)
			printf("Column %d = \"%s\" has SQL type %d, mytype %d\n",
				j, dbc->name, dbc->type, dbc->mytype);
		nk[dbc->mytype]++;
		if (dbc->mytype == 3) {
			nbuf[0] = 0;
			if (prc(&h, "SQLColAttributes(SQL_COLUMN_TYPE_NAME)",
			 SQLColAttributes(hs, u, SQL_COLUMN_TYPE_NAME, nbuf,
						(SWORD)sizeof(nbuf), &len, 0)))
				goto badret;
			L += strlen(s = TI->colnames[i]) + strlen(dsn);
			sprintf(TI->Errmsg = (char*)TM(L + len + 90),
		"Column \"%s\" in table %s of \"%s\" has unknown type %ld\nnamed \"%s\".",
				s, tname, dsn, (long)dbc->type, nbuf);
			goto bailout;
			}
		}
	if (h.sqldb)
		printf("\n");
	dd = 0; /* shut up erroneous warning */
	cd = 0; /* ditto */
	td = 0; /* ditto */
	if ((i = nk[0]))
		dd = (double*)TM(i*sizeof(double));
	if ((i = nk[1]))
		cd = (char**)TM(i*sizeof(char*));
	if ((i = nk[2])) {
		td = (TIMESTAMP_STRUCT*)TM(i*sizeof(TIMESTAMP_STRUCT));
		}
	memset(nk, 0, sizeof(nk));
	for(i = 0; i < nf; i++) {
		j = p[i];
		dbc = dbc0 + j;
		u = j + 1;
		dbc->myoffset = j = nk[k = dbc->mytype]++;
		switch(k) {
		 case 0:
			ptr = (PTR)&dd[j];
			break;
		 case 1:
			ptr = (PTR)(s = cd[j] = (char*)TM(dbc->prec));
			*s = 0;	/* bypass MS Excel bug with empty cells */
			break;
		 case 2: ptr = (PTR)&td[j];
			td[j] = Missing_time;
			break;
		 default: continue;
		 }
		if (prc(&h, "SQLBindCol_3",
				SQLBindCol(hs, u, sqlc[k], ptr, dbc->prec, &dbc->len)))
			goto badret;
		}
	if (h.oldquotes)
		mix = 3;
	else if ((mix = h.nsmix) == 2) {
		if ((ct = h.coltypes)) {
			for(i = mix = 0; i < nf; ++i) {
				if (ct[i] == 3) {
					dbc = dbc0 + p[i];
					if (dbc->mytype == 1)
						dbc->mytype = 3;
					}
				}
			}
		else
			mix = 1;
		}
	if (prc(&h, "SQLExecute", SQLExecute(hs)))
		goto badret;
	db = TI->cols;
	while((i = SQLFetch(hs)) == SQL_SUCCESS || i == SQL_SUCCESS_WITH_INFO) {
		for(i = 0; i < nf; i++, db++) {
			dbc = dbc0 + p[i];
			switch(dbc->mytype) {
			 case 0:
				db->sval[0] = 0;
				db->dval[0] = dd[dbc->myoffset];
				break;
			 case 1:
				db->sval[0] = scrunch(&h, cd[dbc->myoffset], db->dval, mix);
				break;
			 case 2:
				ts = td + dbc->myoffset;
				if (ts->year < 0
				 || !memcmp(ts, &No_time, sizeof(No_time))
				 || !memcmp(ts, &Missing_time, sizeof(No_time))) {
					db->sval[0] = TI->Missing;
					break;
					}
				db->sval[0] = 0;
				t = 100.*ts->year + ts->month;
				t = 100.*t + ts->day;
				t = 100.*t + ts->hour;
				t = 100.*t + ts->minute;
				t = 100.*t + ts->second;
				/*if (ts->fraction)
					t += ts->fraction / 4294967296.;*/ /* 2^32 */
				db->dval[0] = t;
				*ts = Missing_time;
				break;
			  case 3:
				db->sval[0] = scrunch(&h, cd[dbc->myoffset], db->dval, 1);
			 }
			}
		db = TI->cols;
		if ((*TI->AddRows)(TI, db, 1)) {
			cleanup(&h);
			return DB_Error;
			}
		if (a <= 0)
			break;
		for(i = nk[1]; i--;)
			*cd[i] = 0;	/* to bypass MS Excel empty-cell bug */
		}
	cleanup(&h);
	return DB_Done;
	}

 void
funcadd(AmplExports *ae)
{
	static char tname[] = "odbc\n"
	"AMPL ODBC handler: expected 2-8 strings before \":[...]\":\n"
	"  'ODBC', connection_spec ['ext_name'] [option [option...]]\n"
	"Connection_spec gives a connection to the external table.  If the table's\n"
	"external name differs from the AMPL table name, the external name must be\n"
	"given in place of 'ext_name'.  For IN tables, 'SQL=sqlstmt' can appear in\n"
	"place of 'ext_name', where sqlstmt is a SQL statement, such as a SELECT\n"
	"statement.  Possible options, explained below:\n\n"
	"\t'maxlen=nnn'\n"
	"\t'nsmix=...'\n"
	"\t'time=...'\n"
	"\t'verbose' or 'verbose=n'\twith 0 <= n <= 3\n"
	"\t'write=append' or 'write=drop'\n\n"
	"Use 'maxlen=nnn' to limit character strings to nnn characters (discarding\n"
	"any excess characters).\n\n"
	"With 'nsmix=*', columns of string data are treated as containing both\n"
	"strings and numbers, with strings that look like decimal fixed- or floating-\n"
	"point numbers treated as numbers (the default).  With 'nsmix=-', columns of\n"
	"string data are kept as string data.  To give the effect of 'nsmix=*' to\n"
	"some columns and of 'nsmix=-' to others, use syntax\n"
	"\tindex ~ 'Strcol:colname'\n"
	"in the AMPL table declaration for columns that should be treated as with\n"
	"the global 'nsmix=-' string.  For columns to be treated in the default way,\n"
	"optionally you can use syntax \"index ~ 'Mixed:colname'\" in the table\n"
	"declaration.  The \"Strcol:\" and \"Mixed:\" portion will be stripped from\n"
	"the column name presented to the database.   Alternatively, use 'nsmix=!'\n"
	"to use the old inconsistent practice of adding quotes to strings that look\n"
	"like numbers when reading but not removing the quotes when writing.\n\n"
	"The ... in 'time=...' is a comma-separated list of names of columns that\n"
	"should be regard as time data when writing tables.  A better alternative\n"
	"is to use syntax \"index ~ 'Time:colname'\" in the table declaration.\n"
	"As with \"Strcol:\" and \"Mixed:\", the \"Time:\" portion will be stripped\n"
	"from the column name shown to the database.\n\n"
	"For OUT and INOUT tables, 'write=...' specifies how \"write table\" should work:\n"
	"\t'write=drop' ==> drop and completely rewrite an existing table (default)\n"
	/*"\t'write=update' ==> use SQL UPDATE TABLE to update existing rows (default)\n"*/
	"\t'write=append' ==> assume \"write table\" is only appending new rows\n\n"
	"For 'verbose=n', n is the sum of\n"
	"\t\t1 ==> show connection strings and\n"
	"\t\t2 ==> show column types.\nPlain 'verbose' is treated as 'verbose=1'.\n\n"
	"Alternatives for connection_spec:\n"
	"\t'filename.ext', where ext is a registered ODBC extension;\n"
	"\t'filename.dsn' (written by the ODBC control panel's \"File DSN\");\n"
	"\tan explicit connection string of the form 'DSN=...' or 'DRIVER=...';\n"
	"\tor an ODBC data source name (see the ODBC control panel)."
	;

	if (!ae->asl)
	add_table_handler(Read_odbc, Write_odbc, tname, 0, 0);
	}

#ifndef NO_Adjust_ampl_odbc

 static int
apermute(int a, int nf, int nt, int *p, char **namf, int *zf, DBColinfo *dbc,
	 int *zt, int *nn, int *nfn)
{
	int i, j, k, kfn, n, zi, zj;
	char *s;

	for(i = j = kfn = n = 0; i < nf; i++) {
		s = namf[zi = zf[i]];
		for(;;j++) {
			if (j >= nt) {
				if (zi < a)
					return zi + 1;
				goto ret;
				}
			if (!(k = strcmp(dbc[zj = zt[j]].name, s))) {
				j++;
				break;
				}
			if (k > 0) {
				nfn[++kfn] = zi;
				if (zi < a)
					return zi + 1;
				goto break2;
				}
			p[zj] = nf + n++;
			}
		p[zj] = zi;
 break2:	;
		}
	while(j < nt)
		p[zt[j++]] = nf + n++;
 ret:
	*nn = n;
	*nfn = kfn;
	return 0;
	}

 static TIMESTAMP_STRUCT *
timeval(HInfo *h, TIMESTAMP_STRUCT *ts, real t, long nrow)
{
	int bad = 0, i;

	i = (int)(t/1e10);
	if (i*1e10 > t)
		--i;
	ts->year = i;
	t -= i*1e10;
	ts->month = i = (int)(t/1e8);
	if (i <= 0 || i > 12)
		bad++;
	t -= i*1e8;
	ts->day = i = (int)(t/1e6);
	if (i == 0 && ts->month == 0 && ts->year == 0)
		bad = 0;
	else if (i <= 0 || i > 31)
		bad++;
	t -= i*1e6;
	ts->hour = i = (int)(t/1e4);
	if (i < 0 || i > 59)
		bad++;
	t -= i*1e4;
	ts->minute = i = (int)(t/1e2);
	if (i < 0 || i > 59)
		bad++;
	t -= i*1e2;
	ts->second = i = (int)t;
	if (i < 0 || i > 59)
		bad++;
	ts->fraction = 0 /*(t-i)*4294967296.*/;  /* nonzero values fail */
	if (bad) {
		i = h->thiscol;
		if (++h->totbadtimes == 1) {
			h->firstbadtimes = h->totbadtimes = 1;
			h->firstbadcol = i;
			h->firstbadrow = nrow;
			h->ts = ts;
			}
		else if (i == h->firstbadcol)
			h->firstbadtimes++;
		}
	return ts;
	}

 static TIMESTAMP_STRUCT **
time_convert(HInfo *h, int j)
{
	AmplExports *ae = h->AE;
	TIMESTAMP_STRUCT **tsp;
	TableInfo *TI = h->TI;
	DbCol *db = TI->cols + j;
	char *Missing = TI->Missing, **sa = db->sval;
	real *ra = db->dval;
	long m = 0, nrows = TI->nrows;

	h->thiscol = j;
	if (sa) {
		tsp = (TIMESTAMP_STRUCT**)sa;
		if (ra)
			for(; m < nrows; m++)
				tsp[m] = (sa[m] == Missing)
					? &Missing_time
					: timeval(h,(TIMESTAMP_STRUCT *)
						TM(sizeof(TIMESTAMP_STRUCT)),ra[m],m);
		else
			goto no_time;
		}
	else {
		tsp = (TIMESTAMP_STRUCT**)(*TI->ColAlloc)(TI, j, 1);
		if (ra)
			for(; m < nrows; m++)
				tsp[m] = timeval(h,(TIMESTAMP_STRUCT *)
						TM(sizeof(TIMESTAMP_STRUCT)),ra[m],m);
		else {
 no_time:
			while(m < nrows)
				tsp[m++] = &Missing_time;
			}
		}
	return tsp;
	}

 static char *
ascrunch(HInfo *h, char *s, real **dp, long j, int k)
{
	TableInfo *TI;
	char *t;
	real *d;

	if (!s) {
		if (!(d = *dp)) {
			TI = h->TI;
			d = *dp = (real*)(*TI->ColAlloc)(TI, k, 0);
			}
		d[j] = h->dd[k];
		}
	else if (s != h->Missing) {
		strcpy(t = Mem(h,strlen(s)+1), s);
		s = t;
		}
	return s;
	}

 static int
Adjust_ampl_odbc(HInfo *h, char *tname, TIMESTAMP_STRUCT ****tsqp,
		 int **pp, int **pip, int *deltry, int **swp)
{
	AmplExports *ae;
	DBColinfo *dbc0, *dbc, *dbce;
	DbCol *db, *db0, *db1, *dbe;
	HSTMT hs;
	PTR ptr;
	SWORD len, ncols;
	TIMESTAMP_STRUCT *td, *ts, **tsp, ***tsq, **tsx;
	TableInfo *TI = h->TI;
	UWORD u;
	char *Missing, buf[512], **cd, *cs, dbuf[32], *s, *seen, **sa, **sp;
	double *dd, t;
	int a = TI->arity, i, i1, j, k, kfn, mix, n, nc, *nfn, nk[4], nn;
	int nf = a + TI->ncols;
	int nnt, nt, nt0, nt1, ntimes, nts, rc, wantsv;
	int *cc, *ct = h->coltypes, *p, *p_, *pi, *sw, *zf, *zt;
	long m, maxrows, nrows, nrows0, nseen;
	real *ra, *ra0;
	size_t L, Lt;

	ae = h->AE;
	*deltry = 0;
	ntimes = h->ntimes;
	h->totbadtimes = nn = rc = 0;
	tsq = 0;
	tsx = 0;
	p_ = pi = 0;
	if (!(TI->flags & DBTI_flags_IN) && h->wrmode != wr_append)
		goto done1;
	cs = buf;
	Lt = strlen(tname);
	if ((L = Lt + 20) > sizeof(buf))
		cs = (char*)TM(L);
	sprintf(cs, "SELECT ALL * FROM %s", tname);
	if (prc(h, "SELECT ALL * in Adjust_ampl_odbc",
		SQLPrepare(hs = h->hs, UC cs, SQL_NTS)))
		goto done1;
	if ((i = SQLNumResultCols(hs, &ncols)) != SQL_SUCCESS
	 && i != SQL_SUCCESS_WITH_INFO)
		goto done;
	h->s = h->se = 0;
	rc = 1;
	nt = ncols;
	dbc0 = dbc = (DBColinfo *)TM(nt*sizeof(DBColinfo) + (nt+1)*sizeof(int));
	nfn = (int*)(dbc + nt);
	memset(nk, 0, sizeof(nk));
	for(i = 1; i <= nt; i++, dbc++) {
		if (prc(h, "SQLColAttributes(NAME)",
 				SQLColAttributes(hs, (UWORD)i, SQL_COLUMN_NAME, buf,
						(SWORD)sizeof(buf), &len, 0))) {
 badret:
			sprintf(TI->Errmsg = Mem(h,Lt+48),
				"ODBC failure reading old table %s for update.", tname);
 			goto done;
			}
		memcpy(dbc->name = Mem(h,len+1), buf, len);
		dbc->name[len] = 0;
		if (prc(h, "SQLColAttributes(TYPE)",
		 SQLColAttributes(hs, (UWORD)i, SQL_COLUMN_TYPE, 0, 0, 0, &dbc->type)))
			goto badret;
		if (needprec(h, dbc, i)) {
			if (prc(h, "SQLColAttributes(DISPLAY_SIZE)",
			 SQLColAttributes(hs, (UWORD)i, SQL_COLUMN_DISPLAY_SIZE,
				0,0,0, &dbc->prec)))
				goto badret;
			if (dbc->prec <= 0 || ++dbc->prec > h->maxlen)
				dbc->prec = h->maxlen;
			}
		nk[dbc->mytype]++;
		if (dbc->mytype == 3) {
			sprintf(TI->Errmsg = Mem(h, Lt + len + 72),
				"Column %s in table %s has unexpected type %ld.",
				dbc->name, tname, (long)dbc->type);
			goto done;
			}
		}
	dbce = dbc;
	zf = (int*)TM((nf+2*nt+nk[1])*sizeof(int));
	zt = zf + nf;
	p = zt + nt;
	cc = p + nt;
	strsort(ae, nf, TI->colnames, compar, zf);
	strsort(ae, nt, dbc0, dbcompar, zt);
	if ((i = apermute(a, nf, nt, p, TI->colnames, zf, dbc0, zt, &nn, nfn))) {
		L = Lt + strlen(cs = TI->colnames[i-1]) + 80;
		sprintf(TI->Errmsg = Mem(h,L), "%s%s\n\tcolumn \"%s\"%s",
			"Cannot update existing table ", tname, cs,
			 " is not in the existing table.");
		goto done;
		}
	if ((kfn = *nfn)) {
		if (h->wrmode != wr_drop) {
			L = 120;
			for(i = 0; i < kfn; ++i)
				L += strlen(TI->colnames[i]) + 2;
			TI->Errmsg = s = Mem(h,L);
			s += sprintf(s, "Cannot use 'write=append' because write table "
					"%s needs to add\n", tname);
			if (kfn == 1)
				s += sprintf(s, "column \"%s\" to the table.",
					TI->colnames[nfn[1]]);
			else {
				s += sprintf(s, "%d columns to the table:", kfn);
				for(i = 0; i < kfn; )
					s += sprintf(s, "\n\t\"%s\"", TI->colnames[nfn[++i]]);
				}
			goto done;
			}
		}
	else if (h->wrmode == wr_append) {
		rc = 2;
		goto done;
		}
	n = nf + nn;
	if (nn) {
		db = (DbCol*)TM(n*(sizeof(DbCol) + sizeof(char*)));
		memcpy(db, TI->cols, nf*sizeof(DbCol));
		memset(db + nf, 0, nn*sizeof(DbCol));
		TI->cols = db;
		sp = (char**)(db + n);
		memcpy(sp, TI->colnames, nf*sizeof(char*));
		TI->colnames = sp;
		for(i = 0; i < nt; i++)
			if ((j = p[i]) >= nf)
				sp[j] = dbc0[i].name;
		}
	if (!nn && nf == nt) {
		*deltry = 1;
		if (memcmp(zf,zt,nf*sizeof(int)))
			goto p_comp;
		}
	else {
 p_comp:
		pi = (int*)TM(2*sizeof(int)*n);
		p_ = pi + n;
		for(i = 0; i < n; i++)
			pi[i] = -1;
		for(i = 0; i < nt; i++)
			pi[p[i]] = i;
		for(j = nt, i = 0; i < n; i++)
			if (pi[i] < 0)
				pi[i] = j++;
		for(i = 0; i < n; i++)
			p_[pi[i]] = i;
		}
	db0 = TI->cols;
	ra = ra0 = 0;
	sa = 0;
	j = k = 0;
	for(i = 0; i < nt; i++) {
		if ((p[i] < a)) {
			if (dbc0[i].mytype == 1)
				k++;
			else
				j++;
			}
		}
	if (nk[0]) {
		dd = (double*)TM(n*sizeof(double));
		if (j)
			ra0 = dd;
		}
	if ((i = nk[1])) {
		i1 = i;
		if (k)
			i += a;
		cd = (char**)TM(i*sizeof(char*));
		if (k)
			memset(sa = cd + i1, 0, a*sizeof(char*));
		}
	if ((nnt = nts = nt1 = nk[2])) {
		td = (TIMESTAMP_STRUCT*)TM(nts*sizeof(TIMESTAMP_STRUCT));
		for(i = 0; i < nts; i++)
			td[i] = Missing_time;
		for(i = j = 0, dbc = dbc0; i < nt; i++, dbc++)
			if (dbc->mytype == 2) {
				if (p[i] < nf) {
					zf[j] = p[i];
					zt[j++] = i;
					}
				else
					zt[--nt1] = i;
				}
		}
	if (nn) {
		for(i = 0; i < nt; i++)
			if (p[i] >= nf && dbc0[i].mytype != 2)
				zt[nnt++] = i;
		}

	/* Now zt[i], i in [0,nts) are timestamp columns of the existing table,   */
	/* of which zt[i], i in [0,nt1) are in TI, */
	/* and zt[i], i in [nt1, nnt) are columns of the existing table not in TI. */

	nrows = nrows0 = TI->nrows;
	Missing = TI->Missing;
	nc = nf + nn;
	if (ntimes) {
		tsx = (TIMESTAMP_STRUCT **)TM(nc*sizeof(TIMESTAMP_STRUCT*));
		memset(tsx, 0, nc*sizeof(TIMESTAMP_STRUCT**));
		for(i = 0; i < nt; i++)
			if ((j = p[i]) < nf)
				tsx[j] = &Missing_time;

		/* Let the columns of the existing table speak for themselves */
		/* about whether they are time-stamp columns. */
		for(i = 0; i < nf; ++i) {
			if (ct[i] == 2 && tsx[i])
				--j;
			}
		ntimes = j;
		}
	if ((nt0 = nt1)) {
		/* Convert time columns in TI */
		*tsqp = tsq = tsx
			? (TIMESTAMP_STRUCT***)tsx
			: (TIMESTAMP_STRUCT***)TM(nc*sizeof(TIMESTAMP_STRUCT**));
		memset(tsq, 0, nc*sizeof(TIMESTAMP_STRUCT**));
		for(i = nt0 = 0; i < nt1; i++) {
			j = zf[i];
			if (j < a)
				nt0++;
			time_convert(h, j);
			}
		if (nt0 && nt1 > 1)
			qsortv(zt, nt1, sizeof(int), icompar, p);
		}

	memset(nk, 0, sizeof(nk));
	u = 0;
	for(i = 0, dbc = dbc0; dbc < dbce; dbc++, i++) {
		dbc->myoffset = j = nk[k = dbc->mytype]++;
		switch(k) {
		 case 0:
			ptr = (PTR)&dd[p[i]];
			break;
		 case 1:
			cc[j] = i;
			ptr = (PTR)(dbc->val0 = cd[j] = Mem(h, dbc->prec));
			if ((i1 = p[i]) < a)
				sa[i1] = cd[j];
			*dbc->val0 = 0;	/* bypass MS Excel bug with empty cells */
			break;
		 case 2:
			ptr = (PTR)&td[j];
		 }
		if (prc(h, "SQLBindCol_4",
				SQLBindCol(hs, ++u, sqlc[k], ptr, dbc->prec, &dbc->len)))
			goto badret;
		}
	if (prc(h, "SQLExecute", SQLExecute(hs)))
		goto badret;
	h->dd = dd;

	/* Allocate new columns even if the existing table is empty. */

	TI->ncols += nn;
	for(i = nt1; i < nnt; i++) {
		db = db0 + (k = p[j = zt[i]]);
		wantsv = dbc0[j].mytype >= 1;
		if (!(wantsv ? db->sval : (char**)db->dval))
			(*TI->ColAlloc)(TI, k, wantsv);
		}

	if (h->oldquotes)
		mix = 3;
	else if ((mix = h->nsmix) == 2) {
		if (ct) {
			for(i = mix = 0; i < nf; ++i) {
				if (ct[i] == 3 && (j = pi[i]) >= 0) {
					dbc = dbc0 + p[i];
					if (dbc->mytype == 1)
						dbc->mytype = 3;
					}
				}
			}
		else
			mix = 1;
		}
	maxrows = TI->maxrows;
	nseen = 0;
	memset(seen = (char*)TM(nrows), 0, nrows);
	while((i = SQLFetch(hs)) == SQL_SUCCESS || i == SQL_SUCCESS_WITH_INFO) {

		/* Trim trailing blanks from character fields, */
		/* (DbLink can supply gratuitous trailing blanks) */
		/* and deal with '...' */

		ra = ra0;
		for(i = nk[1]; i--; ) {
			dbc = dbc0 + (k = cc[i]);
			k = p[k];
			for(s = cs = dbc->val0; *s; s++);
			while(s > cs && *--s == ' ')
				*s = 0;
			if (*cs == 0)
				cs = Missing;
			else if (*cs == '\'' && s > cs && *s == '\'') {
				cs++;
				*s = 0;
				}
			else if (dbc->mytype == 3) {
				if (!(cs = scrunch(h, cs, dd+k, 1)))
					ra = dd;
				}
			else if (mix && mustquote(cs,0) == 2) {
				dd[k] = strtod(cs,0);
				if (mix == 3)
					goto keep_dd;
				else {
					sprintf(dbuf, "%.g", dd[k]);
					if (!strcmp(cs,dbuf)) {
 keep_dd:
						cs = 0;
						ra = dd;
						}
					}
				}
			dbc->val = cs;
			if (k < a)
				sa[k] = cs;
			}

		/* Convert timestamps for Lookup */
		for(i = 0; i < nt0; i++) {
			dbc = dbc0 + (j = zt[i]);
			ts = td + dbc->myoffset;
			t = 100.*ts->year + ts->month;
			t = 100.*t + ts->day;
			t = 100.*t + ts->hour;
			t = 100.*t + ts->minute;
			t = 100.*t + ts->second;
			if (ts->fraction)
				t += ts->fraction / 4294967296.; /* 2^32 */
			dd[p[j]] = t;
			memset(ts, 0, sizeof(TIMESTAMP_STRUCT));
			}
		j = (*TI->Lookup)(ra, sa, TI);
		if (j >= 0) {
			nseen++;
			seen[j] = 1;
			/* Get entries in new columns... */
			for(i = nt1; i < nnt; i++) {
				dbc = dbc0 + (k = zt[i]);
				db = db0 + (k = p[k]);
				switch(dbc->mytype) {
				  case 0:
					db->dval[j] = dd[k];
					break;
				  case 1:
					db->sval[j] = ascrunch(h, dbc->val,
							&db->dval, (long)j, k);
					break;
				  case 2:
					tsp = (TIMESTAMP_STRUCT**)db->sval;
					tsp[j] = ts = (TIMESTAMP_STRUCT*)
							TM(sizeof(TIMESTAMP_STRUCT));
					*ts = td[dbc->myoffset];
					if (ts->year < 0
					 || !memcmp(ts, &Missing_time, sizeof(No_time))
					 || !memcmp(ts, &No_time, sizeof(No_time)))
						*ts = Missing_time; /* kludge*/
					td[dbc->myoffset] = Missing_time;
				  }
				}
			if (a <= 0)
				break;
			}
		else {
			if (nrows >= maxrows)
				maxrows = (*TI->AdjustMaxrows)(TI, 2*maxrows);
			for(i = 0, dbc = dbc0; i < nt; dbc++, i++) {
				db = db0 + (k = p[i]);
				switch(dbc->mytype) {
				  case 0:
					db->dval[nrows] = dd[k];
					break;
				  case 1:
					db->sval[nrows] = ascrunch(h, dbc->val,
							&db->dval, nrows, k);
					break;
				  case 2:
					tsp = (TIMESTAMP_STRUCT**)db->sval;
					tsp[nrows] = ts = (TIMESTAMP_STRUCT*)
							TM(sizeof(TIMESTAMP_STRUCT));
					*ts = td[dbc->myoffset];
					if (!memcmp(ts, &No_time,
							sizeof(TIMESTAMP_STRUCT)))
						*ts = Missing_time; /* kludge*/;
					td[dbc->myoffset] = Missing_time;
				  }
				}
			TI->nrows = ++nrows;
			}
		for(i = nk[1]; i--;)
			*cd[i] = 0;	/* to bypass MS Excel empty-cell bug */
		}
	if (nseen < nrows0 && nn) {
		db = db1 = db0 + nf;
		dbe = db + nn;
		for(i = nf; db < dbe; i++, db++)
			if (!db->sval)
				(*TI->ColAlloc)(TI, i, 1);
		for(m = 0; m < nrows0; m++)
			if (!seen[m])
				for(db = db1; db < dbe; db++)
					db->sval[m] = dbc->mytype == 2
						? (char*)&Missing_time
						: Missing;
		}
	if (tsq) {
		for(i = 0; i < nts; i++) {
			j = p[zt[i]];
			tsq[j] = (TIMESTAMP_STRUCT**)db0[j].sval;
			}
		}
	if (nk[1]) { /* retain text type if specified (added 20030801) */
		sw = 0;
		if (!nrows) {
			sw = *swp = (int*)TM(nc*sizeof(int));
			memset(sw, 0, nc*sizeof(int));
			}
		for(dbc = dbc0; dbc < dbce; dbc++)
			if (dbc->mytype == 1) {
				i = dbc - dbc0;
				db = db0 + (k = p[i]);
				if (sw)
					sw[k] = dbc->prec - 1;
				if (db->sval)
					continue;
				(*TI->ColAlloc)(TI, k, 1);
				}
		}
	rc = 0;
 done:
	prc(h, "SQLFreeStmt", SQLFreeStmt(hs, SQL_CLOSE));
	prc(h, "SQLFreeStmt", SQLFreeStmt(hs, SQL_UNBIND));
 done1:
	nc = TI->arity + TI->ncols;
	if (!p_) {
		p_ = pi = (int*)TM(nc*sizeof(int));
		for(i = 0; i < nc; i++)
			p_[i] = i;
		}
	*pp = p_;
	*pip = pi;
	if (rc)
		goto done2;
	if (ntimes) {
		nrows = TI->nrows;
		if (!tsq) {
			*tsqp = tsq = tsx
				? (TIMESTAMP_STRUCT***)tsx
				: (TIMESTAMP_STRUCT***)
					TM(nc*sizeof(TIMESTAMP_STRUCT**));
			}
		memset(tsq, 0, nc*sizeof(TIMESTAMP_STRUCT**));
		for(j = 0; j < nf; ++j) {
			if (ct[j] == 2)
				tsq[j] = time_convert(h, j);
			}
		}
	if (h->totbadtimes) {
		rc = 1;
		i = h->firstbadcol;
		TI->Errmsg = cs = (char*)TM(strlen(TI->colnames[i]) + 256);
		cs += sprintf(cs, "%d bad time values", h->totbadtimes);
		if (h->firstbadtimes < h->totbadtimes)
			cs += sprintf(cs, ", including %d bad values",
				h->firstbadtimes);
		cs += sprintf(cs, " in column %d (\"%s\").\n", i+1, TI->colnames[i]);
		cs += sprintf(cs,
			"Values should be of the form YYYYMMDDhhmmss; row %ld has\n",
			h->firstbadrow + 1);
		ts = h->ts;
		if (ts->year || ts->month || ts->day)
			cs += sprintf(cs, "\tYYYY = %d\n\tMM = %d\n\tDD = %d\n",
				ts->year, ts->month, ts->day);
		if (ts->hour || ts->minute || ts->second)
			cs += sprintf(cs, "\thh = %d\n\tmm = %d\n\tss = %d\n",
				ts->hour, ts->minute, ts->second);
		}
 done2:
	return rc;
	}

#endif /* NO_Adjust_ampl_odbc */
