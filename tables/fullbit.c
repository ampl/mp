/****************************************************************
Copyright (C) 2000 Lucent Technologies
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

/* Example of complete table handlers for .bit and .tab files. */
/* Except for their names and descriptions (specified in the   */
/* funcadd call at the end), they are equivalent to AMPL's     */
/* builtin "tab" and "bit" handlers. */

#include "stdlib.h"
#include "string.h"
#include "funcadd.h"
#include "arith.h"	/* for Arith_Kind_ASL and Long */

#define TM(len) (*ae->Tempmem)(TI->TMI,len)

 static int Adjust_ampl_bit(AmplExports*, TableInfo*, FILE*, char*);
 static int Adjust_ampl_tab(AmplExports*, TableInfo*, FILE*, char*);

 static int
cantopen(AmplExports *ae, TableInfo *TI, char *tname)
{
	sprintf(TI->Errmsg = (char*)TM(strlen(tname) + 24),
		"Cannot open \"%s\".", tname);
	return DB_Error;
	}

 static char*
tabnamebit(TableInfo *TI)
{
	char *s, *tname;

	switch(TI->nstrings) {
	  default: return 0;
	  case 1:
		tname = TI->strings[0];
		break;
	  case 2:
		s = TI->strings[0];
		if (strcmp(s, "lib-bit"))
			return 0;
		tname = TI->strings[1];
	   }
	return !(s = strrchr(tname, '.')) || strcmp(s,".bit") ? 0 : tname;
	}

 static char*
tabnametab(AmplExports *ae, TableInfo *TI)
{
	char *s, *tname;
	size_t L;

	switch(TI->nstrings) {
	  default: return 0;
	  case 1:
		tname = TI->strings[0];
		if (strcmp(tname, "lib-tab"))
			break;
	  case 0:
		L = strlen(tname = TI->tname);
		memcpy(s = (char*)TM(L+5), tname, L);
		strcpy(s+L, ".tab");
		return s;
	  case 2:
		s = TI->strings[0];
		if (strcmp(s, "lib-tab"))
			return 0;
		tname = TI->strings[1];
	   }
	return !(s = strrchr(tname, '.')) || strcmp(s,".tab") ? 0 : tname;
	}

 int
Write_ampl_tab(AmplExports *ae, TableInfo *TI)
{
	DbCol *db, *dbe;
	FILE *f;
	char *Missing, **colnames, *s, **sv, *tname;
	int a, c, i, n, n1, nc;
	long j, nrows;
	real *r;

	if (!(tname = tabnametab(ae,TI)))
		return DB_Refuse;
	if (TI->flags & DBTI_flags_IN
	 && (f = fopen(tname,"r"))
	 && Adjust_ampl_tab(ae, TI, f, tname))
		return DB_Error;
	if (!(f = fopen(tname,"w")))
		return cantopen(ae,TI,tname);
	fprintf(f, "ampl.tab %d %d\n", a = TI->arity, nc = TI->ncols);
	n = a + nc;
	n1 = n - 1;
	colnames = TI->colnames;
	for(i = 0; i < n; i++)
		fprintf(f, "%q%c", colnames[i], i < n1 ? '\t' : '\n');
	dbe = TI->cols + n - 1;
	nrows = TI->nrows;
	Missing = TI->Missing;
	for(j = 0; j < nrows; j++) {
		c = '\t';
		for(db = TI->cols; db <= dbe; db++) {
			if (db == dbe)
				c = '\n';
			if (r = db->dval)
				if ((sv = db->sval) && (s = sv[j]))
 sprint:				if (s == Missing)
 dotprint:					fprintf(f,".%c",c);
					else
						fprintf(f,"%q%c",s,c);
				else
					fprintf(f,"%.g%c",r[j],c);
			else if (s = db->sval[j])
				goto sprint;
			else
				goto dotprint;
			}
		}
	fclose(f);
	return DB_Done;
	}

 typedef struct
Tinfo {
	unsigned long buflen;
	long line;
	AmplExports *ae;
	TableInfo *TI;
	char *tname;
	char **fields;
	char *buf;
	FILE *f;
	int *quoted;
	char *s, *se;	/* scratch */
	int nf;
	} Tinfo;

 static int
readerr(Tinfo *ti, char *msg)
{
	AmplExports *ae = ti->ae;
	TableInfo *TI = ti->TI;
	int len = strlen(ti->tname) + strlen(msg);
	sprintf(TI->Errmsg = (char*)TM(len + 32),
		"Error reading file %s:\n\t%s.", ti->tname, msg);
	fclose(ti->f);
	return DB_Error;
	}

 static int
getnum(char *s, int *np, char **se)
{
	*np = (int)strtol(s,se,10);
	return **se > ' ';
	}

 static int
geteol(AmplExports *ae, FILE *f, char *b, size_t blen)
{
	char *be = b + blen - 1;
	int c;

	for(;;) {
		if ((c = getc(f)) == EOF)
			return 1;
		if (b < be)
			*b++ = c;
		if (c == '\n')
			break;
		}
	*b = 0;
	return 0;
	}

 static void
bgrow(char **fi, Tinfo *ti, char **bp, char **bpe)
{
	AmplExports *ae = ti->ae;
	TableInfo *TI = ti->TI;
	size_t olen = *bpe - *fi;
	unsigned long nlen = 2*ti->buflen;
	char *b = (char*)TM(nlen);
	if (olen)
		memcpy(b, *fi, olen);
	*bp = b + ti->buflen;
	ti->buflen = nlen;
	*fi = ti->buf = b;
	*bpe = b + nlen;
	}

#define bput(x) if (b >= be) bgrow(fi,ti,&b,&be); *b++ = x;

 static int
getfields(Tinfo *ti)
{
	AmplExports *ae = ti->ae;
	FILE *f = ti->f;
	TableInfo *TI = ti->TI;
	char *b, *be, **fi, **fie;
	int c, c1, *q;

	++ti->line;
	b = ti->buf;
	be = b + ti->buflen;
	fi = ti->fields;
	fie = fi + ti->nf;
	q = ti->quoted;
	for(;;q++) {
		while((c = getc(f)) <= ' ') {
			if (c == EOF) {
				if (fi == ti->fields && ti->line > 2) {
					fclose(f);
					return 1;
					}
 eof:
				sprintf(b = ti->buf,
					"early end of file on partial line %ld",
					ti->line);
				goto ret2;
				}
			if (c == '\n') {
 eol:
				sprintf(b = ti->buf,
					"too few fields on line %ld", ti->line);
				goto ret2;
				}
			}
		*fi = b;
		if (c == '\'' || c == '"') for(*q = 1;;) {
			if ((c1 = getc(f)) == EOF)
				goto eof;
			if (c == c1) {
				if ((c1 = getc(f)) == EOF)
					goto eof;
				if (c1 != c) {
					if (c1 > ' ') {
						sprintf(b = ti->buf,
					 "Malformed quoted string, line %ld",
							ti->line);
						goto ret2;
						}
					c = c1;
					break;
					}
				}
			if (c1 == '\n')
				ti->line++;
			bput(c1);
			}
		else {
			*q = 0;
			if (c == '.') {
				if ((c = getc(f)) <= ' ') {
					*fi = 0;
					goto no_bput;
					}
				bput('.');
				}
			bput(c);
			while((c = getc(f)) > ' ') {
				bput(c);
				}
			}
		bput(0);
 no_bput:
		if (++fi == fie)
			break;
		if (c == '\n')
			goto eol;
		}
	while(c != '\n') {
		if ((c = getc(f)) == EOF)
			goto eof;
		if (c > ' ') {
			sprintf(b = ti->buf, "Too many fields on line %ld",
				ti->line);
 ret2:
			readerr(ti,b);
			return 2;
			}
		}
	return 0;
	}
#undef bput

 static char **bletch;

 static int
compar(const void *a, const void *b)
{
	return strcmp(bletch[*(int*)a], bletch[*(int*)b]);
	}

 static void
strsort(int n, char **nam, int *z)
{
	int i;
	for(i = 0; i < n; i++)
		z[i] = i;
	if (n > 1) {
		bletch = nam;
		qsort(z, n, sizeof(int), compar);
		}
	}

 static int
permute(int n1, int n2, int *p, char **nam1, char **nam2, int *z, int *ncp)
{
	/* 0 < n1 <= n2 */
	/* z = scratch */

	int i, j, k, *z2, zi, zj;
	char *s;

	strsort(n1, nam1, z);
	strsort(n2, nam2, z2 = z + n1);
	for(i = j = 0; i < n1; i++) {
		s = nam1[zi = z[i]];
		for(;;j++) {
			if (j >= n2) {
 notfound:
				if (!ncp)
					return zi+1;
				zj = -++*ncp;
				break;
				}
			if (!(k = strcmp(nam2[zj = z2[j]], s)))
				break;
			if (k > 0)
				goto notfound;
			}
		p[zi] = zj;
		}
	return 0;
	}

 int
Read_ampl_tab(AmplExports *ae, TableInfo *TI)
{
	DbCol *db;
	Tinfo ti;
	char buf[64], *s, *se, *tbuf;
	int a, i, j, k, nc, ncf, *p, *z;
	real t;

	if (!(ti.tname = tabnametab(ae,TI)))
		return DB_Refuse;
	if (!(ti.f = fopen(ti.tname,"r")))
		return cantopen(ae,TI,ti.tname);
	ti.ae = ae;
	ti.TI = TI;
	if (geteol(ae, ti.f,buf,sizeof(buf)))
		return readerr(&ti,"missing header line");
	if (strncmp(buf,"ampl.tab ",9)
	 || getnum(buf+9,&a,&s)
	 || getnum(s,&ncf,&s))
		return readerr(&ti,"bad header line");
	if (a != TI->arity) {
		sprintf(buf, "header gives arity %d rather than %d",
			a, TI->arity);
		return readerr(&ti, buf);
		}
	nc = TI->ncols;
	if (ncf < nc) {
		sprintf(buf,
			"header says only %d columns; expected at least %d",
			a + ncf, a + nc);
		return readerr(&ti, buf);
		}
	ti.nf = a + ncf;
	ti.fields = (char**)TM(ti.nf*(sizeof(char*) + 4*sizeof(int)));
	ti.quoted = (int*)(ti.fields + ti.nf);
	p = ti.quoted + ti.nf;
	z = p + ti.nf;
	ti.buf = (char*)TM(ti.buflen = 2000);
	ti.line = 1;
	if (getfields(&ti))
		return DB_Error;

	/* Allow file columns to be permuted */

	if (i = permute(a, a, p, TI->colnames, ti.fields, z, 0)) {
 badfield:
		s = TI->colnames[--i];
		sprintf(tbuf = (char*)TM(strlen(s) + 32),
			"Column %Q not found.", s);
		return readerr(&ti, tbuf);
		}

	if (nc) {
		if (i = permute(nc, ncf, p+a, TI->colnames + a,
				ti.fields + a, z, 0)) {
			i += a;
			goto badfield;
			}
		for(i = 0; i < nc; i++)
			p[i+a] += a;
		}
	TI->nrows = 1;
	nc += a;
	db = TI->cols;
	while(!(k = getfields(&ti))) {
		for(i = 0; i < nc; i++, db++) {
			if (!(s = ti.fields[j = p[i]]))
				db->sval[0] = TI->Missing;
			else if (ti.quoted[j])
				db->sval[0] = s;
			else {
				t = strtod(s, &se);
				if (!*se) {/* valid number */
					db->sval[0] = 0;
					db->dval[0] = t;
					}
				else
					db->sval[0] = s;
				}
			}
		db = TI->cols;
		if ((*TI->AddRows)(TI, db, 1)) {
			fclose(ti.f);
			return DB_Error;
			}
		}
	return k == 1 ? DB_Done : DB_Error;
	}

 static char *
ti_copy(Tinfo *ti, char *s)
{
	AmplExports *ae;
	TableInfo *TI;
	char *rv, *t;
	int n = strlen(s);

	if (ti->se - ti->s <= n) {
		n += 4000;
		ae = ti->ae;
		TI = ti->TI;
		ti->s = (char*)TM(n);
		ti->se = ti->s + n;
		}
	rv = t = ti->s;
	while(*t++ = *s++);
	ti->s = t;
	return rv;
	}

 static int
Adjust_ampl_tab(AmplExports *ae, TableInfo *TI, FILE *f, char *tname)
{
	DbCol *db, *db0, *db1, *dbe;
	Tinfo ti;
	char *Missing, buf[64], *s, *se, **sp, **sv, **sv1, **sv2;
	int a, i, j, k, nc, ncf, nh, nn;
	int *h, *p, *pe, *pe1, *z;
	real *dv, *dv1, *dv2, t;

	ti.ae = ae;
	ti.TI = TI;
	ti.f = f;
	ti.s = ti.se = 0;
	ti.tname = tname;
	if (geteol(ae, f,buf,sizeof(buf)))
		return readerr(&ti,"missing header line");
	if (strncmp(buf,"ampl.tab ",9)
	 || getnum(buf+9,&a,&s)
	 || getnum(s,&ncf,&s))
		return readerr(&ti,"bad header line");
	if (a != TI->arity) {
		sprintf(buf, "header gives arity %d rather than %d",
			a, TI->arity);
		return readerr(&ti, buf);
		}
	nc = TI->ncols;
	ti.nf = a + ncf;
	k = a + nc;
	dv = (real*)TM(ti.nf*(2*sizeof(char*) + 2*sizeof(int) + sizeof(real))
				+ a*(sizeof(real) + sizeof(char*))
				+ (ncf + 2*k)*sizeof(int));
	dv2 = dv + ti.nf;
	sv = (char**)(dv2 + a);
	sv2 = sv + ti.nf;
	ti.fields = sv2 + a;
	ti.quoted = (int*)(ti.fields + ti.nf);
	p = ti.quoted + ti.nf;
	pe = p + k;
	z = pe + ncf;
	ti.buf = (char*)TM(ti.buflen = 2000);
	ti.line = 1;
	if (getfields(&ti))
		return DB_Error;

	/* Allow file columns to be permuted */

	if (i = permute(a, a, p, TI->colnames, ti.fields, z, 0)) {
		s = TI->colnames[--i];
		sprintf(TI->Errmsg = (char*)TM(strlen(s) + 32),
			"Column %Q not found.", s);
		return DB_Error;
		}

	nh = 0;
	h = z + ti.nf;
	pe1 = p + a;
	if (nc) {
		permute(nc, ncf, pe1, TI->colnames + a, ti.fields + a, z, &nh);
		for(i = j = 0; i < nc; i++)
			if (pe1[i] >= 0)
				pe1[i] += a;
			else
				h[j++] = i + a;
		}
	for(i = 0; i < ti.nf; i++)
		z[i] = -1;
	nc += a;
	for(i = 0; i < nc; i++)
		if ((j = p[i]) >= 0)
			z[j] = i;
	nn = 0;
	for(i = a; i < ti.nf; i++)
		if (z[i] < 0)
			nn++;
	sp = TI->colnames;
	s = 0;	/* shut up warning of non-initialization */
	if (nn) {
		k = nc + nn;
		db = (DbCol*)TM(k*(sizeof(DbCol) + sizeof(char*)));
		memcpy(db, TI->cols, nc*sizeof(DbCol));
		memset(db + nc, 0, nn*sizeof(DbCol));
		TI->cols = db;
		sp = (char**)(db + k);
		memcpy(sp, TI->colnames, nc*sizeof(char*));
		TI->colnames = sp;
		TI->ncols += nn;
		}
	j = 0;
	db0 = TI->cols;
	db1 = db0 + nc;
	for(i = a; i < ti.nf; i++)
		if (z[i] < 0) {
			sp[nc] = ti_copy(&ti, ti.fields[i]);
			pe[j++] = i;
			z[i] = nc++;
			}
	dbe = db0 + nc;
	Missing = TI->Missing;
	while(!(k = getfields(&ti))) {
		dv1 = 0;
		sv1 = 0;
		for(i = 0; i < ti.nf; i++) {
			if (!(s = ti.fields[i]))
				s = Missing;
			else if (!ti.quoted[i]) {
				t = strtod(s, &se);
				if (!*se) {/* valid number */
					dv[i] = t;
					s = 0;
					}
				}
			sv[i] = s;
			}
		t = 0.;
		for(i = 0; i < a; i++) {
			j = z[i];
			if (s = sv[i]) {
				if (s == Missing)
					s = ".";
				sv1 = sv2;
				}
			else {
				dv1 = dv2;
				dv2[j] = dv[i];
				}
			sv2[j] = s;
			}
		j = (*TI->Lookup)(dv1, sv1, TI);
		if (j >= 0) {
			pe1 = pe;
			for(db = db1; db < dbe; db++) {
				i = *pe1++;
				if (s = sv[i]) {
					if (!db->sval)
						(*TI->ColAlloc)(TI, db-db0, 1);
					if (s != Missing)
						s = ti_copy(&ti, s);
					db->sval[j] = s;
					}
				else {
					if (!db->dval)
						(*TI->ColAlloc)(TI, db-db0, 0);
					db->dval[j] = dv[i];
					}
				}
			}
		else {
			if (TI->nrows >= TI->maxrows)
				(*TI->AdjustMaxrows)(TI, 2*TI->maxrows);
			j = (int)TI->nrows++;
			for(i = 0; i < ti.nf; i++) {
				db = db0 + z[i];
				if (s = sv[i]) {
					if (!db->sval)
						(*TI->ColAlloc)(TI, z[i], 1);
					if (s != Missing)
						s = ti_copy(&ti, s);
					db->sval[j] = s;
					}
				else {
					if (!db->dval)
						(*TI->ColAlloc)(TI, z[i], 0);
					db->dval[j] = dv[i];
					}
				}
			for(i = 0; i < nh; i++) {
				db = db0 + h[i];
				if (!db->sval)
					(*TI->ColAlloc)(TI, h[i], 1);
				db->sval[j] = Missing;
				}
			}
		if (!a) {
			fclose(f);
			return 0;
			}
		}
	return k != 1;
	}

#ifndef Long /* 32-bit int */
#define Long int
#endif

 typedef struct
BIntaryTable_header {
	char magic[12]; /* \0\nampl.bit\n\0 */
	char arkind[4];	/* kind of arithmetic: "1\0\0\0" or "2\0\0\0" */
			/* for little- or big-endian IEEE arithmetic */
	Long arity;
	Long ncols;	/* not counting first arity columns */
	Long nrows;
	Long nrcols;	/* number of real columns (nonzero DbCol.dval) */
	Long nscols;	/* number of symbolic columns (nonzero DbCol.sval) */
			/* nrcols + nscols >= arity + ncols */
	Long strtablen;	/* length of string table following this header */
			/* String table starts with an array of 1-character */
			/* column descriptors: */
			/* '1' = just real, '2' = just string, '3' = both. */
			/* Then come arity + ncols column names. */
			/* Then come the columns; for type '3', the real */
			/* column comes first.  String columns are arrays */
			/* of Longs, 0 = NULL, 1 = Missing, >= 2 ==> next */
			/* string table entry. */

	} BIntaryTable_header;

 static int
fwr(void *v, size_t vlen, Tinfo *ti, char *what)
{
	AmplExports *ae;
	TableInfo *TI;
	int len;

	ae = ti->ae;
	if (fwrite((const char*)v,vlen, 1, ti->f) == 1)
		return 0;
	fclose(ti->f);
	TI = ti->TI;
	len = strlen(ti->tname) + strlen(what);
	sprintf(TI->Errmsg = (char*)TM(len + 32),
		"could not write %s to file %s.", what, ti->tname);
	return DB_Error;
	}

static char magic[12] = "\0\nampl.bit\n";

 int
Write_ampl_bit(AmplExports *ae, TableInfo *TI)
{
	BIntaryTable_header bh;
	DbCol *db, *dbe;
	FILE *f;
	Long *Lt, *Lt0, ns;
	Tinfo ti;
	char *Missing, *s, *s0, *s1, *se, **sp, **sp1;
	int i, nc;
	size_t Ld, Ls;
	unsigned long nr1;

	if (!(ti.tname = tabnamebit(TI)))
		return DB_Refuse;
	if (TI->flags & DBTI_flags_IN
	 && (f = fopen(ti.tname,"rb"))
	 && Adjust_ampl_bit(ae, TI, f, ti.tname))
		return DB_Error;
	if (!(ti.f = fopen(ti.tname,"wb")))
		return cantopen(ae,TI,ti.tname);
	ti.ae = ae;
	ti.TI = TI;

	memset(&bh, 0, sizeof(bh));
	memcpy(bh.magic, magic, 12);
	i = sprintf(bh.arkind, "%d", Arith_Kind_ASL);
	while(i < 4)
		bh.arkind[i++] = 0;
	bh.arity = TI->arity;
	bh.ncols = TI->ncols;
	bh.nrows = (Long)TI->nrows;
	bh.nrcols = bh.nscols = 0;
	nc = (int)(bh.arity + bh.ncols);
	db = TI->cols;
	for(dbe = db + nc; db < dbe; db++) {
		if (db->dval)
			bh.nrcols++;
		if (db->sval)
			bh.nscols++;
		}
	bh.strtablen = 2*nc + 1;
	sp = TI->colnames;
	sp1 = sp + nc;
	while(sp < sp1)
		bh.strtablen += strlen(*sp++);
	Missing = TI->Missing;
	if (bh.nscols) {
		for(db = TI->cols; db < dbe; db++)
			if (sp = db->sval)
				for(sp1 = sp + bh.nrows; sp < sp1;)
					if ((s = *sp++) && s != Missing)
						bh.strtablen += strlen(s) + 1;
		}
	nr1 = bh.nrows;
	if (nr1 & 1)
		nr1++;
	bh.strtablen = (bh.strtablen + 7) & ~7;
	if (bh.nscols) {
		Lt0 = (Long*)TM(nr1*sizeof(Long));
		Lt0[nr1-1] = 0; /* in case it's not otherwise set */
		}
	s = s0 = (char*)TM(bh.strtablen);
	se = s + bh.strtablen;
	for(db = TI->cols; db < dbe; db++) {
		i = db->dval ? '1' : '0';
		if (db->sval)
			i += 2;
		*s++ = i;
		}
	*s++ = 0;
	sp = TI->colnames;
	for(db = TI->cols; db < dbe; db++)
		for(s1 = *sp++; *s++ = *s1; s1++);
	if (bh.nscols) {
		for(db = TI->cols; db < dbe; db++)
			if (sp = db->sval)
				for(sp1 = sp + bh.nrows; sp < sp1;)
					if ((s1 = *sp++) && s1 != Missing)
						while(*s++ = *s1)
							s1++;
		}
	while(s < se)
		*s++ = 0;
	if (fwr(&bh, sizeof(bh), &ti, "header")
	 || fwr(s0, bh.strtablen, &ti, "strings"))
		return DB_Error;
	Ld = bh.nrows * sizeof(real);
	Ls = (size_t)(nr1 * sizeof(Long));
	ns = nc;
	for(db = TI->cols; db < dbe; db++) {
		if (db->dval)
			if (fwr(db->dval, Ld, &ti, "real column"))
				return DB_Error;
		if (sp = db->sval) {
			Lt = Lt0;
			for(sp1 = sp + bh.nrows; sp < sp1; Lt++) {
				*Lt = 0;
				if (s = *sp++)
					*Lt = s == Missing ? 1 : ++ns;
				}
			if (fwr(Lt0, Ls, &ti, "string column"))
				return DB_Error;
			}
		}
	fclose(ti.f);
	return DB_Done;
	}

 static int
frd(void *v, size_t vlen, Tinfo *ti, char *what)
{
	AmplExports *ae;
	TableInfo *TI;
	int len;

	ae = ti->ae;
	if (fread(v, vlen, 1, ti->f) == 1)
		return 0;
	TI = ti->TI;
	len = strlen(ti->tname) + strlen(what);
	sprintf(TI->Errmsg = (char*)TM(len + 32),
		"could not read %s from file %s.", what, ti->tname);
	fclose(ti->f);
	return DB_Error;
	}

 static void
Lswap(Long *x, size_t L)
{
	char *s = (char*)x;
	char *se = s + L;
	char *s1, *s2;
	int i;

	while(s < se) {
		s1 = s;
		s2 = s + 3;
		s = s2 + 1;
		do {
			i = *s1;
			*s1++ = *s2;
			*s2-- = i;
			}
			while(s1 < s2);
		}
	}

 static void
rswap(real *x, size_t L)
{
	char *s = (char*)x;
	char *se = s + L;
	char *s1, *s2;
	int i;

	while(s < se) {
		s1 = s;
		s2 = s + 7;
		s = s2 + 1;
		do {
			i = *s1;
			*s1++ = *s2;
			*s2-- = i;
			}
			while(s1 < s2);
		}
	}

 static void
headswap(BIntaryTable_header *bh)
{
	Lswap(&bh->arity, 6*sizeof(Long));
	}

 int
Read_ampl_bit(AmplExports *ae, TableInfo *TI)
{
	BIntaryTable_header bh;
	DbCol *db;
	Long *si, *sie, *sinfo;
	Tinfo ti;
	char buf[64], *ct, *s, **sp, **sp1, *stab, *tbuf;
	int a, bswap, i, j, je, k, nc, nr, ns, *p, rs, *z, *zs;
	real *r, *r1, x;
	size_t Lr, Ls, nr1;

	if (!(ti.tname = tabnamebit(TI)))
		return DB_Refuse;
	if (!(ti.f = fopen(ti.tname,"rb")))
		return cantopen(ae,TI,ti.tname);
	ti.ae = ae;
	ti.TI = TI;

	if (frd(&bh, sizeof(bh), &ti, "header"))
		return DB_Error;
	bswap = 0;
	x = strtod(bh.arkind,&s);
	if (x != Arith_Kind_ASL) {
		if (x == 3. - Arith_Kind_ASL) {
			bswap = 1;
			headswap(&bh);
			}
		else {
			sprintf(buf, "header gives unexpected arkind = %g", x);
			return readerr(&ti,buf);
			}
		}
	ti.nf = (int)(bh.arity + bh.ncols);
	a = (int)bh.arity;
	if (a != TI->arity) {
		sprintf(buf, "header gives arity %d rather than %d",
			a, TI->arity);
		return readerr(&ti, buf);
		}
	nc = TI->ncols;
	if (ti.nf < nc) {
		sprintf(buf,
			"header says only %d columns; expected at least %d",
			a + ti.nf, a + nc);
		return readerr(&ti, buf);
		}
	ti.fields = (char**)TM(bh.strtablen + ti.nf*sizeof(char*));
	s = ct = stab = (char*)(ti.fields + ti.nf);
	if (frd(stab, bh.strtablen, &ti, "string table"))
		return DB_Error;
	while(*s++);
	for(i = 0; i < ti.nf; i++) {
		ti.fields[i] = s;
		while(*s++);
		}

	/* Allow file columns to be permuted */

	p = (int*)TM((3*ti.nf + 1)*sizeof(int));
	z = p + ti.nf;
	if (i = permute(a, a, p, TI->colnames, ti.fields, z, 0)) {
 badfield:
		s = TI->colnames[--i];
		sprintf(tbuf = (char*)TM(strlen(s) + 32),
			"Column %Q not found.", s);
		return readerr(&ti, tbuf);
		}

	if (nc) {
		if (i = permute(nc, ti.nf - a, p+a, TI->colnames + a,
				ti.fields + a, z, 0)) {
			i += a;
			goto badfield;
			}
		for(i = 0; i < nc; i++)
			p[i+a] += a;
		}
	nc += a;
	zs = z + ti.nf;
	for(i = 0; i < ti.nf; i++)
		zs[i] = 0;
	/* This would allow duplicate column names. */
	for(i = 0; i < nc; i++)
		zs[p[i]]++;
	for(i = j = 0; i < ti.nf; i++)
		j = zs[i] += j;
	zs[ti.nf] = nc;
	for(i = 0; i < nc; i++)
		z[--zs[p[i]]] = i;
	nr = ns = 0;
	for(i = j = 0; i < ti.nf;) {
		je = zs[++i];
		while(j < je) {
			k = z[j++];
			k = ct[k] - '0';
			if (k & 1)
				nr++;
			if (k & 2)
				ns++;
			}
		}
	nr1 = (size_t) bh.nrows;
	if (bh.nrcols) {
		Lr = nr1 * sizeof(real);
		if (nr < bh.nrcols)
			nr++;	/* elbow room for omitting columns */
		r = (real*)TM(Lr * nr);
		}
	if (bh.nscols) {
		if (nr1 & 1)
			nr1++;
		Ls = nr1 * sizeof(Long);
		sinfo = (Long*)TM(nr1*sizeof(Long));
		sie = sinfo + bh.nrows;
		if (ns)
			sp = (char**)TM(bh.nrows*ns*sizeof(char**));
		}
	for(i = j = k = 0;;) {
		je = zs[++i];
		rs = *ct++ - '0';
		r1 = 0;
		sp1 = 0;
		if (rs & 1) {
			if (frd(r1 = r, Lr, &ti, "real data"))
				return DB_Error;
			if (bswap)
				rswap(r1, Lr);
			}
		if (rs & 2) {
			if (frd(si = sinfo, Ls, &ti, "symbol pointers"))
				return DB_Error;
			if (j >= je)
				continue;
			if (bswap)
				Lswap(si, Ls);
			sp1 = sp;
			while(si < sie)
				switch(*si++) {
				 case 0:
					*sp++ = 0;
					break;
				 case 1:
					*sp++ = TI->Missing;
					break;
				 default:
					*sp++ = s;
					while(*s++);
				 }
			}
		if (r1 && j < je)
			r += bh.nrows;
		while(j < je) {
			db = TI->cols + z[j++];
			db->dval = r1;
			db->sval = sp1;
			if (++k >= nc)
				goto done;
			}
		}
 done:
	fclose(ti.f);
	if ((*TI->AddRows)(TI, TI->cols, bh.nrows))
		return DB_Error;
	return DB_Done;
	}

 static int
Adjust_ampl_bit(AmplExports*ae, TableInfo *TI, FILE *f, char *tname)
{
	BIntaryTable_header bh;
	DbCol *db, *db1, *dbh;
	Long *si, *sie, *sinfo;
	Tinfo ti;
	char *Missing, buf[64];
	char *ct, *s, *s1, **sp, **sp0, **sp00, **sp1, **spi, *stab;
	int a, bswap, i, j, j1, k, nc, nc0, nh, nn, nnr;
	int nor, nr, nr0, nrf, ns, rs;
	int *h, *p, *pe1, *q, *qh, *qn, *z, *zs;
	real *r, *r0, *r00, *r1, x;
	size_t Lr, Ls, nr1;

	ti.ae = ae;
	ti.TI = TI;
	ti.tname = tname;
	ti.f = f;

	if (frd(&bh, sizeof(bh), &ti, "header"))
		return DB_Error;
	bswap = 0;
	x = strtod(bh.arkind,&s);
	if (x != Arith_Kind_ASL) {
		if (x == 3. - Arith_Kind_ASL) {
			bswap = 1;
			headswap(&bh);
			}
		else {
			sprintf(buf, "header gives unexpected arkind = %g", x);
			return readerr(&ti,buf);
			}
		}
	ti.nf = (int)(bh.arity + bh.ncols);
	a = (int)bh.arity;
	if (a != TI->arity) {
		sprintf(buf, "header gives arity %d rather than %d",
			a, TI->arity);
		return readerr(&ti, buf);
		}
	nc = TI->ncols;
	if (ti.nf < nc) {
		sprintf(buf,
			"header says only %d columns; expected at least %d",
			a + ti.nf, a + nc);
		return readerr(&ti, buf);
		}
	ti.fields = (char**)TM(bh.strtablen + ti.nf*sizeof(char*));
	s = ct = stab = (char*)(ti.fields + ti.nf);
	if (frd(stab, bh.strtablen, &ti, "string table"))
		return DB_Error;
	while(*s++);
	for(i = 0; i < ti.nf; i++) {
		ti.fields[i] = s;
		while(*s++);
		}

	/* Allow file columns to be permuted */

	k = ti.nf + nc;
	dbh = (DbCol*)TM(a*sizeof(DbCol) + (3*k + 1 + nc)*sizeof(int));
	p = (int*)(dbh + a);
	h = p + k;
	z = h + nc;
	if (i = permute(a, a, p, TI->colnames, ti.fields, z, 0)) {
		s = TI->colnames[--i];
		sprintf(TI->Errmsg = (char*)TM(strlen(s) + 32),
			"Column %Q not found.", s);
		return DB_Error;
		}

	nh = 0;
	pe1 = p + a;
	nn = 0;
	if (nc) {
		permute(nc, ti.nf - a, pe1, TI->colnames + a,
			ti.fields + a, z, &nh);
		for(i = 0; i < nc; i++)
			if (pe1[i] >= 0)
				pe1[i] += a;
			else
				h[nn++] = i + a;
		}
	nc0 = nc += a;
	for(i = 0; i < ti.nf; i++)
		z[i] = -1;
	for(i = 0; i < nc; i++)
		if ((j = p[i]) >= 0)
			z[j] = i;
	zs = z + ti.nf;
	for(i = a; i < nc; i++)
		zs[i] = -1;
	for(i = a; i < ti.nf; i++)
		if ((j = z[i]) < 0)
			z[i] = nc++;	/* columns in bh but not in TI */
	nr = ns = 0;
	for(i = 0; i < a;) {
		k = z[i++];
		k = ct[k] - '0';
		if (k & 1)
			nr++;
		if (k & 2)
			ns++;
		}
	nr1 = (size_t) bh.nrows;
	if (bh.nrcols) {
		Lr = nr1 * sizeof(real);
		if (!nr)
			nr = 1;
		r = r00 = (real*)TM(Lr * nr);
		}
	if (bh.nscols) {
		if (nr1 & 1)
			nr1++;
		Ls = nr1 * sizeof(Long);
		sinfo = (Long*)TM(nr1*sizeof(Long));
		sie = sinfo + bh.nrows;
		if (!ns)
			ns = 1;
		sp = sp00 = (char**)TM(bh.nrows*ns*sizeof(char**));
		}
	Missing = TI->Missing;
	for(i = 0; i < a; i++) {
		rs = *ct++ - '0';
		r1 = 0;
		sp1 = 0;
		if (rs & 1) {
			if (frd(r1 = r, Lr, &ti, "real data"))
				return DB_Error;
			if (bswap)
				rswap(r1, Lr);
			}
		if (rs & 2) {
			if (frd(si = sinfo, Ls, &ti, "symbol pointers"))
				return DB_Error;
			if (bswap)
				Lswap(si, Ls);
			sp1 = sp;
			while(si < sie)
				switch(*si++) {
				 case 0:
					*sp++ = 0;
					break;
				 case 1:
					*sp++ = Missing;
					break;
				 default:
					*sp++ = s;
					while(*s++);
				 }
			}
		db = dbh + p[i];
		if (db->dval = r1)
			r += bh.nrows;
		db->sval = sp1;
		}
	nrf = (int)bh.nrows;
	nr = nr0 = (int)TI->nrows;
	r0 = (real*)TM(a*(sizeof(real) + sizeof(char*))
			+ (nrf + nr)*sizeof(int));
	memset(r0, 0, a*sizeof(real));
	sp0 = (char**)(r0 + a);
	q = (int*)(sp0 + a);
	qh = q + nrf;
	for(i = 0; i < nr; i++)
		qh[i] = -1;
	for(i = 0; i < nrf; i++) {
		r1 = 0;
		sp1 = 0;
		db = dbh;
		for(j = 0; j < a; j++, db++) {
			if ((spi = db->sval)
			 && (s1 = spi[i])
			 && s1 != Missing) {
				if (!sp1)
					memset(sp1 = sp0, 0, a*sizeof(char*));
				sp1[j] = s1;
				}
			if (db->dval) {
				r1 = r0;
				r1[j] = db->dval[i];
				}
			}
		if ((j = (*TI->Lookup)(r1, sp1, TI)) >= 0)
			qh[q[i] = j] = i;
		else
			q[i] = nr++;
		}
	if (nr > TI->maxrows)
		(*TI->AdjustMaxrows)(TI, nr);
	if (nnr = nr - nr0) {
		qn = (int*)TM(nnr*sizeof(int));
		for(i = j = 0; i < nrf; i++)
			if (q[i] >= nr0) {
				qn[j++] = i;
				db = dbh;
				db1 = TI->cols;
				for(k = 0; k < a; k++, db++, db1++) {
					if ((spi = db->sval)
					 && (s1 = spi[i])
					 && s1 != Missing) {
						if (!(sp1 = db1->sval))
						 sp1 = (char**)
						 (*TI->ColAlloc)(TI, k, 1);
						sp1[q[i]] = s1;
						continue;
						}
					if (!(r1 = db1->dval))
						r1 = (real*)
						 (*TI->ColAlloc)(TI, k, 0);
					r1[q[i]] = db->dval[i];
					}
				}
		TI->nrows = nr;
		dbh = TI->cols;
		for(i = 0; i < nn; i++) {
			db = dbh + (j = h[i]);
			if (!(sp1 = db->sval))
				sp1 = (char**)(*TI->ColAlloc)(TI, j, 1);
			for(j = nr0; j < nr; j++)
				sp1[j] = Missing;
			}
		}
	if (nc > nc0) {
		db1 = (DbCol*)TM(nc*(sizeof(DbCol)+sizeof(char*)));
		sp1 = (char**)(db1 + nc);
		memcpy(db1, TI->cols, nc0*sizeof(DbCol));
		memset(db1+nc0, 0, (nc-nc0)*sizeof(DbCol));
		memcpy(sp1, TI->colnames, nc0*sizeof(char*));
		TI->cols = db1;
		TI->colnames = sp1;
		j = 0;
		for(i = nc0; i < nc; i++) {
			while(z[j] < nc0)
				j++;
			sp1[i] = ti.fields[j++];
			}
		TI->ncols = nc - a;
		}
	else if (!nnr)
		return 0;
	nor = 0;
	for(i = 0; i < nr0; i++)
		if (qh[i] < 0)
			qh[nor++] = i;
	dbh = TI->cols;
	if (nor)
		for(i = nc0; i < nc; i++) {
			db = dbh + i;
			if (!(sp1 = db->sval))
				sp1 = (char**)(*TI->ColAlloc)(TI, i, 1);
			for(j = 0; j < nor; j++)
				sp1[qh[j]] = Missing;
			}
	r = r00;
	sp = sp00;
	for(i = a; i < ti.nf; i++) {
		rs = *ct++ - '0';
		j = z[i];
		db = dbh + j;
		if (rs & 1) {
			if (frd(r, Lr, &ti, "real data"))
				return DB_Error;
			if (bswap)
				rswap(r, Lr);
			if (j >= nc0 || nnr) {
				if (!(r1 = db->dval))
					r1 = (real*)(*TI->ColAlloc)(TI, j, 0);
				k = 0;
				if (j >= nc0)
					for(; k < nrf; k++)
						r1[q[k]] = r[k];
				else
					for(; k < nnr; k++) {
						j1 = qn[k];
						r1[q[j1]] = r[j1];
						}
				}
			}
		if (rs & 2) {
			if (frd(si = sinfo, Ls, &ti, "symbol pointers"))
				return DB_Error;
			if (j < nc0 && !nnr)
				continue;
			if (bswap)
				Lswap(si, Ls);
			if (!(sp1 = db->sval))
				sp1 = (char**)(*TI->ColAlloc)(TI, j, 1);
			sp = sp00;
			while(si < sie)
				switch(*si++) {
				 case 0:
					*sp++ = 0;
					break;
				 case 1:
					*sp++ = TI->Missing;
					break;
				 default:
					*sp++ = s;
					while(*s++);
				 }
			k = 0;
			if (j >= nc0)
				for(; k < nrf; k++)
					sp1[q[k]] = sp00[k];
			else
				for(; k < nnr; k++) {
					j1 = qn[k];
					sp1[q[j1]] = sp00[j1];
					}
			}
		}
	return 0;
	}

 void
funcadd(AmplExports *ae)
{
	static char bitdesc[] = "lib-bit\n"
	"Library file.bit (binary table) handler: one or two strings\n"
	"(an option 'lib-bit' and the file name, ending in \".bit\")\n"
	"expected before \":[...]\".";

	static char tabdesc[] = "lib-tab\n"
	"Library file.tab (ASCII table) handler: at most one or two strings\n"
	"(an optional 'lib-tab' and the file name, ending in \".tab\")\n"
	"expected before \":[...]\";\n"
	"table_name.tab is assumed if there are no strings or only 'lib-tab'.";

	add_table_handler(Read_ampl_tab, Write_ampl_tab, tabdesc, 0, 0);
	add_table_handler(Read_ampl_bit, Write_ampl_bit, bitdesc, 0, 0);
	}
