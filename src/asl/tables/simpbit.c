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

/* Example of simplified table handlers for .bit and .abt files,   */
/* both of which have the same binary format, which coincides with */
/* the one used by AMPL's builtin "bit" (binary table) handlers.   */

/* The simplified reader does not permit column names to be   */
/* reordered, does not permit "extra" columns, and only deals */
/* with the native binary format.  */

/* The simplified writer makes no attempt to preserve any rows */
/* or columns it is not writing. */

/*
 * The essence of a table handler is to provide one or more table
 * readers and writers.  This example provides one table reader, called
 * Read_ampl_abt, and one table writer, called Write_ampl_ab.  They are
 * static functions (i.e., their names are known only in this file)
 * that are made known to AMPL in the call on add_table_handler in
 * funcadd(), which is defined at the end of this file.
 *
 * Aside from the routines just mentioned, the other routines defined
 * below are auxiliary routines that are helpful in this example.
 *
 * Recall that a table column may contain only numeric values, only
 * symbolic (string) values, or a mixture of numeric and string values.
 * Some of the complexity below is for dealing with string values.
 *
 * All string data are recorded in a single string table, which starts
 * with column names, followed by an array of column types (see comments
 * in the declaration of struct BIntaryTable_header below) and ends with
 * string values for each column.
 *
 * Numeric columns are simply arrays of numeric values.
 * String columns
 */

#include <stdlib.h>
#include <string.h>
#include "funcadd.h"
#include "arith.h"	/* for Arith_Kind_ASL and Long */

#define TM(len) (*ae->Tempmem)(TI->TMI,len)

 static int
cantopen(AmplExports *ae, TableInfo *TI, char *tname)
{
	/* Report failure to open a file. */

	sprintf(TI->Errmsg = (char*)TM(strlen(tname) + 24),
		"Cannot open \"%s\".", tname);
	return DB_Error;
	}

 static char*
tabname(TableInfo *TI)
{
	/* Recognize .abt and .bit table names */
	/* and reject everything else. */

	char *s, *tname;

	switch(TI->nstrings) {
	 case 2:
		if (strcmp(TI->strings[0], "simple-bit"))
			return 0;
		tname = TI->strings[1];
		break;
	 case 1:
		tname = TI->strings[0];
		break;
	 default:
		return 0; /* reject */
	  }
	return !(s = strrchr(tname,'.'))
		|| strcmp(s,".abt") && strcmp(s,".bit") ? 0 : tname;
	}

 typedef struct		/* information for auxiliary routines */
Tinfo {
	AmplExports *ae;
	TableInfo *TI;
	char *tname;
	char **fields;
	char *buf;
	FILE *f;
	unsigned long buflen;
	} Tinfo;

 static int
readerr(Tinfo *ti, char *msg)
{
	/* report an error to AMPL */

	AmplExports *ae = ti->ae;
	TableInfo *TI = ti->TI;
	int len = strlen(ti->tname) + strlen(msg);
	sprintf(TI->Errmsg = (char*)TM(len + 32),
		"Error reading file %s: %s.", ti->tname, msg);
	return DB_Error;
	}

#ifndef Long /* 32-bit int */
#define Long int
#endif

 typedef struct		/* Start of .abt file */
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
	/* write to .abt file, complaining about write failures */

	AmplExports *ae;
	TableInfo *TI;
	int len;

	ae = ti->ae;
	if (fwrite(v, vlen, 1, ti->f) == 1)
		return 0;	/* all went well */

	/* Something went wrong; report failure... */

	TI = ti->TI;
	len = strlen(ti->tname) + strlen(what);
	sprintf(TI->Errmsg = (char*)TM(len + 32),
		"could not write %s to file %s.", what, ti->tname);
	return DB_Error;
	}

static char magic[12] = "\0\nampl.bit\n";

 static int
Write_ampl_abt(AmplExports *ae, TableInfo *TI)
{
	/* .abt table writer */

	BIntaryTable_header bh;
	DbCol *db, *dbe;
	Long *Lt, *Lt0, ns;
	Tinfo ti;
	char *Missing, *s, *s0, *s1, *se, **sp, **sp1;
	int i, nc;
	size_t Ld, Ls;
	unsigned long nr1;

	/* check for name.abt */

	if (!(ti.tname = tabname(TI)))
		return DB_Refuse;

	/* Try to open the file */

	if (!(ti.f = fopen(ti.tname,"wb")))
		return cantopen(ae,TI,ti.tname);

	/* initialize ti and  header (bh) */

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

	/* Compute length of string table */

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

	/* nr1 := smallest even integer >= bh.nrows */

	nr1 = bh.nrows;
	if (nr1 & 1)
		nr1++;

	/* round up strtablen */

	bh.strtablen = (bh.strtablen + 7) & ~7;

	/* allocate string table */

	if (bh.nscols) {
		Lt0 = (Long*)TM(nr1*sizeof(Long));
		Lt0[nr1-1] = 0; /* in case it's not otherwise set */
		}
	s = s0 = (char*)TM(bh.strtablen);
	se = s + bh.strtablen;

	/* record column types -- numeric, string, or both */

	for(db = TI->cols; db < dbe; db++) {
		i = db->dval ? '1' : '0';
		if (db->sval)
			i += 2;
		*s++ = i;
		}
	*s++ = 0;

	/* copy column names to string table */

	sp = TI->colnames;
	for(db = TI->cols; db < dbe; db++)
		for(s1 = *sp++; *s++ = *s1; s1++);

	if (bh.nscols) {

		/* copy string values to string table */

		for(db = TI->cols; db < dbe; db++)
			if (sp = db->sval)
				for(sp1 = sp + bh.nrows; sp < sp1;)
					if ((s1 = *sp++) && s1 != Missing)
						while(*s++ = *s1)
							s1++;
		}

	/* zero fill end of string table */

	while(s < se)
		*s++ = 0;

	/* write the header and string table */

	if (fwr(&bh, sizeof(bh), &ti, "header")
	 || fwr(s0, bh.strtablen, &ti, "strings"))
		return DB_Error;

	/* write the data columns */

	Ld = bh.nrows * sizeof(real);
	Ls = (size_t)(nr1 * sizeof(Long));
	ns = nc;
	for(db = TI->cols; db < dbe; db++) {

		if (db->dval)	/* numeric values for this column */
			if (fwr(db->dval, Ld, &ti, "real column"))
				return DB_Error;

		if (sp = db->sval) {

			/* string values for this column */

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

	/* finish */

	fclose(ti.f);
	return DB_Done;
	}

 static int
frd(void *v, size_t vlen, Tinfo *ti, char *what)
{
	/* Read from .abt file, complaining if the read fails. */

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
	return DB_Error;
	}

 static int
Read_ampl_abt(AmplExports *ae, TableInfo *TI)
{
	/* .abt table reader */

	BIntaryTable_header bh;
	DbCol *db, *dbe;
	Long *si, *sie, *sinfo;
	Tinfo ti;
	char buf[64], *ct, *s, **sp, *stab;
	int a, i, nc, rs;
	real *r;
	size_t Lr, Ls, nr1;

	/* check for name.abt */

	if (!(ti.tname = tabname(TI)))
		return DB_Refuse;

	/* Try to open the file */

	if (!(ti.f = fopen(ti.tname,"rb")))
		return cantopen(ae,TI,ti.tname);

	ti.ae = ae;
	ti.TI = TI;

	/* read header */

	if (frd(&bh, sizeof(bh), &ti, "header"))
		return DB_Error;

	/* check header */

	a = (int)bh.arity;
	if (a != TI->arity) {
		sprintf(buf, "header gives arity %d rather than %d",
			a, TI->arity);
		return readerr(&ti, buf);
		}
	nc = (int)bh.ncols;
	if (nc != TI->ncols) {
		sprintf(buf,
			"header says %d columns; expected %d",
			a + nc, a + TI->ncols);
		return readerr(&ti, buf);
		}

	/* read string table */

	nc += a;
	ct = stab = (char*)TM(bh.strtablen);
	if (frd(stab, bh.strtablen, &ti, "string table"))
		return DB_Error;

	/* allocate memory for numeric and string columns */

	nr1 = (size_t) bh.nrows;
	Lr = nr1 * sizeof(real);
	if (nr1 & 1)
		nr1++;
	Ls = nr1 * sizeof(Long);
	if (bh.nrcols)
		r = (real*)TM(Lr * bh.nrcols);
	if (bh.nscols) {
		sinfo = (Long*)TM(nr1*sizeof(Long));
		sie = sinfo + bh.nrows;
		sp = (char**)TM(bh.nrows*bh.nscols*sizeof(char**));
		}
	db = TI->cols;

	s = stab;

	/* skip past column types */

	while(*s++);

	/* check column names */

	for(i = 0; i < nc; i++) {
		if (strcmp(s,TI->colnames[i]))
			return readerr(&ti, "wrong column names");
		while(*s++);
		}

	/* read the columns */

	for(dbe = db + nc; db < dbe; db++) {
		rs = *ct++ - '0';
		if (rs & 1) {
			/* read numeric values for this column */

			if (frd(db->dval = r, Lr, &ti, "real data"))
				return DB_Error;
			r += bh.nrows;
			}
		else
			db->dval = 0;
		if (rs & 2) {
			/* read string table pointers for this column */

			if (frd(si = sinfo, Ls, &ti, "symbol pointers"))
				return DB_Error;
			db->sval = sp;
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
		else
			db->sval = 0;
		}

	/* send all rows to AMPL at once */

	if ((*TI->AddRows)(TI, TI->cols, bh.nrows))
		return DB_Error;

	/* indicate that we're done... */

	return DB_Done;
	}

 void
funcadd(AmplExports *ae)
{
	/* description of .abt handlers */

	static char info[] = "simple-bit\n"
	"Library file.bit and file.abt (binary table) handler:\n"
	"one or two strings (an optional 'simple-bit' and the file name,\n"
	"ending in \".bit\" or \".abt\") expected before \":[...]\".";

	/* Inform AMPL about the .abt handlers */

	add_table_handler(Read_ampl_abt, Write_ampl_abt, info, 0, 0);
	}
