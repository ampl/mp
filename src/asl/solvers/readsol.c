/****************************************************************
Copyright (C) 1997, 1999,2000 Lucent Technologies
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

#include "nlp.h"

#define MSGGULP 1024

 typedef struct
msginfo {
	char *msg;
	char *msg0;
	char *msgend;
	ftnlen msglen;
	} msginfo;

 static void
msgput(msginfo *m, const char *b, int n)
{
	char *msg, *msg0, *msgend;
	const char *be;
	ftnlen msglen0;

	msg = m->msg;
	msgend = m->msgend;
	be = b + n;
	while(b < be) {
		if (msg >= msgend) {
			msglen0 = m->msglen;
			msg0 = m->msg0 = (char*)
				Realloc(m->msg0, m->msglen += MSGGULP);
			msg = msg0 + msglen0;
			m->msgend = msg0 + m->msglen;
			}
		*msg++ = *b++;
		}
	m->msg = msg;
	}

 static void
badnumber(ASL *asl, fint a, fint b, const char *what)
{
	fprintf(Stderr, "%s indicates %ld rather than %ld %s\n",
		filename, (long)a, (long)b, what);
	fflush(Stderr);
	}

 static int
decstring(char *buf, real *val)
{
	char *be;
	int c;

	*val = strtod(buf, &be);
	return be <= buf || (((c = be[-1]) < '0' || c > '9') && c != '.');
	}

 char *
read_sol_ASL(ASL *asl, real **xp, real **yp)
{
	int binary, flag1, i, j, je, n, need_vbtol;
	FILE *f;
	char buf[512], *s, *se;
	real vbtol, *y;
	ftnlen L, L1, L2;
	fint Objno[2], Options[14], nOpts, *z;
	msginfo mi;

	if (!asl || asl->i.ASLtype < 1 || asl->i.ASLtype > 5)
		badasl_ASL(asl,0,"read_soln");
	strcpy(stub_end, ".sol");
	f = fopen(filename, "rb");
	if (!f) {
		fprintf(Stderr, "Can't open %s\n", filename);
		fflush(Stderr);
		return 0;
		}
	if (fread(&L, sizeof(ftnlen), 1, f) && L == 6) {
		/* binary files may be written by Fortran unformatted writes */
		binary = 1;
		if (!fread(buf, 6, 1, f)
		 || strncmp(buf,"binary",6)
		 || !fread(&L, sizeof(ftnlen), 1, f)
		 || L != 6) {
 badbinary:
			fprintf(Stderr, "bad binary file %s\n", filename);
			fflush(Stderr);
			goto done;
			}
		}
	else {
		binary = 0;
		rewind(f);
		}

	L1 = 0;
	z = 0;

	/* Read termination msg */
	nOpts = i = need_vbtol = 0;
	mi.msg = mi.msg0 = (char *)Malloc(mi.msglen = MSGGULP);
	mi.msgend = mi.msg0 + MSGGULP;
	if (binary) {
		for(;; i++) {
			if (!fread(&L,sizeof(ftnlen),1,f))
				goto early_eof;
			if ((L1 = L)) {
				do {
					n = L < sizeof(buf) ? (int)L : (int)sizeof(buf);
					L -= n;
					if (!fread(buf, n, 1, f))
						goto early_eof;
					if (!L) {
						while(--n >= 0 && buf[n] == ' ');
						n++;
						}
					msgput(&mi, buf, n);
					}
					while(L);
				msgput(&mi, "\n", 1);
				}
			if (!fread(&L, sizeof(ftnlen), 1, f))
				goto early_eof;
			if (L != L1)
				goto badbinary;
			if (!L)
				break;
			}
		L1 = n_con * sizeof(real);
		if (!fread(&L, sizeof(ftnlen), 1, f))
			goto badbinary;
		if (L >= 3*sizeof(fint) + 7) {
			/* check for Options */
			if (!fread(buf, 7, 1, f))
				goto badbinary;
			if (strncmp(buf, "Options", 7))
				goto badbinary;
			if (!fread(&Options, sizeof(fint), 3, f))
				goto badbinary;
			nOpts = Options[0];
			if (nOpts < 3 || nOpts > 9) {
 bad_nOpts:
				fprintf(Stderr,
				"expected nOpts between 3 and 9; got %ld: ",
					(long)nOpts);
				goto badbinary;
				}
			if (Options[2] == 3) {
				nOpts -= 2;
				need_vbtol = 1;
				}
			L2 = (nOpts + 5)*sizeof(fint) + 7;
			if (need_vbtol)
				L2 += sizeof(real);
			if (L != L2)
				goto badbinary;
			if (!fread(Options+3, sizeof(fint), nOpts + 2, f))
				goto badbinary;
			if (need_vbtol
			 && !fread(&vbtol, sizeof(real), 1, f))
				goto badbinary;
			if (!fread(&L2, sizeof(ftnlen), 1, f)
				|| L != L2)
				goto badbinary;
			}
		else if (L != L1)
			goto badbinary;
		}
	else {
		for(;; i++) {
			if (!fgets(buf, sizeof(buf), f)) {
 early_eof:
				fprintf(Stderr,
					"early end of file reading %s\n",
					filename);
				fflush(Stderr);
 done:
				fclose(f);
				return 0;
				}
			if (*buf == '\n' || (*buf == '\r' && buf[1] == '\n'))
				break;
			msgput(&mi, buf, strlen(buf));
			}
		while((j = getc(f)) == '\n');
		if (j != 'O')
			ungetc(j,f);
		else {
			if (!fgets(buf, sizeof(buf), f))
				goto early_eof;
			if (!strncmp(buf, "ptions", 6)) {
				for(j = 0; j <3; j++) {
					if (!fgets(buf, sizeof(buf), f))
						goto early_eof;
					Options[j] = strtol(buf,&se,10);
					if (se == buf)
						goto badline;
					}
				nOpts = Options[0];
				if (nOpts < 3 || nOpts > 9)
					goto bad_nOpts;
				if (Options[2] == 3) {
					nOpts -= 2;
					need_vbtol = 1;
					}
				je = (int)nOpts + 4;
				for(j = 3; j <= je; j++) {
					if (!fgets(buf, sizeof(buf), f))
						goto early_eof;
					Options[j] = strtol(buf,&se,10);
					if (se == buf)
						goto badline;
					}
				if (need_vbtol && !fgets(buf, sizeof(buf), f))
					goto early_eof;
				/* We don't do anything here with vbtol, */
				/* so we don't bother converting it. */
				}
			}
		}
	memcpy(ampl_options, Options, (nOpts+1)*sizeof(fint));
	msgput(&mi, "", 1);	/* add null to end */

	if (i)
		fflush(stdout);

	if (nOpts) {
		z = Options + nOpts + 1;
		j = (int)z[3];
		if (j > n_var || j < 0) {
			badnumber(asl, j, n_var, "variables");
			goto done;
			}
		j = (int)z[1];
		if (j > n_con || j < 0) {
			badnumber(asl, j, n_con, "constraints");
			goto done;
			}
		if (binary) {
			L1 = j * sizeof(real);
			if (!fread(&L, sizeof(ftnlen), 1, f))
				goto badbinary;
			if (L != L1)
				goto badbinary;
			}
		}
	else
		j = n_con;
	if (!j) {
		*yp = 0;
		goto get_x;
		}
	y = *yp = (real *)Malloc(n_con * sizeof(real));
	if (binary) {
		if (fread(y, sizeof(real), j, f) != j)
			goto early_eof;
		if (!fread(&L, sizeof(ftnlen), 1, f) || L != L1)
			goto badbinary;
		y += j;
		}
	else for(i = 0; i < j; i++) {
		if (!fgets(buf, sizeof(buf), f))
			goto early_eof;
		if (!decstring(buf, y++))
			continue;
 badline:
		fprintf(Stderr, "bad line in %s: %s", filename, buf);
		fflush(Stderr);
		goto done;
		}
	y = *yp;
	while(j < n_con)
		y[j++] = 0;
 get_x:
	Objno[0] = 0;
	Objno[1] = -1;
	flag1 = asl->i.flags & ~1;;
	if (!(j = nOpts ? (int)z[3] : n_var)) {
		*xp = 0;
		goto ret;
		}
	y = *xp = (real *)Malloc(n_var*sizeof(real));

	if (binary) {
		L1 = j * sizeof(real);
		if (!fread(&L, sizeof(ftnlen), 1, f) || L != L1)
			goto badbinary;
		if (fread(y, sizeof(real), j, f) != j)
			goto early_eof;
		y += j;
		/* do we have an obj_no ? */
		if (!fread(&L, sizeof(ftnlen), 1, f) || L != L1)
			goto badbinary;
		if (fread(&L, sizeof(ftnlen), 1, f)) {
			i = 1;
			if (L == 2*sizeof(fint)) {
				i = 2;
				flag1 |= 1;
				}
			else if (L != sizeof(fint))
				goto badbinary;
			if (!fread(Objno, i*sizeof(fint), 1, f))
				goto badbinary;
			}
		}
	else {
		for(i = j; i > 0; i--) {
			if (!fgets(buf, sizeof(buf), f))
				goto early_eof;
			if (decstring(buf, y++))
				goto badline;
			}
		if (fgets(buf,sizeof(buf), f)) {
			if (strncmp(buf,"objno ",6)) {
 extra_line:
				fprintf(Stderr, "Bug: extra line in %s:\n%s",
					filename, buf);
				fflush(Stderr);
				}
			else {
				Objno[0] = strtol(buf+6, &se, 10);
				if (se == buf+6 || *se > ' ')
					goto extra_line;
				if (*se == ' ') {
					Objno[1] = strtol(se,&s,10);
					if (s == se || *s > ' ')
						goto extra_line;
					flag1 |= 1;
					}
				}
			}
		}
	asl->i.flags = flag1;
	obj_no = (int)Objno[0];
	solve_result_num = (int)Objno[1];
	y = *xp;
	while(j < n_var)
		y[j++] = 0;
 ret:
	fclose(f);
	return (char*)Realloc(mi.msg0, mi.msglen);
	}
