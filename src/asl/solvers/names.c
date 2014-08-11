/****************************************************************
Copyright (C) 1997, 1999-2001 Lucent Technologies
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

 static char **
get_names(ASL *asl, const char *suf, int no, int n0, int n1, int n)
{
	FILE *f;
	char buf[512], **nam, **ne, **rv, *s;
	int nt;

	nt = n1 + no;
	nam = rv = (char**)mem(nt*sizeof(char*));
	ne = nam + nt;
	strcpy(stub_end, suf);
	if ((f = fopen(filename, "r"))) {
		while(nam < ne && fgets(buf,sizeof(buf),f)) {
			for(s = buf; *s && *s != '\n'; s++);
			*s = 0;
			strcpy(*nam++ = (char*)mem(s-buf+1), buf);
			}
		fclose(f);
		}
	while(nam < ne)
		*nam++ = 0;
	return rv;
	}

 static void
get_row_names(ASL *asl)
{
	asl->i.connames = get_names(asl, ".row", n_obj + n_lcon,
				asl->i.n_con0, asl->i.n_con1, n_con);
	asl->i.lconnames = asl->i.connames + n_con;
	asl->i.objnames = asl->i.lconnames + n_lcon;
	}

static char badconname[] = "**con_name(bad n)**";
static char badvarname[] = "**var_name(bad n)**";

 char *
con_name_nomap_ASL(ASL *asl, int n, int *p)
{
	char buf[32], **np, *rv;
	const char *fmt;

	if (n < 0 || n >= asl->i.n_con1)
		return badconname;
	if (!asl->i.connames)
		get_row_names(asl);
	np = asl->i.connames + n;
	if (!(rv = *np)) {
		fmt = "_scon[%d]";
		if (p && p[n] < 0)
			fmt = "_scon_aux[%d]";
		*np = rv = (char*)mem(Sprintf(buf,fmt,n+1)+1);
		strcpy(rv, buf);
		}
	return rv;
	}

 char *
con_name_ASL(ASL *asl, int n)
{
	int k, *p;

	if (n < 0 || n >= n_con)
		return badconname;
	if ((p = asl->i.cmap)) {
		if ((k = p[n]) >= 0 && k < asl->i.n_con1)
			n = k;
		}
	return con_name_nomap_ASL(asl, n, p);
	}

 char *
lcon_name_ASL(ASL *asl, int n)
{
	char buf[32], **np, *rv;
	static char badlconname[] = "**lcon_name(bad n)**";

	if (n < 0 || n >= n_lcon)
		return badlconname;
	if (!asl->i.lconnames)
		get_row_names(asl);
	np = asl->i.lconnames + n;
	if (!(rv = *np)) {
		*np = rv = (char*)mem(Sprintf(buf,"_slogcon[%d]",n+1)+1);
		strcpy(rv, buf);
		}
	return rv;
	}

 char *
obj_name_ASL(ASL *asl, int n)
{
	char buf[32], **np, *rv;
	static char badobjname[] = "**obj_name(bad n)**";

	if (n < 0 || n >= n_obj)
		return badobjname;
	if (!asl->i.objnames)
		get_row_names(asl);
	np = asl->i.objnames + n;
	if (!(rv = *np)) {
		*np = rv = (char*)mem(Sprintf(buf,"_sobj[%d]",n+1)+1);
		strcpy(rv, buf);
		}
	return rv;
	}

 char *
var_name_nomap_ASL(ASL *asl, int n, int *p /* not used */)
{
	char buf[32], **np, *rv;

	if (n < 0 || n >= asl->i.n_var1)
		return badvarname;
	if (!asl->i.varnames)
		asl->i.varnames = get_names(asl, ".col", 0,
				asl->i.n_var0, asl->i.n_var1, n_var);
	np = asl->i.varnames + n;
	if (!(rv = *np)) {
		*np = rv = (char*)mem(Sprintf(buf,"_svar[%d]",n+1)+1);
		strcpy(rv, buf);
		}
	return rv;
	}

 char *
var_name_ASL(ASL *asl, int n)
{
	int k, *p;

	if (n < 0 || n >= n_var)
		return badvarname;
	if ((p = asl->i.vmap)) {
		if ((k = p[n]) >= 0 && k < asl->i.n_var1)
			n = k;
		}
	return var_name_nomap_ASL(asl, n, 0);
	}
