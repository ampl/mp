/****************************************************************
Copyright (C) 1997, 1999, 2000 Lucent Technologies
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

#include "getstub.h"

 typedef struct
SufHead {
	char sufid[8];
	fint kind;
	fint n;
	fint namelen;
	fint tablen;
	} SufHead;

typedef char *(*Name)(ASL*,int,int*);

 static void
getsufhead(ASL *asl, SufDesc *d, SufHead *sh, int *np, int **zp)
{
	int i, *ip, *ipe, *map, n, nz;
	real *rp, *rpe;

	memcpy(sh->sufid, "\nSuffix\n", 8);
	sh->kind = d->kind &
		(ASL_Sufkind_mask | ASL_Sufkind_real | ASL_Sufkind_iodcl);
	*np = n = (&asl->i.n_var_)[i = d->kind & ASL_Sufkind_mask];
	*zp = map = i < 2 ? (&asl->i.vmap)[i] : 0;
	nz = 0;
	if (d->kind & ASL_Sufkind_real) {
		rp = d->u.r;
		if (map) {
			for(i = 0; i < n; ++i)
				if (rp[i] && map[i] >= 0)
					++nz;
			}
		else {
			rpe = rp + n;
			while(rp < rpe)
				if (*rp++)
					++nz;
			}
		}
	else {
		ip = d->u.i;
		if (map) {
			for(i = 0; i < n; ++i)
				if (ip[i] && map[i] >= 0)
					++nz;
			}
		else {
			ipe = ip + n;
			while(ip < ipe)
				if (*ip++)
					++nz;
			}
		}
	sh->n = nz;
	sh->namelen = strlen(d->sufname) + 1;
	sh->tablen = 0;
	if (d->table)
		sh->tablen = strlen(d->table) + 1;
	}

 static long
tablines(const char *s)
{
	long n;
	if (!s)
		return 0;
	n = 1;
	while(*s)
		if (*s++ == '\n')
			n++;
	return n;
	}

 static void
showsol(ASL *asl, real *x, int n, int n0, Name name, const char *what, const char *pfix)
{
	int i, j, k, k0;

	if (!x || n <= 0)
		return;
	k0 = k = strlen(what);
	for(i = 0; i < n0; ++i)
		if ((j = strlen((*name)(asl,i,0))) > k)
			k = j;
	k += 2;
	printf("\n%s%*s%svalue\n", what, k-k0, "", pfix);
	for(i = 0; i < n0; ++i)
		printf("%-*s%.g\n", k, (*name)(asl,i,0), x[i]);
	}

 static real *
scale(real *x, real *s, real **yp, int n)
{
	real *xe, *y, *y0;

	y0 = y = *yp;
	xe = x + n;
	while(x < xe)
		*y++ = *s++ * *x++;
	*yp = y;
	return y0;
	}

 static real*
copy(int n, int n1, real *x, real **yp, int *z, int *zap)
{
	int i, j;
	real *y;

	y = *yp;
	*yp = y + n;
	for(i = 0; i < n1; ++i)
		if ((j = z[i]) >= 0 && j < n)
			y[j] = x[i];
	if (zap)
		for(n = zap[0], i = 1; i <= n; ++i)
			y[zap[i]] = 0.;
	return y;
	}

 static void
equ_adjust1(int *ip, real *LU, real *U, int n)
{
	int i = 0;
	if (U) {
		for(; i < n; i++)
			if (LU[i] == U[i] && (ip[i] == 3 || ip[i] == 4))
				ip[i] = 5;
		}
	else if (LU)
		for(; i < n; i++, LU += 2)
			if (LU[0] == LU[1] && (ip[i] == 3 || ip[i] == 4))
				ip[i] = 5;
	}

 void
equ_adjust_ASL(ASL *asl, int *cstat, int *rstat)
{
	if (cstat)
		equ_adjust1(cstat, LUv, Uvx, n_var);
	if (rstat)
		equ_adjust1(rstat, LUrhs, Urhsx, n_con);
	}

 static long
AMPL_version_ASL(ASL *asl)
{
	char *s;
	if (ampl_options[0] >= 5)
		return ampl_options[5];
	if (!(s = getenv("version")))
		return 0;
	for(;;) {
		switch(*s++) {
		 case 0: return 0;
		 case 'V': if (!strncmp(s,"ersion ", 7))
				goto break2;
		 }
		}
 break2:
	return strtol(s+7,0,10);
	}

 int
write_solfx_ASL(ASL *asl, const char *msg, double *x, double *y, Option_Info *oi,
		Fwrite fw_d, Fwrite fw_i, Fwrite fw_s, const char *solfname)
{
	FILE *f;
	SufDesc *d;
	SufHead sh;
	char *bsmsg, buf[80];
	const char *s, *s1, *s2;
	fint J[2], m, z[4];
	ftnlen L[6];
	int N, binary, i, i1, *ip, j, k, n, nlneed, rv, tail, wantsol, *zz;
	real *rp, *x0, *y0, *y1, *xycopy;
	size_t nn;
	static const char *wkind[] = {"w", "wb"};

	if (!asl || asl->i.ASLtype < 1 || asl->i.ASLtype > 5)
		badasl_ASL(asl,0,"write_sol");

	rv = 0;
	bsmsg = 0;
	if ((nlneed = need_nl) > 0) {
		if (oi && oi->bsname && (i = nlneed-2) > 0
		 && amplflag
		 && !strncmp(msg,oi->bsname,i)
		 && !strncmp(msg+i,": ",2)
		 && AMPL_version_ASL(asl) >= 20020401L) {
			bsmsg = (char*)Malloc(nlneed + strlen(msg) + 1);
			memset(bsmsg, '\b', nlneed);
			strcpy(bsmsg+nlneed, msg);
			msg = bsmsg;
			nlneed = 0;
			}
		}
	xycopy = 0;
	wantsol = 1;
	if (!solfname)
		solfname = asl->i.solfile;
	if (oi && !solfname)
		wantsol = oi->wantsol;
	if (wantsol || amplflag) {
		k = 0;
		y1 = 0;
		if ((x0 = x)) {
			if (asl->i.vmap) {
				k = asl->i.n_var0;
				if (asl->i.vscale)
					k += n_var;
				}
			else if (asl->i.n_var0 > n_var)
				k = asl->i.n_var0;
			else if (asl->i.vscale)
				k = n_var;
			}
		if (y) {
			if (asl->i.cmap) {
				k += asl->i.n_con0;
				if (asl->i.lscale)
					k += n_con;
				}
			else if (asl->i.n_con0 > n_con)
				k += asl->i.n_con0;
			else if (asl->i.lscale)
				k += n_con;
			}
		if (k)
			y1 = xycopy = (real*)Malloc(k*sizeof(real));
		if (x) {
			if (asl->i.vscale)
				x = scale(x, asl->i.vscale, &y1, n_var);
			if ((ip = asl->i.vmap))
				x = copy(asl->i.n_var0, n_var, x, &y1, ip, asl->i.vzap);
			else if (x == x0 && asl->i.n_var0 > n_var) {
				memcpy(y1, x, n_var*sizeof(real));
				x = y1;
				y1 += asl->i.n_var0;
				}
			}
		z[0] = m = asl->i.n_con0;
		if (!y)
			m = 0;
		else {
			y0 = y;
			if (asl->i.lscale)
				y = scale(y, asl->i.lscale, &y1, n_con);
			if ((ip = asl->i.cmap))
				y = copy(asl->i.n_con0, n_con, y, &y1, ip, asl->i.czap);
			else if (y0 == y && asl->i.n_con0 > n_con) {
				memcpy(y1, y, n_con*sizeof(real));
				y = y1;
				}
			}
		if (asl->i.Or && x)
			obj_adj_xy_ASL(asl, x, x0, y);
		}
	if (!amplflag && !(wantsol & 1))
		goto write_done;
	tail = 0;
	if (obj_no || solve_code != -1)
		tail = 1;
	else  {
		for(i1 = 0; i1 < 4; i1++)
		    for(d = asl->i.suffixes[i1]; d; d = d->next)
			if (d->kind & ASL_Sufkind_output
			 && (d->kind & ASL_Sufkind_real
					? (int*)d->u.r : d->u.i)) {
				tail = 1;
				goto break2;
				}
		}
 break2:
	binary = binary_nl & 1;
	if (!solfname) {
		strcpy(stub_end, ".sol");
		solfname = filename;
		}
	f = fopen(solfname, wkind[binary]);
	if (!f) {
		fprintf(Stderr, "can't open %s\n", solfname);
		rv = 1;
		goto ret;
		}
	z[1] = m;
	z[2] = n = asl->i.n_var0;
	if (!x)
		n = 0;
	z[3] = n;
	k = (int)ampl_options[0];
	if (binary) {
		L[0] = 6;
		L[1] = strlen(msg);
		L[2] = 0;
		L[3] = (ampl_options[0] + 5)*sizeof(fint) + 7;
		L[4] = m*sizeof(double);
		L[5] = n*sizeof(double);
		(*fw_i)(L, sizeof(ftnlen), 1, f);
		fwrite("binary", 6, 1, f);
		(*fw_i)(L, sizeof(ftnlen), 2, f);
		if (L[1]) {
			fwrite(msg, L[1], 1, f);
			(*fw_i)(L+1, sizeof(ftnlen), 2, f);
			}
		if (k) {
			(*fw_i)(L+2, sizeof(ftnlen), 2, f);
			fwrite("Options",7,1,f);
			nn = (size_t)ampl_options[0]+1;
			if (ampl_options[2] == 3)
				ampl_options[0] += 2;
			(*fw_i)(ampl_options, sizeof(fint), nn, f);
			(*fw_i)(z, sizeof(fint), 4, f);
			if (ampl_options[2] == 3)
				(*fw_d)(&ampl_vbtol, sizeof(real), 1, f);
			(*fw_i)(L+3, sizeof(ftnlen), 2, f);
			}
		else {
			(*fw_i)(L+2, sizeof(ftnlen), 1, f);
			(*fw_i)(L+4, sizeof(ftnlen), 1, f);
			}
		if (y)
			(*fw_d)(y, sizeof(double), m, f);
		(*fw_i)(L+4, sizeof(ftnlen), 2, f);
		if (x)
			(*fw_d)(x, sizeof(double), n, f);
		(*fw_i)(L+5, sizeof(ftnlen), 1, f);
		if (tail)
		  switch(asl->i.flags & 1) {
		    case 0:
			if (obj_no) {
				L[0] = L[2] = sizeof(fint);
				L[1] = obj_no;
				(*fw_i)(L, sizeof(fint), 3, f);
				}
			break;
		    case 1:
			L[0] = L[3] = 2*sizeof(fint);
			L[1] = obj_no;
			L[2] = solve_code;
			(*fw_i)(L, sizeof(fint), 4, f);
			for(i1 = 0; i1 < 4; i1++)
			  for(d = asl->i.suffixes[i1]; d; d = d->next)
			    if (d->kind & ASL_Sufkind_output
			     && (d->kind & ASL_Sufkind_real
					? (int*)d->u.r : d->u.i)) {
				getsufhead(asl, d, &sh, &N, &zz);
				L[0] = sizeof(sh) + sh.namelen + sh.tablen
					+ sh.n*(sizeof(int) +
						(d->kind & ASL_Sufkind_real
						? sizeof(real) : sizeof(int)));
				(*fw_i)(L, sizeof(fint), 1, f);
				(*fw_s)(&sh, sizeof(sh), 1, f);
				fwrite(d->sufname, sh.namelen, 1, f);
				if (sh.tablen)
					fwrite(d->table, sh.tablen, 1, f);
				i = j = 0;
				if (d->kind & ASL_Sufkind_real)
				    for(rp = d->u.r; i < N; i++) {
					if (rp[i]) {
						if (zz) {
							if ((J[0] = zz[i]) < 0)
								continue;
							}
						else
							J[0] = i;
						(*fw_i)(J, sizeof(fint), 1, f);
						(*fw_d)(rp+i,sizeof(real),1,f);
						}
					}
				else
				    for(ip = d->u.i; i < N; i++) {
					if ((J[1] = ip[i])) {
						if (zz) {
							if ((J[0] = zz[i]) < 0)
								continue;
							}
						else
							J[0] = i;
						(*fw_i)(J, sizeof(fint), 2, f);
						}
					}
				(*fw_i)(L, sizeof(fint), 1, f);
				}
			}
		}
	else {
		if (*(s = msg)) {
			for(s2 = s + strlen(s); s2 > s && s2[-1] == '\n'; --s2);
			while (s < s2) {
				for(s1 = s; *s1 != '\n' && ++s1 < s2;);
				fprintf(f, s == s1 ? " \n" : "%.*s\n",s1-s,s);
				s = s1 + 1;
				}
			}
		fprintf(f, "\n");
		if ((k = (int)ampl_options[0])) {
			if (ampl_options[2] == 3)
				ampl_options[0] += 2;
			fprintf(f, "Options\n");
			for(i = 0; i <= k; i++)
				fprintf(f,"%ld\n",(long)ampl_options[i]);
			fprintf(f,"%ld\n%ld\n%ld\n%ld\n",
				(long)z[0],(long)z[1],(long)z[2],(long)z[3]);
			if (ampl_options[2] == 3) {
				g_fmtp(buf, ampl_vbtol, 0);
				fprintf(f, "%s\n", buf);
				}
			}
		y1 = y;
		while(--m >= 0) {
			g_fmtp(buf, *y1++, 0);
			fprintf(f,"%s\n", buf);
			}
		y1 = x;
		while(--n >= 0) {
			g_fmtp(buf, *y1++, 0);
			fprintf(f, "%s\n", buf);
			}
		if (tail)
		  switch(asl->i.flags & 1) {
		    case 0:
			if (obj_no)
				fprintf(f, "objno %d\n", obj_no);
			break;
		    case 1:
			fprintf(f, "objno %d %d\n", obj_no, solve_code);
			for(i1 = 0; i1 < 4; i1++)
			  for(d = asl->i.suffixes[i1]; d; d = d->next)
			    if (d->kind & ASL_Sufkind_output
			     && (d->kind & ASL_Sufkind_real
					? (int*)d->u.r : d->u.i)) {
				getsufhead(asl, d, &sh, &N, &zz);
				fprintf(f, "suffix %ld %ld %ld %ld %ld\n%s\n",
					(long)sh.kind, (long)sh.n,
					(long)sh.namelen, (long)sh.tablen,
					tablines(d->table), d->sufname);
				if (sh.tablen)
					fprintf(f, "%s\n", d->table);
				i = j = 0;
				if (d->kind & ASL_Sufkind_real)
				    for(rp = d->u.r; i < N; i++) {
					if (rp[i]) {
						if (zz) {
							if ((j = zz[i]) < 0)
								continue;
							}
						else
							j = i;
						fprintf(f, "%d %.g\n",
							j, rp[i]);
						}
					}
				else
				    for(ip = d->u.i; i < N; i++) {
					if (ip[i]) {
						if (zz) {
							if ((j = zz[i]) < 0)
								continue;
							}
						else
							j = i;
						fprintf(f, "%d %d\n",
							j, ip[i]);
						}
					}
				}
			}
		}
	fclose(f);
 write_done:
	if ((i = nlneed)) {
		if (i > sizeof(buf)-1 || i < 0)
			printf("\n");
		else {
			buf[i] = 0;
			do buf[--i] = '\b';
				while(i > 0);
			printf(buf);
			}
		}
	if (!amplflag) {
		if (!oi || !(oi->wantsol & 8))
			printf("%s\n", msg);
		if (wantsol & 2)
			showsol(asl, x, n_var, asl->i.n_var0,
				var_name_nomap_ASL, "variable", "");
		if (wantsol & 4)
			showsol(asl, y, n_con, asl->i.n_con0,
				con_name_nomap_ASL, "constraint", "dual ");
		}
 ret:
	if (xycopy)
		free(xycopy);
	if (bsmsg)
		free(bsmsg);
	return rv;
	}

 void
write_sol_ASL(ASL *asl, const char *msg, double *x, double *y, Option_Info *oi)
{
	if (write_solfx_ASL(asl, msg, x, y, oi, fwrite, fwrite, fwrite, 0))
		exit(2);
	}

 int
write_solf_ASL(ASL *asl, const char *msg, double *x, double *y, Option_Info *oi, const char *fname)
{
	return write_solfx_ASL(asl, msg, x, y, oi, fwrite, fwrite, fwrite, fname);
	}

/* Affected by ASL update of 20020503 */
