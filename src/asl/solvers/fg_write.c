/****************************************************************
Copyright (C) 1999,2000 Lucent Technologies
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

#ifdef __cplusplus
extern "C" {
#endif

typedef int Pf ANSI((FILE*, const char*, ...));

 typedef struct Staticfgw
{
	Pf *pf_;
	FILE *nl_;
	efunc **r_ops_;
	jmp_buf wjb;
	cexp1 *cexps1_;
	expr_v *v;
	int com1off;
	int nv0;
	} Staticfgw;

#undef r_ops
#define r_ops S->r_ops_

#include "r_opn.hd"

 static int
aprintf(FILE *fd, const char *fmt, ...)
{
	char *s;
	char buf[32];
	va_list ap;
	int i, j;
	double x;
	int rc = 0;

	va_start(ap, fmt);
	if (*fmt != '%')
		rc++;
	for(;;) {
		for(;;) {
			switch(i = *fmt++) {
				case 0:	  goto done;
				case '%': break;
				default:  putc(i,fd);
					  continue;
				}
			break;
			}
		rc++;
		switch(*fmt++) {
			case 'c':
				i = va_arg(ap, int);
				putc(i,fd);
				continue;
			case 'd':
				i = va_arg(ap, int);
				if (i < 0) {
					putc('-',fd);
					i = -i;
					}
				s = buf;
				do {
					j = i / 10;
					*s++ = i - 10*j + '0';
					}
					while((i = j));
				do {
					i = *--s;
					putc(i,fd);
					}
					while(s > buf);
				continue;
			case '.':
				while(*fmt++ != 'g');
			case 'g':
				x = va_arg(ap, double);
				i = g_fmt(s = buf, x);
				goto have_s;
			case 's':
				s = va_arg(ap, char*);
 have_s:
				while((i = *s++))
					putc(i,fd);
				continue;
			default:
				fprintf(Stderr, "aprintf bug: unexpect fmt %s\n",
					fmt-1);
				exit(1);
			}
		}
 done:
	va_end(ap);
	return rc;
	}

#ifndef Int
#define Int Long
#endif

 static int
bprintf(FILE *fd, const char *fmt, ...)
{
	char *s;
	int i, rc;
	size_t len;
	union U { double x; short sh; Long L; Int i; char c; } u;
	va_list ap;

	va_start(ap, fmt);

	rc = 0;
	len = 0; /* silence buggy "not-initialized" warning */
	if ((i = *fmt) != '%') {
		fmt++;
#ifdef DMG
		if (i != 'i')
#endif
		{
		u.c = i;
		fwrite(&u.c, 1, 1, fd);
		rc++;
		}}

	for(;;) {
		while(*fmt == ' ')
			fmt++;
		if (*fmt++ != '%')
			break;
		switch(*fmt++) {
			case 'c':
				u.c = va_arg(ap, int);
				len = 1;
				break;
			case 'd':
				u.i = va_arg(ap, int);
				len = sizeof(Int);
				break;
			case '.':
				while(*fmt++ != 'g');
			case 'g':
				u.x = va_arg(ap, double);
				len = sizeof(double);
				break;
			case 'h':
				u.sh = va_arg(ap, int);
				len = sizeof(short);
				if (*fmt == 'd')
					fmt++;
				break;
			case 'l':
				u.L = (Long)va_arg(ap, long);
				len = sizeof(Long);
				if (*fmt == 'd')
					fmt++;
				break;
			case 's':
				s = va_arg(ap, char*);
				u.i = strlen(s);
				fwrite((char *)&u.i, sizeof(Int), 1, fd);
				fwrite(s, u.i, 1, fd);
				goto s_written;
			default:
				fprintf(Stderr, "bprintf bug: unexpect fmt %s\n",
					fmt-1);
				exit(1);
			}
		fwrite((char *)&u.L, len, 1, fd);
 s_written:
		rc++;
		}
	va_end(ap);
	return rc;
	}

#define pf S->pf_
#define nl S->nl_

 static void
eput(Staticfgw *S, expr *e)
{
	arglist *al;
	de *d, *dee;
	efunc *op;
	expr **ap, **ape, etemp;
	expr_if *eif;
	expr_f *ef;
	expr_h *eh;
	expr_n *en, entemp;
	expr_v *v;
	expr_va *va;
	func_info *fi;
	int i, nop;
	plterm *p;
	real *r, *re;

 top:
	op = e->op;
	if (op == f_OP1POW || op == f_OPCPOW)
		op = f_OPPOW;
	else if (op == f_OP2POW) {
		op = f_OPPOW;
		etemp.L.e = e->L.e;
		etemp.R.en = &entemp;
		entemp.op = (efunc_n*)f_OPNUM;
		entemp.v = 2.;
		e = &etemp;
		}
	nop = Intcast op;
	if ((i = optypeb[nop]) < 7)
		(*pf)(nl, "o%d\n", nop = Intcast op);
	switch(i) {
	 case 1:
		e = e->L.e;
		goto top;
	 case 2:
		eput(S, e->L.e);
		e = e->R.e;
		goto top;
	 case 3:
		va = (expr_va *)e;
		d = dee = va->L.d;
		while(dee->e)
			dee++;
		(*pf)(nl, "%d\n", (int)(dee-d));
		for(; d < dee; d++)
			eput(S, d->e);
		break;
	 case 4:
		p = e->L.p;
		(*pf)(nl, "%d\n", p->n);
		r = p->bs;
		re = r + 2*p->n - 1;
		while(r < re)
			(*pf)(nl, "n%g\n", *r++);
		e = e->R.e;
		goto top;
	 case 5:
		eif = (expr_if*)e;
		eput(S, eif->e);
		eput(S, eif->T);
		e = eif->F;
		goto top;
	 case 6:
		ap = e->L.ep;
		ape = e->R.ep;
		(*pf)(nl, "%d\n", (int)(ape - ap));
		while(ap < ape)
			eput(S, *ap++);
		break;
	 case 7:
		ef = (expr_f*)e;
		fi = ef->fi;
		al = ef->al;
		(*pf)(nl, "f%d %d\n", fi->findex, al->n);
		ap = ef->args;
		ape = ap + al->n;
		while(ap < ape)
			eput(S, *ap++);
		break;
	 case 8:
		eh = (expr_h*)e;
		(*pf)(nl, "h%d:%d\n", strlen(eh->sym), eh->sym);
		break;
	 case 9:
		en = (expr_n*)e;
		(*pf)(nl, "n%g\n", en->v);
		break;
	 case 10:
		v = (expr_v*)e;
		i = (int)(v - S->v);
		(*pf)(nl, "v%d\n", i);
		break;
	 default:
		fprintf(Stderr, "fg_write: unexpected type %d in eput.\n",
			optypeb[nop]);
		longjmp(S->wjb, 1);
	 }
	}

#define offset_of(t,c) ((size_t)(char *)&((t*)0)->c)

 static void
coput(Staticfgw *S, int c, cde *de, int n, int *cexp1st, char *ot, int voff,
	int nn, real *oc, char *Not)
{
	cexp1 *ce;
	expr_v *v, *v0;
	int i, i1, j, je, k;
	linpart *L, *Le;
	real t;

	if (cexp1st) {
		j = cexp1st[0];
		ce = S->cexps1_ + j;
		v0 = S->v;
		}
	else /* silence buggy "not-initialized" warnings */
		{ ce = 0; v0 = 0; j = 0; }

	for(i = 0; i < n; i++) {
		if (cexp1st) {
			je = cexp1st[i1 = i + 1];
			i1 += voff;
			while(j < je) {
				k = S->com1off + j++;
				(*pf)(nl, "V%d %d %d\n", k, ce->nlin, i1);
				L = ce->L;
				for(Le = L + ce->nlin; L < Le; L++) {
					v = (expr_v*)((char*)L->v.rp
							- offset_of(expr_v,v));
					k = (int)(v - v0);
					(*pf)(nl, "%d %g\n", k, L->fac);
					}
				eput(S, ce->e);
				ce++;
				}
			}
		if (ot)
			(*pf)(nl, "%c%d %d\n", c, i, ot[i]);
		else
			(*pf)(nl, "%c%d\n", c, i);
		eput(S, de[i].e);
		}
	t = 0.;
	for(n += nn; i < n; i++) {
		if (ot) {
			(*pf)(nl, "%c%d %d\n", c, i, Not ? *Not++ : 0);
			if (oc)
				t = *oc++;
			}
		else
			(*pf)(nl, "%c%d\n", c, i);
		(*pf)(nl, "n%g\n", t);
		}
	}

#undef pf
#undef nl

 static void
iguess(Pf *pf, FILE *nl, int c, real *x, char *havex, int n, int nn, real *y)
{
	int i, k;

	if (n + nn <= 0)
		return;
	i = k = 0;
	if (x) {
		if (havex) {
			while(i < n)
				if (havex[i++])
					k++;
			}
		else {
			while(i < n)
				if (x[i++])
					k++;
			}
		}
	if (y)
		for(i = 0; i < nn; i++)
			if (y[i])
				k++;
	if (!k)
		return;
	(*pf)(nl, "%c%d\n", c, k);
	if (x) {
		if (havex) {
			for(i = 0; i < n; i++)
				if (havex[i])
					(*pf)(nl, "%d %g\n", i, x[i]);
			}
		else {
			for(i = 0; i < n; i++)
				if (x[i])
					(*pf)(nl, "%d %g\n", i, x[i]);
			}
		}
	if (y) {
		for(i = 0; i < nn; i++)
			if (y[i])
				(*pf)(nl, "%d %g\n", i+n, y[i]);
		}
	}

 static void
br(Pf *pf, FILE *nl, int c, real *Lb, real *Ub, int n)
{
	int i;
	real L, U;

	if (n <= 0)
		return;
	if (c)
		(*pf)(nl, "%c\n", c);
	for(i = 0; i < n; i++) {
		L = *Lb++;
		U = Ub ? *Ub++ : *Lb++;
		if (L <= negInfinity)
			(*pf)(nl, U >= Infinity ? "3\n" : "1 %g\n", U);
		else
			(*pf)(nl, U >= Infinity ? "2 %g\n"
				: L == U ? "4 %g\n"
					 : "0 %g %g\n",
				L, U);
		}
	}

 static void
Gput(Pf *pf, FILE *nl, int c, int i, int n, ograd **ogp)
{
	ograd *og;
	int k;

	if (n <= 0)
		return;
	for(n += i; i < n; i++, ogp++) {
		if (!(og = *ogp))
			continue;
		k = 0;
		do k++;
			while((og = og->next));
		(*pf)(nl, "%c%d %d\n", c, i, k);
		for(og = *ogp; og; og = og->next)
			(*pf)(nl, "%d %g\n", og->varno, og->coef);
		}
	}

 static void
k2put(Pf *pf, FILE *nl, cgrad **cgp, int nc, int n, int k, int nnv,
	int nnc, ograd **ogp)
{
	cgrad *cg;
	ograd *og;
	int i, n1, *z;

	if (k) {
		n1 = n + nnv;
		z = (int*)Malloc(n1*sizeof(int));
		memset(z, 0, n1*sizeof(int));
		for(i = 0; i < nc; i++)
			for(cg = cgp[i]; cg; cg = cg->next)
				z[cg->varno]++;
		for(i = 0; i < nnc; i++)
			for(og = ogp[i]; og; og = og->next)
				z[og->varno]++;
		(*pf)(nl, "k%d\n", --n1);
		for(i = k = 0; i < n1; i++)
			(*pf)(nl, "%d\n", k += z[i]);
		free(z);
		}
	for(i = 0; i < nc; i++) {
		if (!(cg = cgp[i]))
			continue;
		k = 0;
		do k++;
			while((cg = cg->next));
		(*pf)(nl, "J%d %d\n", i, k);
		for(cg = cgp[i]; cg; cg = cg->next)
			(*pf)(nl, "%d %g\n", cg->varno, cg->coef);
		}
	Gput(pf, nl, 'J', nc, nnc, ogp);
	}

 static void
k1put(Pf *pf, FILE *nl, int *cs, real *a, int *rn, int nc, int n,
	int nnv, int nnc, ograd **ogp)
{
	int *cs1, i, j, j1, k, ftn, nz;
	cgrad *cg, *cg0, *cg1, **cgp, **cgq;
	ograd *og;

	ftn = cs[0];
	nz = cs[n] - ftn;
	k = n;
	if (nnc) {
		k += nnv;
		if (nz <= k)
			nz = k + 1;
		}
	cg0 = cg1 = (cgrad*)Malloc(nz*sizeof(cgrad) + nc*sizeof(cgrad*));
	cs1 = cs;
	if (nnc) {
		cs1 = (int*)cg1;
		for(i = 0; i < n; i++)
			cs1[i] = cs[i+1] - cs[i];
		while(i < k)
			cs1[i++] = 0;
		for(i = 0; i < nnc; i++)
			for(og = ogp[i]; og; og = og->next)
				cs1[og->varno]++;
		j = ftn;
		for(i = 0; i < k; i++) {
			j1 = j + cs1[i];
			cs1[i] = j;
			j = j1;
			}
		cs1[k] = j;
		}
	(*pf)(nl, "k%d\n", k - 1);
	for(i = 1; i < k; i++)
		(*pf)(nl, "%d\n", cs1[i] - ftn);
	memset(cgp = (cgrad**)(cg0 + nz), 0, nc*sizeof(cgrad*));
	k = cs[n] - ftn;
	for(i = n; --i >= 0;) {
		for(j = cs[i] - ftn; --k >= j;) {
			cg = cg1++;
			cg->coef = a[k];
			cgq = cgp + (cg->varno = rn[k] - ftn);
			cg->next = *cgq;
			*cgq = cg;
			}
		}
	k2put(pf, nl, cgp, nc, n, 0, nnv, nnc, ogp);
	free(cg0);
	}

 static int
LUcheck(int n, real *LU, real *u, int *nnep, int *nnrp)
{
	int i, nne, nnr;
	real L, U;

	nnr = nne = 0;
	if (!LU)
		return 1;
	for(i = 0; i < n; i++) {
		L = *LU++;
		U = u ? *u++ : *LU++;
		if (L < U) {
			if (L > negInfinity && U < Infinity)
				nnr++;
			}
		else if (U <= negInfinity
		 || L >= Infinity
		 || L > U
		 || L != L /* NaN */
		 || U != U)
			return 1;
		else
			nne++;
		}
	if (nnep) {
		*nnep = nne;
		*nnrp = nnr;
		}
	return 0;
	}

 static int
ogcheck(int n, int nn, ograd **ogp, int *nzp)
{
	int nz;
	ograd *og;

	if (!ogp)
		return 1;
	nz = 0;
	n += nn;
	while(nn--)
		for(og = *ogp++; og; og = og->next) {
			++nz;
			if (og->varno < 0
			 || og->varno >= n
			 || og->coef != og->coef
			 || og->coef == Infinity
			 || og->coef == negInfinity)
				return 1;
			}
	*nzp = nz;
	return 0;
	}

 static SufDesc*
reverse(SufDesc *sd)
{
	SufDesc *sn, *sp;
	sp = 0;
	while(sd) {
		sn = sd->next;
		sd->next = sp;
		sp = sd;
		sd = sn;
		}
	return sp;
	}

 int
fg_write_ASL(ASL *a, const char *stub, NewVCO *nu, int flags)
{
	ASL_fg *asl = (ASL_fg*)a;
	FILE *nl;
	Pf *pf;
	Staticfgw S;
	SufDesc *sd, *sd0;
	cexp *ce, *cee;
	char buf[256], *nbuf, *ts;
	const char *eol, *name, *obase, *s;
	efunc *rops[N_OPS];
	expr_v *v;
	func_info *fi;
	int ak, c, i, j, *ip, *ipe, n, nnc, nne, nno, nnr, nnv, nnzc, nnzo;
	int nx, oblen, rflag;
	linpart *L, *Le;
	real *r, *re, t;
	static NewVCO nu0;

	ASL_CHECK(a, ASL_read_fg, "fg_write");
	if ((comc1 && !c_cexp1st) || (como1 && !o_cexp1st))
		return ASL_writeerr_badcexp1st;
	nnc = nne = nno = nnr = nnv = nnzc = nnzo = 0;
	if (!nu || (nu->nnv == 0 && nu->nnc == 0 && nu->nno == 0))
		nu = &nu0;
	else {
		nnc = nu->nnc;
		nno = nu->nno;
		nnv = nu->nnv;
		if ((nnv <= 0
		  || nnc < 0
		  || nno < 0
		  || nnc + nno <= 0
		  || nnc > 0) && !nu->LUnc)
			return ASL_writeerr_badNewVCO;
		if (LUcheck(nnv, nu->LUnv, nu->Unv, 0, 0))
			return ASL_writeerr_badNewVCO;
		n = n_var + nnv;
		if (nnc) {
			if (LUcheck(nnc, nu->LUnc, nu->Unc, &nnr, &nne))
				return ASL_writeerr_badNewVCO;
			if (ogcheck(n, nnc, nu->newc, &nnzc))
				return ASL_writeerr_badNewVCO;
			}
		if (nno) {
			if (ogcheck(n, nno, nu->newo, &nnzo))
				return ASL_writeerr_badNewVCO;
			if ((s = nu->ot))
			    for(i = 0; i < nno; i++)
				if (s[i] & ~1)
					return ASL_writeerr_badNewVCO;
			if ((r = nu->oc))
			    for(re = r + nno; r < re; r++) {
				if ((t = *r) <= negInfinity
				 || t >= Infinity
				 || t != t)
					return ASL_writeerr_badNewVCO;
				}
			}
		}

	S.r_ops_ = rops;
	for(i = 0; i < N_OPS; i++)
		rops[i] = (efunc*)(size_t)i;

	s = name = obase = stub;
	while(*s)
	  switch(*s++) {
		case '/':
		case '\\':
			obase = s;
		}
	c = s - stub;
	nbuf = 0;
	oblen = s - obase;
	if (c <= 3 || strcmp(s - 3, ".nl")) {
		ts = buf;
		if (c + 4 > sizeof(buf))
			ts = nbuf = (char*)Malloc(c+4);
		memcpy(ts, stub, c);
		strcpy(ts+c, ".nl");
		name = ts;
		}
	else
		oblen -= 3;
	nl = fopen(name, "wb");
	if (nbuf)
		free(nbuf);
	if (!nl)
		return ASL_writeerr_openfail;
	i = setjmp(S.wjb);
	if (i) {
		fclose(nl);
		return ASL_writeerr_badrops;
		}
	if (flags & ASL_write_ASCII) {
		ak = 0;
		c = 'g';
		pf = aprintf;
		}
	else {
		ak = Arith_Kind_ASL;
		c = 'b';
		pf = bprintf;
		}
	S.nl_ = nl;
	S.pf_ = pf;
	eol = (char*)(flags & ASL_write_CR ? "\r\n" : "\n");
	fprintf(nl, "%c%d", c, n = ampl_options[0]);
	for(i = 1; i <= n; i++)
		fprintf(nl, " %d", ampl_options[i]);
	if (ampl_options[2] == 3)
		fprintf(nl, " %.g", ampl_vbtol);
	fprintf(nl, "\t# problem %.*s%s", oblen, obase, eol);
	fprintf(nl, " %d %d %d %d", n_var + nnv, n_con + nnc,
		n_obj + nno, nranges + nnr);
	s = "";
	if ((n = n_eqn + nne) >= 0) {
		fprintf(nl, " %d", n);
		s = ", eqns";
		}
	fprintf(nl, "\t# vars, constraints, objectives, ranges%s%s", s, eol);
	if (n_cc | nlcc)
		fprintf(nl, " %d %d %d %d%s%s", nlc, nlo, n_cc, nlcc,
		"\t# nonlinear constrs, objs; ccons: lin, nonlin", eol);
	else
		fprintf(nl, " %d %d\t# nonlinear constraints, objectives%s",
			nlc, nlo, eol);
	fprintf(nl, " %d %d\t# network constraints: nonlinear, linear%s",
		nlnc, lnc, eol);
	fprintf(nl, " %d %d %d%s%s", nlvc, nlvo, nlvb,
		"\t# nonlinear vars in constraints, objectives, both", eol);
	s = "";
	fprintf(nl, " %d %d", nwv, nfunc);
	if (ak | asl->i.flags) {
		fprintf(nl, " %d %d", ak, asl->i.flags);
		s = "; arith, flags";
		}
	fprintf(nl, "\t# linear network variables; functions%s%s", s, eol);
	fprintf(nl, " %d %d %d %d %d%s%s", nbv, niv, nlvbi, nlvci, nlvoi,
		"\t# discrete variables: binary, integer, nonlinear (b,c,o)",
		eol);
	fprintf(nl, " %d %d\t# nonzeros in Jacobian, gradients%s",
		nzc + nnzc, nzo + nnzo, eol);
	fprintf(nl, " 0 0\t# max name lengths: constraints, variables%s", eol);
	fprintf(nl, " %d %d %d %d %d\t# common exprs: b,c,o,c1,o1%s",
		comb, comc, como, comc1, como1, eol);

	for(i = 0; i < nfunc; i++) {
		fi = funcs[i];
		fi->findex = i; /* for eput */
		(*pf)(nl, "F%d %d %d %s\n", i, fi->ftype, fi->nargs, fi->name);
		}

	for(i = 0; i < 4; i++) {
		if (!(sd = asl->i.suffixes[i]))
			continue;
		nx = (&asl->i.n_var_)[i];
		for(sd = sd0 = reverse(sd); sd; sd = sd->next) {
			n = rflag = 0;
			if (sd->kind & ASL_Sufkind_real) {
				rflag = ASL_Sufkind_real;
				r = sd->u.r;
				re = r + nx;
				while(r < re)
					if (*r++)
						n++;
				}
			else {
				ip = sd->u.i;
				ipe = ip + nx;
				while(ip < ipe)
					if (*ip++)
						n++;
				}
			if (!n)
				continue;
			(*pf)(nl, "S%d %d %s\n", i | rflag, n, sd->sufname);
			j = 0;
			if (rflag) {
				r = sd->u.r;
				for(; j < nx; j++)
					if (r[j])
						(*pf)(nl, "%d %g\n", j, r[j]);
				}
			else {
				ip = sd->u.i;
				for(; j < nx; j++)
					if (ip[j])
						(*pf)(nl, "%d %d\n", j, ip[j]);
				}
			}
		reverse(sd0);
		}
	ce = cexps;
	n = n_var + nnv;
	S.v = var_e;
	for(cee = ce + comb + comc + como; ce < cee; ce++) {
		(*pf)(nl, "V%d %d %d\n", n++, ce->nlin, 0);
		L = ce->L;
		for(Le = L + ce->nlin; L < Le; L++) {
			v = (expr_v*)((char*)L->v.rp - offset_of(expr_v,v));
			(*pf)(nl, "%d %g\n", (int)(v - S.v), L->fac);
			}
		eput(&S, ce->e);
		}
	S.cexps1_ = asl->I.cexps1_;
	S.nv0 = n_var;
	S.com1off = S.nv0 + comb + comc + como;
	coput(&S, 'C', con_de, n_con, c_cexp1st, 0, 0, nnc, 0, 0);
	coput(&S, 'O', obj_de, n_obj, o_cexp1st, objtype, n_con,
		nno, nu->oc, nu->ot);
	iguess(pf, nl, 'd', pi0, havepi0, n_con, nnc, nu->d0);
	iguess(pf, nl, 'x', X0, havex0, n_var, nnv, nu->x0);
	br(pf, nl, 'r', LUrhs, Urhsx, n_con);
	br(pf, nl, 0, nu->LUnc, nu->Unc, nnc);
	br(pf, nl, 'b', LUv, Uvx, n_var);
	br(pf, nl, 0, nu->LUnv, nu->Unv, nnv);
	if (A_vals)
		k1put(pf, nl, A_colstarts, A_vals, A_rownos, n_con, n_var,
			nnv, nnc, nu->newc);
	else
		k2put(pf, nl, Cgrad, n_con, n_var, 1, nnv, nnc, nu->newc);

	Gput(pf, nl, 'G', 0, n_obj, Ograd);
	Gput(pf, nl, 'G', n_obj, nno, nu->newo);

	fclose(nl);
	return 0;
	}

 int
fg_wread_ASL(ASL *asl, FILE *f, int flags)
{
	want_xpi0 = 7;
	if (comc1)
		c_cexp1st = (int*)M1zapalloc((n_con + 1)*sizeof(int));
	if (como1)
		o_cexp1st = (int*)M1zapalloc((n_obj + 1)*sizeof(int));
	if (!(flags & ASL_keep_derivs)) {
		maxfwd = 0;
		want_deriv = 0;
		}
	if (!(flags & ASL_omit_all_suffixes))
		flags |= ASL_keep_all_suffixes;
	if (!(flags & ASL_forbid_missing_funcs))
		flags |= ASL_allow_missing_funcs;
	return qp_read_ASL(asl, f, flags);
	}

#ifdef __cplusplus
}
#endif
