/****************************************************************
Copyright (C) 1997-2001 Lucent Technologies
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

#define Egulp 400
#define GAP_MAX 10

#ifdef Just_Linear
#define who "f_read"
#define fg_read_ASL f_read_ASL
#else
#define who "fg_read"
#ifdef __cplusplus
extern "C" {
 static real Missing_func(arglist*);
 static int compar(const void*, const void*, void*);
	}
#endif /* __cplusplus */
#endif /* Just_Linear */

#undef nzc
#undef r_ops

 typedef struct
Static {
	int _k_seen, _nv0;
	ASL *a;
	ASL_fg *asl;
	efunc **_r_ops;
#ifndef Just_Linear
	derp *_last_d;
	expr *(*_holread) ANSI((EdRead*));
	expr_if *_iflist, *_if2list, *_if2list_end;
	expr_va *_varglist, *_varg2list, *_varg2list_end;
	relo *_relolist, *_relo2list;
	int *_imap, *_vrefnext, *_vrefx, *_zc, *_zci;
	int _amax1, _co_first, _firstc1, _imap_len;
	int _last_cex, _lasta, _lasta0, _lasta00, _lastc1, _lastj;
	int _max_var, _ncom_togo, _nderp, _nocopy;
	int _nv01, _nv011, _nv0b, _nv0c, _nv1, _nvref, _nzc, _nzclim;
	int nvar0, nvinc;
#endif /* Just_Linear */
	} Static;

#define amax1		S->_amax1
#define co_first	S->_co_first
#define firstc1		S->_firstc1
#define holread		S->_holread
#define if2list		S->_if2list
#define if2list_end	S->_if2list_end
#define iflist		S->_iflist
#define imap		S->_imap
#define imap_len	S->_imap_len
#define k_seen		S->_k_seen
#define last_cex	S->_last_cex
#define last_d		S->_last_d
#define lasta		S->_lasta
#define lasta0		S->_lasta0
#define lasta00		S->_lasta00
#define lastc1		S->_lastc1
#define lastj		S->_lastj
#define max_var		S->_max_var
#define ncom_togo	S->_ncom_togo
#define nderp		S->_nderp
#define nocopy		S->_nocopy
#define nv0		S->_nv0
#define nv01		S->_nv01
#define nv011		S->_nv011
#define nv0b		S->_nv0b
#define nv0c		S->_nv0c
#define nv1		S->_nv1
#define nvref		S->_nvref
#define nzc		S->_nzc
#define nzclim		S->_nzclim
#define r_ops		S->_r_ops
#define relo2list	S->_relo2list
#define relolist	S->_relolist
#define varg2list	S->_varg2list
#define varg2list_end	S->_varg2list_end
#define varglist	S->_varglist
#define vrefnext	S->_vrefnext
#define vrefx		S->_vrefx
#define zc		S->_zc
#define zci		S->_zci

#ifdef Just_Linear

 static void
sorry_nonlin(EdRead *R)
{
	fprintf(Stderr,
		"Sorry, %s can't handle nonlinearities.\n",
		progname ? progname : "");
	exit_ASL(R,ASL_readerr_nonlin);
	}
#else /* Just_Linear */

#include "r_opn.hd"
#endif /* Just_Linear */

 static void
sorry_CLP(EdRead *R, const char *what)
{
	fprintf(Stderr,
		"Sorry, %s cannot handle %s.\n",
		progname ? progname : "", what);
	exit_ASL(R,ASL_readerr_CLP);
	}

 static Static *
ed_reset(Static *S, ASL *a)
{
	memset(S, 0, sizeof(Static));
	S->asl = (ASL_fg*)a;
	S->a = a;
	a->i.memLast = a->i.memNext = 0;
#ifndef Just_Linear
	co_first = 1;
#endif
	return S;
	}

#ifdef Double_Align
#define memadj(x) x
#else
#define memadj(x) (((x) + (sizeof(long)-1)) & ~(sizeof(long)-1))
#endif

#ifndef Just_Linear

 static void
fscream(EdRead *R, const char *name, int nargs, const char *kind)
{
	scream(R, ASL_readerr_argerr,
		"line %ld: attempt to call %s with %d %sargs\n",
		R->Line, name, nargs, kind);
	}

 static void
new_derp(Static *S, int a, int b, real *c)
{
	derp *d;
	if (a == nv1)
		return;
	nderp++;
	d = (derp *)mem_ASL(S->a, sizeof(derp));
	d->next = last_d;
	last_d = d;
	d->a.i = a;
	d->b.i = b;
	d->c.rp = c;
	}

 static derp *
new_relo(Static *S, expr *e, derp *Dnext, int *ap)
{
	relo *r;
	derp *d;

	r = (relo *)mem_ASL(S->a, sizeof(relo));
	r->next = relolist;
	r->next2 = relo2list;
	relo2list = relolist = r;
	if (last_d != Dnext) {
		*ap = e->a;
		for(d = last_d; d->next != Dnext; d = d->next);
		d->next = 0;
		}
	else {
		last_d = 0;
		new_derp(S, e->a, *ap = lasta++, &edagread_one);
		}
	r->D = r->Dcond = last_d;
	r->Dnext = Dnext;
	return r->D;
	}

 static relo *
new_relo1(Static *S, derp *Dnext)
{
	relo *r;

	r = (relo *)mem_ASL(S->a, sizeof(relo));
	r->next = relolist;
	relolist = r;
	r->D = 0;
	r->Dnext = Dnext;
	return r;
	}

 static expr *
new_expr(Static *S, int opcode, expr *L, expr *R, int deriv)
{
	expr *rv;
	efunc *o;
	int L1, R1;

	o = r_ops[opcode];
	if (o == f_OPPOW) {
		if (R->op == f_OPNUM)
			if (((expr_n *)R)->v == 2.) {
				o = f_OP2POW;
				R = 0;
				}
			else
				o = f_OP1POW;
		else if (L->op == f_OPNUM)
			o = f_OPCPOW;
		}
	rv = (expr *)mem_ASL(S->a, sizeof(expr));
	rv->op = o;
	rv->L.e = L;
	rv->R.e = R;
	rv->a = nv1;
	if (deriv) {
		L1 = L && L->op != f_OPNUM && L->a != nv1;
		R1 = R && R->op != f_OPNUM && R->a != nv1;
		if (L1 | R1) {
			rv->a = lasta++;
			if (L1)
				new_derp(S, L->a, rv->a, &rv->dL);
			if (R1)
				new_derp(S, R->a, rv->a, &rv->dR);
			}
		}
	return rv;
	}

#endif /* Just_Linear */

 static expr *
eread(EdRead *R, int deriv)
{
	ASL_fg *asl;
	Static *S;
	expr_n *rvn;
	fint L1;
	int (*Xscanf)(EdRead*, const char*, ...);
	real r;
#ifndef Just_Linear
	char **sa;
	int a0, a1, *at, i, j, k, kd, ks, numargs, symargs;
	real *b, *ra;
	expr *L, *arg, **args, **argse, *rv;
	expr_va *rva;
	plterm *p;
	de *d;
	derp *dsave;
	efunc *op;
	expr_if *rvif;
	expr_f *rvf;
	func_info *fi;
	arglist *al;
	argpair *ap, *sap;
	char *dig;
	static real dvalue[] = {
#include "dvalue.hd"
		};
#else  /* Just_Linear */
	Not_Used(deriv);
#endif /* Just_Linear */

	S = (Static *)R->S;
	asl = S->asl;
	Xscanf = xscanf;
	switch(edag_peek(R)) {
#ifdef Just_Linear
		case 'f':
		case 'h':
		case 'o':
		case 'v':
			sorry_nonlin(R);
#else
		case 'f':
			if (Xscanf(R, "%d %d", &i, &j) != 2)
				badline(R);
			fi = funcs[i];
			if (fi->nargs >= 0) {
				if (fi->nargs != j) {
 bad_nargs:
					fscream(R, fi->name, j, "");
					}
				}
			else if (-(1+fi->nargs) > j)
				goto bad_nargs;
			rvf = (expr_f *)mem(sizeof(expr_f)
					+ (j-1)*sizeof(expr *));
			rvf->op = f_OPFUNCALL;
			rvf->fi = fi;
			args = rvf->args;
			argse = args + j;
			k = ks = symargs = numargs = 0;
			while(args < argse) {
				arg = *args++ = eread(R, deriv);
				if ((op = arg->op) == f_OPHOL)
					symargs++;
				else if (op == f_OPIFSYM)
					ks++;
				else {
					numargs++;
					if (op != f_OPNUM)
						k++;
					}
				}
			symargs += ks;
			if (symargs && !(fi->ftype & 1))
				fscream(R, fi->name, symargs, "symbolic ");
			if (deriv && k) {
				kd = numargs;
				rvf->a = lasta++;
				}
			else {
				kd = 0;
				rvf->a = nv1;
				}
			ra = (real *)mem(sizeof(arglist)
					+ (k+ks)*sizeof(argpair)
					+ (numargs+kd)*sizeof(real)
					+ symargs*sizeof(char *)
					+ j*sizeof(int));
			dig = 0;
			if (k < kd)
				dig = (char*)mem(numargs);
			b = ra + kd;
			al = rvf->al = (arglist *)(b + numargs);
			al->n = numargs + symargs;
			al->nr = numargs;
			al->ra = ra;
			if (kd)
				memset(al->derivs = b, 0, kd*sizeof(real));
			else
				al->derivs = 0;
			al->hes = 0;
			al->dig = dig;
			al->funcinfo = fi->funcinfo;
			al->AE = asl->i.ae;
			al->sa = (Const char**)(sa = (char **)(al + 1));
			ap = rvf->ap = (argpair *)(sa + symargs);
			sap = rvf->sap = ap + k;
			at = al->at = (int *)(sap + ks);
			symargs = numargs = 0;
			for(args = rvf->args; args < argse; at++) {
				arg = *args++;
				if ((op = arg->op) == f_OPHOL) {
					*at = --symargs;
					*sa++ = ((expr_h *)arg)->sym;
					}
				else if (op == f_OPIFSYM) {
					*at = --symargs;
					sap->e = arg;
					(sap++)->u.s = sa++;
					}
				else {
					*at = numargs++;
					j = 1;
					if (op == f_OPNUM)
						*ra = ((expr_n *)arg)->v;
					else  {
						ap->e = arg;
						(ap++)->u.v = ra;
						if (kd) {
							j = 0;
							new_derp(S, arg->a,
								rvf->a, b);
							*b = 0;
							}
						}
					if (dig)
						*dig++ = j;
					b++;
					ra++;
					}
				}
			rvf->ape = ap;
			rvf->sape = sap;
			return (expr *)rvf;

		case 'h':
			return holread(R);
#endif /* Just_Linear */
		case 's':
			if (Xscanf(R, "%hd", &L1) != 1)
				badline(R);
			r = L1;
			goto have_r;

		case 'l':
			if (Xscanf(R, "%ld", &L1) != 1)
				badline(R);
			r = L1;
			goto have_r;

		case 'n':
			if (Xscanf(R, "%lf", &r) != 1)
				badline(R);
 have_r:
			rvn = (expr_n *)mem(size_expr_n);
			rvn->op = (efunc_n*)f_OPNUM;
			rvn->v = r;
			return (expr *)rvn;
#ifndef Just_Linear

		case 'o':
			break;

		case 'v':
			if (Xscanf(R,"%d",&k) != 1 || k < 0)
				badline(R);
			if (k >= S->nvar0)
				k += S->nvinc;
			if (k > max_var)
				badline(R);
			if (k < nv01 && deriv && !zc[k]++)
				zci[nzc++] = k;
			return (expr *)(var_e + k);

#endif /* Just_Linear */
		default:
			badline(R);
		}

#ifndef Just_Linear

	if (Xscanf(R, asl->i.opfmt, &k) != 1 || k < 0 || k >= N_OPS)
		badline(R);
	switch(optype[k]) {

		case 1:	/* unary */
			rv = new_expr(S, k, eread(R, deriv), 0, deriv);
			rv->dL = dvalue[k];	/* for UMINUS, FLOOR, CEIL */
			return rv;

		case 2:	/* binary */
			if (dvalue[k] == 11)
				deriv = 0;
			L = eread(R, deriv);
			rv = new_expr(S, k, L, eread(R, deriv), deriv);
			rv->dL = 1.;
			rv->dR = dvalue[k];	/* for PLUS, MINUS, REM */
			return rv;

		case 3:	/* vararg (min, max) */
			i = -1;
			Xscanf(R, "%d", &i);
			if (i <= 0)
				badline(R);
			rva = (expr_va *)mem(sizeof(expr_va));
			rva->op = r_ops[k];
			rva->L.d = d = (de *)mem(i*sizeof(de) + sizeof(expr *));
			rva->next = varglist;
			varglist = varg2list = rva;
			if (!last_d && deriv) {
				new_derp(S, lasta, lasta, &edagread_one);
				lasta++;
				}
			rva->d0 = dsave = last_d;
			a0 = a1 = lasta;
			for(j = 0; i > 0; i--, d++) {
				last_d = dsave;
				d->e = L = eread(R, deriv);
				if (L->op == f_OPNUM || L->a == nv1 || !deriv) {
					d->d = dsave;
					d->dv.i = nv1;
					}
				else {
					d->d = new_relo(S, L, dsave, &d->dv.i);
					j++;
					if (a1 < lasta)
						a1 = lasta;
					lasta = a0;
					}
				}
			d->e = 0;	/* sentinnel expr * */
			last_d = dsave;
			if (j) {
				rva->a = lasta = a1;
				new_derp(S, 0, lasta++, &edagread_one);
				/* f_MINLIST or f_MAXLIST will replace the 0 */
				rva->R.D = last_d;
				nocopy = 1;
				}
			else {
				rva->a = nv1;
				rva->R.D = 0;
				}
			return (expr *)rva;

		case 4: /* piece-wise linear */
			i = -1;
			Xscanf(R, "%d", &i);
			if (i <= 1)
				badline(R);
			plterms++;
			j = 2*i - 1;
			p = (plterm *)mem(sizeof(plterm) + (j-1)*sizeof(real));
			p->n = i;
			b = p->bs;
			do {
				switch(edag_peek(R)) {
					case 's':
						if (Xscanf(R,"%hd",&L1) != 1)
							badline(R);
						r = L1;
						break;
					case 'l':
						if (Xscanf(R,"%ld",&L1) != 1)
							badline(R);
						r = L1;
						break;
					case 'n':
						if (Xscanf(R,"%lf",&r) == 1)
							break;
					default:
						badline(R);
					}
				*b++ = r;
				}
				while(--j > 0);
			if (b[-2] <= 0.)
				p->z = 2*i - 2;
			else {
				b = p->bs + 1;
				while(*b <= 0.)
					b += 2;
				p->z = (b - p->bs) - 1;
				}
			rv = (expr *)mem(sizeof(expr));
			rv->op = f_OPPLTERM;
			rv->L.p = p;
			rv->R.e = L = eread(R, deriv);
			if (deriv)
				new_derp(S, L->a, rv->a = lasta++, &rv->dL);
			return rv;

		case 5: /* if */
			rvif = (expr_if *)mem(sizeof(expr_if));
			rvif->op = r_ops[k];
			rvif->next = iflist;
			iflist = if2list = rvif;
			if (!last_d && deriv) {
				new_derp(S, lasta, lasta, &edagread_one);
				lasta++;
				}
			rvif->d0 = dsave = last_d;
			rvif->e = eread(R, 0);
			a0 = lasta;
			rvif->T = L = eread(R, deriv);
			j = 0;
			if (L->op == f_OPNUM || L->a == nv1 || !(j = deriv)) {
				rvif->dT = dsave;
				rvif->Tv.i = nv1;
				}
			else
				rvif->dT = new_relo(S, L, dsave, &rvif->Tv.i);
			a1 = lasta;
			lasta = a0;
			last_d = dsave;
			rvif->F = L = eread(R, deriv);
			if (L->op == f_OPNUM || L->a == nv1 || !(j = deriv)) {
				rvif->dF = dsave;
				rvif->Fv.i = nv1;
				}
			else
				rvif->dF = new_relo(S, L, dsave, &rvif->Fv.i);
			if (lasta < a1)
				lasta = a1;
			last_d = dsave;
			if (j) {
				new_derp(S, 0, rvif->a = lasta++, &edagread_one);
				rvif->D = last_d;
				nocopy = 1;
				}
			else {
				rvif->a = nv1;
				rvif->D = 0;
				}
			return (expr *)rvif;

		case 11: /* OPCOUNT */
			deriv = 0;
			/* no break */
		case 6: /* sumlist */
			i = j = 0;
			Xscanf(R, "%d", &i);
			if (i <= 2 && (optype[k] == 6 || i < 1))
				badline(R);
			rv = (expr *)mem(sizeof(expr) - sizeof(real)
					+ i*sizeof(expr *));
			rv->op = r_ops[k];
			a0 = lasta;
			rv->L.ep = args = (expr **)&rv->dR;
			if (deriv) {
				rv->a = lasta++;
				do {
					*args++ = L = eread(R, deriv);
					if (L->op != f_OPNUM && L->a != nv1) {
						new_derp(S, L->a, rv->a,
							&edagread_one);
						j++;
						}
					}
					while(--i > 0);
				}
			else do
				*args++ = eread(R, deriv);
				while(--i > 0);
			rv->R.ep = args;
			if (!j) {
				rv->a = nv1;
				lasta = a0;
				}
			return rv;
			}
#endif /* Just_Linear */
	badline(R);
	return 0;
	}

#ifndef Just_Linear

 static list *
new_list(ASL *asl, list *nxt)
{
	list *rv = (list *)mem(sizeof(list));
	rv->next = nxt;
	return rv;
	}

 static list *
crefs(Static *S)
{
	int i;
	list *rv = 0;

	while(nzc > 0) {
		if ((i = zci[--nzc]) >= nv0) {
			rv = new_list(S->a, rv);
			rv->item.i = i;
			}
		zc[i] = 0;
		}
	return rv;
	}

 static funnel *
funnelfix(funnel *f)
{
	cexp *ce;
	funnel *fnext, *fprev;

	for(fprev = 0; f; f = fnext) {
		fnext = f->next;
		f->next = fprev;
		fprev = f;
		ce = f->ce;
		ce->z.i = ce->d->b.i;
		}
	return fprev;
	}

 static derp *
derpadjust(Static *S, derp *d0, int a, derp *dnext)
{
	derp *d, *d1;
	int *r, *re;
	relo *rl;
	expr_if *il, *ile;
	expr_va *vl, *vle;
	de *de1;
	ASL_fg *asl;

	if (!(d = d0))
		return dnext;
	asl = S->asl;
	r = imap + lasta0;
	re = imap + lasta;
	while(r < re)
		*r++ = a++;
	if (amax < a)
		amax = a;
	r = imap;
	for(;; d = d1) {
		d->a.i = r[d->a.i];
		d->b.i = r[d->b.i];
		if (!(d1 = d->next))
			break;
		}
	d->next = dnext;
	if ((rl = relo2list)) {
		relo2list = 0;
		do {
			d = rl->Dcond;
			do {
				d->a.i = r[d->a.i];
				d->b.i = r[d->b.i];
				}
				while((d = d->next));
			}
			while((rl = rl->next2));
		}
	if (if2list != if2list_end) {
		ile = if2list_end;
		if2list_end = il = if2list;
		do {
			il->Tv.i = r[il->Tv.i];
			il->Fv.i = r[il->Fv.i];
			}
			while((il = il->next) != ile);
		}
	if (varg2list != varg2list_end) {
		vle = varg2list_end;
		varg2list_end = vl = varg2list;
		do {
			for(de1 = vl->L.d; de1->e; de1++)
				de1->dv.i = r[de1->dv.i];
			}
			while((vl = vl->next) != vle);
		}
	return d0;
	}

 static derp *
derpcopy(Static *S, cexp *ce, derp *dnext)
{
	derp	*d, *dprev;
	int	*map;
	derp	d00;

	if (!(d = ce->d))
		return dnext;
	map = imap;
	for(dprev = &d00; d; d = d->next) {
		new_derp(S, map[d->a.i], map[d->b.i], d->c.rp);
		dprev = dprev->next = last_d;
		}
	dprev->next = dnext;
	return d00.next;
	}

 static void
imap_alloc(Static *S)
{
	int i, *r, *re;

	if (imap) {
		imap_len += lasta;
		imap = (int *)Realloc(imap, imap_len*Sizeof(int));
		return;
		}
	imap_len = amax1 > lasta ? amax1 : lasta;
	imap_len += 100;
	r = imap = (int *)Malloc(imap_len*Sizeof(int));
	for(i = 0, re = r + nv1+1; r < re;)
		*r++ = i++;
	}

 static int
compar(const void *a, const void *b, void *v)
{
	Not_Used(v);
	return *(int*)a - *(int*)b;
	}

 static void
comsubs(Static *S, int alen, cde *d, int **z)
{
	list *L;
	int a, i, j, k;
	int *r, *re, *z1;
	cexp *ce;
	derp *D, *dnext;
	relo *R;
	ASL_fg *asl = S->asl;

	D = last_d;
	a = lasta00;
	dnext = 0;
	R = 0;
	for(i = j = 0; i < nzc; i++)
		if ((k = zci[i]) >= nv0)
			zci[j++] = k;
		else
			zc[k] = 0;
	if ((nzc = j)) {
		for(i = 0; i < nzc; i++)
			for(L = cexps[zci[i]-nv0].cref; L; L = L->next)
				if (!zc[L->item.i]++)
					zci[nzc++] = L->item.i;
		if (nzc > 1) {
			if (nzc < nzclim)
				qsortv(zci, nzc, sizeof(int), compar, NULL);
			else for(i = nv0, j = 0; i < max_var; i++)
				if (zc[i])
					zci[j++] = i;
			}
		}
	z1 = 0; /* shut up erroneous warning */
	if (z && (k = lastc1 - firstc1 + nzc)) {
		i = (2*k + 1)*sizeof(int);
		*z = z1 = k > 20 ? (int *)M1alloc(i) : (int *)mem(i);
		*z1++ = k;
		}
	if (nzc > 0) {
		R = new_relo1(S, dnext);
		i = 0;
		do {
			j = zci[i];
			zc[j] = 0;
			ce = &cexps[j - nv0];
			if (ce->funneled)
				imap[var_e[j].a] = a++;
			else {
				r = imap + ce->z.i;
				re = r + ce->zlen;
				while(r < re)
					*r++ = a++;
				}
			if (z) {
				*z1++ = j;
				*z1++ = a - 1;
				}
			dnext = R->D = derpcopy(S, ce, R->D);
			}
			while(++i < nzc);
		nzc = 0;
		}
	if (D || R) {
		if (!R)
			R = new_relo1(S, dnext);
		D = R->D = derpadjust(S, D, a, R->D);
		if (d->e->op != f_OPVARVAL)
			d->e->a = imap[d->e->a];
		}
	if (z)
		for(i = firstc1 + nv01, j = lastc1 + nv01; i < j; i++) {
			*z1++ = i;
			*z1++ = imap[var_e[i].a];
			}
	d->d = D;
	a += alen;
	d->zaplen = (a > lasta00 ? a - nv1 : 0)*sizeof(real);
	if (amax < a)
		amax = a;
	}
#endif /* Just_Linear */

 static void
co_read(EdRead *R, cde *d, int *cexp1_end, int k, int **z, int wd)
{
#ifdef Just_Linear
	Not_Used(cexp1_end);
	Not_Used(z);
#else
	Static *S = (Static *)R->S;
	ASL_fg *asl = S->asl;
	lastc1 = last_cex - nv011;
	if (cexp1_end)
		cexp1_end[k+1] = lastc1;
	if (amax1 < lasta)
		amax1 = lasta;
	if (co_first) {
		co_first = 0;
		if (imap_len < lasta)
			imap_alloc(S);
		f_b = funnelfix(f_b);
		f_c = funnelfix(f_c);
		f_o = funnelfix(f_o);
		}
	if (!lastj) {
		lasta = lasta0;
		last_d = 0;
		}
	lastj = 0;
#endif /* Just_Linear */
	d += k;
	d->e = eread(R, wd);
#ifndef Just_Linear
	{	int alen;
		alen = lasta - lasta0;
		if (imap_len < lasta)
			imap_alloc(S);
		if (z) {
			z += k;
			*z = 0;
			}
		comsubs(S, alen, d, z);
		firstc1 = lastc1;
		}
#endif /* Just_Linear */
	}

#ifndef Just_Linear
 static linpart *
linpt_read(EdRead *R, int nlin)
{
	int (*Xscanf)(EdRead*, const char*, ...);
	linpart *L, *rv;
	ASL *asl = R->asl;

	if (nlin <= 0)
		return 0;
	L = rv = (linpart *)mem(nlin*sizeof(linpart));
	Xscanf = xscanf;
	do {
		if (Xscanf(R, "%d %lf", &L->v.i, &L->fac) != 2)
			badline(R);
		L++;
		}
		while(--nlin > 0);
	return rv;
	}

 static int
funnelkind(Static *S, cexp *ce, int *ip)
{
	int i, j, k, nzc0, rv;
	int *vr, *vre;
	ASL_fg *asl = S->asl;

	ce->vref = 0;
	if (!(nzc0 = nzc) || maxfwd <= 0)
		return 0;
	for(i = k = rv = 0; i < nzc; i++)
		if ((j = zci[i]) < nv0) {
			if (k >= maxfwd)
				goto done;
			vrefx[k++] = j;
			}
		else  {
			if (!(vr = cexps[j-nv0].vref))
				goto done;
			vre = vr + *vr;
			while(++vr <= vre)
				if (!zc[*vr]++)
					zci[nzc++] = *vr;
			}
	if (k >= nvref) {
		nvref = (maxfwd + 1)*(ncom_togo < vrefGulp
					? ncom_togo : vrefGulp);
		vrefnext = (int *)M1alloc(nvref*Sizeof(int));
		}
	vr = ce->vref = vrefnext;
	*vr = k;
	vrefnext += i = k + 1;
	nvref -= i;
	for(i = 0; i < k; i++)
		*++vr = vrefx[i];
	if (nderp > 3*k && !nocopy) {
		*ip = k;
		return 2;
		}
	else {
 done:
		if (nocopy || nderp > 3*nzc0)
			rv = 1;
		}
	while(nzc > nzc0)
		zc[zci[--nzc]] = 0;
	return rv;
	}

 static void
cexp_read(EdRead *R, int k, int nlin)
{
	int a, fk, i, j, la0, na;
	int *z1, **zp;
	funnel *f, **fp;
	linpart *L, *Le;
	expr *e;
	cplist *cl, *cl0;
	cexp *ce;
	Static *S = (Static *)R->S;
	ASL_fg *asl = S->asl;

	ce = cexps + k - nv0;
	L = ce->L = linpt_read(R, ce->nlin = nlin);
	nocopy = 0;
	last_d = 0;
	ce->z.i = la0 = lasta;
	nderps += nderp;
	nderp = 0;
	e = ce->e = eread(R, want_derivs);
	if (la0 == lasta) {
		a = lasta++;
		if (e->op != f_OPNUM)
			new_derp(S, e->a, a, &edagread_one);
		}
	else
		a = e->a;
	var_e[k].a = a;
	ce->zlen = lasta - la0;
	for(Le = L + nlin; L < Le; L++) {
		new_derp(S, i = L->v.i, a, &L->fac);
		if (!zc[i]++)
			zci[nzc++] = i;
		}
	if ((zp = zaC))
		*(zp += k - nv0) = 0;
	i = 0; /* only needed to shut up an erroneous warning */
	if ((fk = funnelkind(S, ce, &i))) {
		/* arrange to funnel */
		fp = k < nv0b ? &f_b : k < nv0c ? &f_c : &f_o;
		ce->funneled = f = (funnel *)mem(sizeof(funnel));
		f->next = *fp;
		*fp = f;
		f->ce = ce;
		if (imap_len < lasta)
			imap_alloc(S);
		if (fk == 1) {
			f->fulld = last_d;
			a = lasta00;
			z1 = 0;
			if (zp) {
				for(i = j = 0; i < nzc; i++)
					if (zci[i] >= nv0)
						j++;
				if (j) {
					i = (2*j + 1)*sizeof(int);
					*zp = z1 = j > 20
						? (int *)M1alloc(i)
						: (int *)mem(i);
					*z1++ = j;
					}
				}
			for(i = nzc; --i >= 0; )
				if ((j = zci[i]) >= nv0) {
					if (z1) {
						*z1++ = j;
						*z1++ = a;
						}
					imap[var_e[j].a] = a++;
					}
			if ((na = ce->zlen) || a > lasta00)
				na += a - nv1;
			f->fcde.zaplen = na*sizeof(real);
			i = nzc;
			derpadjust(S, last_d, a, 0);
			}
		else {
			f->fulld = 0;
			f->fcde.e = e;
			comsubs(S, ce->zlen, &f->fcde, zp);
			memcpy(zci, vrefx, i*sizeof(int));
			}
		last_d = 0;
		cl0 = 0;
		while(i > 0)
			if ((a = var_e[zci[--i]].a) != nv1) {
				new_derp(S, a, lasta0, 0);
				cl = (cplist *)mem(sizeof(cplist));
				cl->next = cl0;
				cl0 = cl;
				cl->ca.i = imap[last_d->a.i];
				cl->cfa = last_d->c.rp =
						(real *)mem(sizeof(real));
				}
		f->cl = cl0;
		var_e[k].a = lasta0++;
		lasta = lasta0;
		}
	lasta0 = lasta;
	ce->d = last_d;
	ce->cref = crefs(S);
	--ncom_togo;
	}

 static void
cexp1_read(EdRead *R, int j, int k, int nlin)
{
	Static *S = (Static *)R->S;
	ASL_fg *asl = S->asl;
	linpart *L, *Le;
	cexp1 *ce = cexps1 + (k - nv01);
	expr *e;
	int la0;

	L = ce->L = linpt_read(R, ce->nlin = nlin);

	if (!lastj) {
		last_d = 0;
		if (amax1 < lasta)
			amax1 = lasta;
		lasta = lasta0;
		lastj = j;
		}
	la0 = lasta;
	e = ce->e = eread(R, want_derivs);
	if (la0 == lasta) {
		j = lasta++;
		if (e->op != f_OPNUM)
			new_derp(S, e->a, j, &edagread_one);
		}
	else
		j = e->a;
	var_e[k].a = j;
	for(Le = L + nlin; L < Le; L++)
		new_derp(S, L->v.i, j, &L->fac);
	last_cex = k;
	}

 static void
ifadjust(real *a, expr_if *e)
{
	for(; e; e = e->next) {
		e->Tv.rp = &a[e->Tv.i];
		e->Fv.rp = &a[e->Fv.i];
		}
	}

 static void
vargadjust(real *a, expr_va *e)
{
	de *d;

	for(; e; e = e->next) {
		for(d = e->L.d; d->e; d++)
			d->dv.rp = &a[d->dv.i];
		}
	}

 static void
funneladj1(real *a, funnel *f)
{
	derp	*d;
	cplist	*cl;

	for(; f; f = f->next) {
		if ((d = f->fulld)) {
			f->fcde.d = d;
			do {
				d->a.rp = &a[d->a.i];
				d->b.rp = &a[d->b.i];
				}
				while((d = d->next));
			}
		for(cl = f->cl; cl; cl = cl->next)
			cl->ca.rp = &a[cl->ca.i];
		}
	}

 static void
funneladjust(ASL_fg *asl)
{
	cexp *c, *ce;
	linpart *L, *Le;
	real *a = adjoints;
	c = cexps;
	for(ce = c + ncom0; c < ce; c++)
		if ((L = c->L))
			for(Le = L + c->nlin; L < Le; L++)
				L->v.rp = &var_e[L->v.i].v;

	funneladj1(a, f_b);
	funneladj1(a, f_c);
	funneladj1(a, f_o);
	}

 static void
com1adjust(ASL_fg *asl)
{
	cexp1 *c, *ce;
	linpart *L, *Le;
	expr_v *v = var_e;

	for(c = cexps1, ce = c + ncom1; c < ce; c++)
		for(L = c->L, Le = L + c->nlin; L < Le; L++)
			L->v.rp = &v[L->v.i].v;
	}
#endif /* Just_Linear */

 static void
adjust_compl_rhs(ASL_fg *asl, efunc *opnum)
{
	cde *C;
	expr *e;
	int *Cvar, i, j, nc, stride;
	real *L, *U, t, t1;

	L = LUrhs;
	if ((U = Urhsx))
		stride = 1;
	else {
		U = L + 1;
		stride = 2;
		}
	C = con_de;
	Cvar = cvar;
	nc = n_con;
	for(i = nlc; i < nc; i++)
		if (Cvar[i] && (e = C[i].e) && e->op == opnum
		&& (t = ((expr_n*)e)->v) != 0.) {
			t1 = t;
			if (L[j = stride*i] > negInfinity) {
				L[j] -= t;
				t1 = 0.;
				}
			if (U[j] < Infinity) {
				U[j] -= t;
				t1 = 0.;
				}
			((expr_n*)e)->v = t1;
			}
	}

 static void
adjust(Static *S, int flags)
{
	ASL_fg *asl = S->asl;
#ifndef Just_Linear
	derp *d, **dp;
	real *a = adjoints;
	relo *r;

	for(r = relolist; r; r = r->next) {
		dp = &r->D;
		while((d = *dp)) {
			d->a.rp = &a[d->a.i];
			d->b.rp = &a[d->b.i];
			dp = &d->next;
			}
		*dp = r->Dnext;
		}
	ifadjust(adjoints, iflist);
	vargadjust(adjoints, varglist);
	if (ncom0)
		funneladjust(asl);
	com1adjust(asl);
#endif /* Just_Linear */
	if (k_seen) {
		if (!A_vals)
			goff_comp_ASL((ASL*)asl);
		else if (Fortran)
			colstart_inc_ASL((ASL*)asl);
		}
	if (n_cc > nlcc && nlc < n_con
	 && !(flags & ASL_no_linear_cc_rhs_adjust))
		adjust_compl_rhs(asl, f_OPNUM);
	}

 static void
br_read(EdRead *R, int nc, real *L, real *U, int *Cvar, int nv)
{
	int i, inc, j, k;
	int (*Xscanf)(EdRead*, const char*, ...);
	ASL *asl = R->asl;

	if (U)
		inc = 1;
	else {
		U = L + 1;
		inc = 2;
		}
	Xscanf = xscanf;
	Xscanf(R, ""); /* purge line */
	for(i = 0; i < nc; i++, L += inc, U += inc) {
		switch(edag_peek(R) - '0') {
		  case 0:
			if (Xscanf(R,"%lf %lf",L,U)!= 2)
				badline(R);
			break;
		  case 1:
			if (Xscanf(R, "%lf", U) != 1)
				badline(R);
			*L = negInfinity;
			break;
		  case 2:
			if (Xscanf(R, "%lf", L) != 1)
				badline(R);
			*U = Infinity;
			break;
		  case 3:
			*L = negInfinity;
			*U = Infinity;
			Xscanf(R, ""); /* purge line */
			break;
		  case 4:
			if (Xscanf(R, "%lf", L) != 1)
				badline(R);
			*U = *L;
			break;
		  case 5:
			if (Cvar) {
				if (Xscanf(R, "%d %d", &k, &j) == 2
				 && j > 0 && j <= nv) {
					Cvar[i] = j;
					*L = k & 2 ? negInfinity : 0.;
					*U = k & 1 ?    Infinity : 0.;
					break;
					}
				}
		  default:
			badline(R);
		  }
		}
	}

#ifndef Just_Linear

 static expr *
aholread(EdRead *R)
{
	int i, k;
	expr_h *rvh;
	char *s1;
	FILE *nl = R->nl;
	Static *S = (Static *)R->S;

	k = getc(nl);
	if (k < '1' || k > '9')
		badline(R);
	i = k - '0';
	while((k = getc(nl)) != ':') {
		if (k < '0' || k > '9')
			badline(R);
		i = 10*i + k - '0';
		}
	rvh = (expr_h *)mem_ASL(R->asl, memadj(sizeof(expr_h) + i));
	for(s1 = rvh->sym;;) {
		if ((k = getc(nl)) < 0) {
			fprintf(Stderr,
				 "Premature end of file in aholread, line %ld of %s\n",
					R->Line, R->asl->i.filename_);
				exit_ASL(R,1);
			}
		if (k == '\n') {
			R->Line++;
			if (!i)
				break;
			}
		if (--i < 0)
			badline(R);
		*s1++ = k;
		}
	*s1 = 0;
	rvh->op = f_OPHOL;
	rvh->a = nv1;
	return (expr *)rvh;
	}

 static expr *
bholread(EdRead *R)
{
	int i;
	expr_h *rvh;
	char *s;
	Static *S = (Static *)R->S;
	ASL_fg *asl = S->asl;

	if (xscanf(R, "%d", &i) != 1)
		badline(R);
	rvh = (expr_h *)mem(memadj(sizeof(expr_h) + i));
	s = rvh->sym;
	if (fread(s, i, 1, R->nl) != 1)
		badline(R);
	s[i] = 0;
	rvh->op = f_OPHOL;
	rvh->a = nv1;
	for(;;) switch(*s++) {
			case 0: goto break2; /* so we return at end of fcn */
			case '\n': R->Line++;
			}
 break2:
	return (expr *)rvh;
	}

 static void
nlvzap(Static *S, int i, int j)
{
	int n = nv1;
	expr_v *v = S->var_e;

	i -= j;
	while(--j >= 0)
		v[i+j].a = n;
	}

 static real
Missing_func(arglist *al)
{
	AmplExports *ae = al->AE;
	func_info *fi = (func_info*)al->funcinfo;

	char *s = (char*)(*ae->Tempmem)(al->TMI, strlen(fi->name) + 64);
	(*ae->SprintF)(al->Errmsg = s,
		"Attempt to call unavailable function %s.",
		fi->name);
	return 0.;
	}
#endif /* Just_Linear */

 int
fg_read_ASL(ASL *a, FILE *nl, int flags)
{
	cgrad *cg, **cgp;
	expr_v *e;
	int i, i1, j, k, *ka, kseen, nc, nc0, nco, nlcon, no, nv, nvc, nvo, nvr, nxv, readall;
	int (*Xscanf)(EdRead*, const char*, ...);
	ograd *og, **ogp;
	real t;
	size_t *kaz, nz;
	unsigned x;
#ifdef Just_Linear
#undef nv1
	int nv1;
#define ASL_readtype ASL_read_f
#else /* Just_Linear */
#define ASL_readtype ASL_read_fg
	int maxfwd1, ncom;
	func_info *fi;
	char fname[128];
	int nlin;
#endif /* Just_Linear */
	EdRead ER, *R;
	Static SS, *S;
	Jmp_buf JB;

	ASL_CHECK(a, ASL_readtype, who);
	flagsave_ASL(a, flags); /* includes allocation of LUv, LUrhs, A_vals or Cgrad, etc. */
	S = ed_reset(&SS, a);
	R = EdReadInit_ASL(&ER, a, nl, S);
	if (flags & ASL_return_read_err) {
		a->i.err_jmp_ = &JB;
		i = setjmp(JB.jb);
		if (i) {
			a->i.err_jmp_ = 0;
			return i;
			}
		}
	nlcon = a->i.n_lcon_;
	if (nlcon && !(flags & ASL_allow_CLP)) {
		if (a->i.err_jmp_)
			return ASL_readerr_CLP;
		sorry_CLP(R, "logical constraints");
		}
	readall = flags & ASL_keep_all_suffixes;
	Xscanf = a->i.xscanf_;
#define asl ((ASL_fg*)a)
	ER.lineinc = 1;
	if (!size_expr_n)	size_expr_n = sizeof(expr_n);
#ifndef Just_Linear
	if (!(SS._r_ops = asl->I.r_ops_))
		SS._r_ops = r_ops_ASL;
	if (c_cexp1st)
		*c_cexp1st = 0;
	if (o_cexp1st)
		*o_cexp1st = comc1;
	if (nfunc)
		func_add(a);
	if (binary_nl)
		holread = bholread;
	else
		holread = aholread;

	ncom = comb + comc + como + comc1 + como1;
#endif /* Just_Linear */
	nc0 = n_con;
	nc = nc0 + a->i.nsufext[ASL_Sufkind_con];
	no = n_obj;
	nvc = c_vars;
	nvo = o_vars;
	nco = nc + no + nlcon;
	if (no < 0 || nco <= 0)
		scream(R, ASL_readerr_corrupt,
			"edagread: nc = %d, no = %d, nlcon = %d\n",
			nc0, no, nlcon);
	if (pi0) {
		memset(pi0, 0, nc*sizeof(real));
		if (havepi0)
			memset(havepi0, 0, nc);
		}
	nxv = a->i.nsufext[ASL_Sufkind_var];
	nvr = n_var; /* nv for reading */
	nv0 = nv1 = nvr + nxv;
#ifdef Just_Linear
	nv = nv1;
	x = nco*sizeof(cde) + no*sizeof(ograd *) + nv*sizeof(expr_v) + no;
#else
	max_var = nv = nv1 + ncom;
	combc = comb + comc;
	ncom0 = ncom_togo = combc + como;
	nzclim = ncom0 >> 3;
	ncom1 = comc1 + como1;
	nv0b = nv1 + comb;
	nv0c = nv0b + comc;
	nv01 = nv1 + ncom0;
	last_cex = nv011 = nv01 - 1;
	amax = lasta = lasta0 = lasta00 = nv1 + 1;
	if ((maxfwd1 = maxfwd + 1) > 1)
		nvref = maxfwd1*((ncom0 < vrefGulp ? ncom0 : vrefGulp) + 1);
	x = nco*sizeof(cde) + no*sizeof(ograd *)
		+ nv*(sizeof(expr_v) + 2*sizeof(int))
		+ ncom0*sizeof(cexp)
		+ ncom1*sizeof(cexp1)
		+ nfunc*sizeof(func_info *)
		+ nvref*sizeof(int)
		+ no;
	SS.nvar0 = a->i.n_var0;
	if (!(SS.nvinc = a->i.n_var_ - SS.nvar0 + nxv))
		SS.nvar0 += ncom0 + ncom1;
#endif /* Just_Linear */
	if (X0)
		memset(X0, 0, nv1*sizeof(real));
	if (havex0)
		memset(havex0, 0, nv1);
	e = var_e = (expr_v *)M1zapalloc(x);
	con_de = (cde *)(e + nv);
	lcon_de = con_de + nc;
	obj_de = lcon_de + nlcon;
	Ograd = (ograd **)(obj_de + no);
#ifdef Just_Linear
	objtype = (char *)(Ograd + no);
#else
	var_ex = e + nv1;
	var_ex1 = var_ex + ncom0;
	for(k = 0; k < nv; e++) {
		e->op = f_OPVARVAL;
		e->a = k++;
		}
	if (skip_int_derivs) {
		if (nlvbi)
			nlvzap(S, nlvb, nlvbi);
		if (nlvci)
			nlvzap(S, nlvb+nlvc, nlvci);
		if (nlvoi)
			nlvzap(S, nlvb+nlvc+nlvo, nlvoi);
		}
	cexps = (cexp *)(Ograd + no);
	cexps1 = (cexp1 *)(cexps + ncom0);
	funcs = (func_info **)(cexps1 + ncom1);
	zc = (int*)(funcs + nfunc);
	zci = zc + nv;
	vrefx = zci + nv;
	objtype = (char *)(vrefx + nvref);
	if (nvref) {
		vrefnext = vrefx + maxfwd1;
		nvref -= maxfwd1;
		}
	last_d = 0;
#endif
	if (n_cc && !cvar)
		cvar = (int*)M1alloc(nc*sizeof(int));
	if (cvar)
		memset(cvar, 0, nc*sizeof(int));
	ka = 0;
	kaz = 0;
	nz = 0;
	j = kseen = 0;
	for(;;) {
		ER.can_end = 1;
		i = edag_peek(R);
		if (i == EOF) {
#ifndef Just_Linear
			free(imap);
			/* Make amax long enough for nlc to handle */
			/* var_e[i].a for common variables i. */
			if (ncom0) {
				i = comb + como;
				if (i < combc)
					i = combc;
				if ((i += nv1 + 1) > amax)
					amax = i;
				}
			adjoints = (real *)M1zapalloc(amax*Sizeof(real));
			adjoints_nv1 = &adjoints[nv1];
			nderps += nderp;
#endif /* Just_Linear */
			adjust(S, flags);
			nzjac = nz;
			if (!Lastx)
				Lastx = (real *)M1alloc(nv1*sizeof(real));
			fclose(nl);
#ifndef Just_Linear
			a->p.Objval  = a->p.Objval_nomap  = obj1val_ASL;
			a->p.Objgrd  = a->p.Objgrd_nomap  = obj1grd_ASL;
			a->p.Conval  = con1val_ASL;
			a->p.Jacval  = jac1val_ASL;
			a->p.Conival = a->p.Conival_nomap = con1ival_ASL;
			a->p.Congrd  = a->p.Congrd_nomap  = con1grd_ASL;
			a->p.Lconval = lcon1val_ASL;
			a->p.Xknown  = x1known_ASL;
#endif /* Just_Linear */
			return prob_adj_ASL(a);
			}
		ER.can_end = 0;
		k = -1;
		switch(i) {
			case 'C':
				Xscanf(R, "%d", &k);
				if (k < 0 || k >= nc0)
					badline(R);
				co_read(R,con_de,c_cexp1st,k,zac,want_derivs);
				break;
#ifdef Just_Linear
			case 'F':
			case 'V':
				sorry_nonlin(R);
#else
			case 'F':
				if (Xscanf(R, "%d %d %d %127s",
						&i, &j, &k, fname) != 4
				|| i < 0 || i >= nfunc)
					badline(R);
				if ((fi = func_lookup(a, fname,0))) {
					if (fi->nargs != k && fi->nargs >= 0
					 && (k >= 0 || fi->nargs < -(k+1)))
						scream(R, ASL_readerr_argerr,
				"function %s: disagreement of nargs: %d and %d\n",
					 		fname,fi->nargs, k);
					}
				else {
					fi = (func_info *)mem(sizeof(func_info));
					fi->ftype = j;
					fi->nargs = k;
					fi->funcp = 0;
					fi->name = (Const char *)strcpy((char*)
						mem(memadj(strlen(fname)+1)),
						fname);
					}
				if (!fi->funcp && !(fi->funcp = dynlink(fname))){
					if (!(flags & ASL_allow_missing_funcs))
					    scream(R, ASL_readerr_unavail,
						"function %s not available\n",
						fname);
					fi->funcp = Missing_func;
					fi->funcinfo = (Char*)fi;
					}
				funcs[i] = fi;
				break;
			case 'L':
				Xscanf(R, "%d", &k);
				if (k < 0 || k >= nlcon)
					badline(R);
				co_read(R, lcon_de, 0, k, 0, 0);
				break;
			case 'V':
				if (Xscanf(R, "%d %d %d", &k, &nlin, &j) != 3)
					badline(R);
				if (k >= SS.nvar0)
					k += SS.nvinc;
				if (k < nvr || k >= nv)
					badline(R);
				if (j)
					cexp1_read(R, j, k, nlin);
				else
					cexp_read(R, k, nlin);
				break;
#endif /* Just_Linear */
			case 'G':
				if (Xscanf(R, "%d %d", &j, &k) != 2
				|| j < 0 || j >= no || k <= 0 || k > nvo)
					badline(R);
				ogp = Ograd + j;
				while(k--) {
					*ogp = og = (ograd *)mem(sizeof(ograd));
					ogp = &og->next;
					if (Xscanf(R, "%d %lf", &i, &og->coef) != 2)
						badline(R);
					og->varno = i;
					}
				*ogp = 0;
				break;
			case 'J':
				if (Xscanf(R, "%d %d", &j, &k) != 2
				|| j < 0 || j >= nc0 || k <= 0 || k > nvc)
					badline(R);
				nz += k;
				if (A_vals) {
					j += Fortran;
					if (ka) {
						while(k--) {
							if (Xscanf(R, "%d %lf",
								&i, &t) != 2)
								badline(R);
							i1 = ka[i]++;
							A_vals[i1] = t;
							A_rownos[i1] = j;
							}
						}
					else {
						while(k--) {
							if (Xscanf(R, "%d %lf",
								&i, &t) != 2)
								badline(R);
							i1 = kaz[i]++;
							A_vals[i1] = t;
							A_rownos[i1] = j;
							}
						}
					break;
					}
				cgp = Cgrad + j;
				while(k--) {
					*cgp = cg = (cgrad *)mem(sizeof(cgrad));
					cgp = &cg->next;
					if (kseen) {
						if (Xscanf(R, "%d %lf", &i, &cg->coef) != 2)
							badline(R);
						}
					else
						if (Xscanf(R, "%d %d %lf", &i, &j,
							    &cg->coef) != 3)
							badline(R);
					cg->varno = i;
					cg->goff = j;
					}
				*cgp = 0;
				break;
			case 'O':
				if (Xscanf(R, "%d %d", &k, &j) != 2
				 || k < 0 || k >= no)
					badline(R);
				objtype[k] = j;
				co_read(R,obj_de,o_cexp1st,k,zao,want_derivs);
				break;
			case 'S':
				Suf_read_ASL(R, readall);
				break;
			case 'r':
				br_read(R, asl->i.n_con0, LUrhs, Urhsx, cvar, nvr);
				break;
			case 'b':
				br_read(R, asl->i.n_var0, LUv, Uvx, 0, 0);
				break;
			case 'K':
			case 'k':
				k_seen = ++kseen;
				if (ka_read_ASL(a, R, i, &ka, &kaz))
					badline(R);
				break;
			case 'x':
				if (!Xscanf(R,"%d",&k)
				|| k < 0 || k > nvr)
					badline(R);
				if (!X0 && want_xpi0 & 1) {
					x = nv1*sizeof(real);
					if (want_xpi0 & 4)
						x += nv1;
					X0 = (real *)M1zapalloc(x);
					if (want_xpi0 & 4)
						havex0 = (char*)(X0 + nv1);
					}
				while(k--) {
					if (Xscanf(R, "%d %lf", &j, &t) != 2
					 || j < 0 || j >= nvr)
						badline(R);
					if (X0) {
						X0[j] = t;
						if (havex0)
							havex0[j] = 1;
						}
					}
				break;
			case 'd':
				if (!Xscanf(R,"%d",&k)
				|| k < 0 || k > nc0)
					badline(R);
				if (!pi0 && want_xpi0 & 2) {
					x = nc*sizeof(real);
					if (want_xpi0 & 4)
						x += nc;
					pi0 = (real *)M1zapalloc(x);
					if (want_xpi0 & 4)
						havepi0 = (char*)(pi0 + nc);
					}
				while(k--) {
					if (Xscanf(R, "%d %lf", &j, &t) != 2
					 || j < 0 || j >= nc0)
						badline(R);
					if (pi0) {
						pi0[j] = t;
						if (havepi0)
							havepi0[j] = 1;
						}
					}
				break;
			default:
				badline(R);
			}
		}
	}
