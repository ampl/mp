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

#include "jac2dim.h"
#include "opnos.hd"

#ifdef __cplusplus
extern "C" {
#endif

#define Egulp 400

 static int com11, n_com1, ncom_togo, nvar0, nvinc, nvref;
 static int *vrefnext, *vrefx;
 static expr_if *if2list_end;
 static expr_va *varg2list_end;

#define NFHASH 23

 static derp *last_d;
 static expr *last_e;

#define efunc efunc2

 static relo *relolist, *relo2list;
 static expr_if *iflist, *if2list;
 static expr_va *varglist, *varg2list;
 static real one = 1.;
 static int nderp;
#undef nzc
 static int amax1, imap_len, k_seen, lasta, lasta0, lasta00, lastj,
	max_var, nocopy, nv0, nv01, nv0b, nv0c, nv1, nzc, nzclim;
 static int co_first = 1;
 static int *imap, *zc, *zci;

 static expr *(*holread) ANSI((EdRead*));
 extern real f_OPCPOW ANSI((expr* A_ASL));

 static ASL_fgh *asl;

 static void
ed_reset(ASL *a)
{

	a->i.memLast = a->i.memNext = 0;
	k_seen = 0;
	nvref = 0;
	vrefnext = vrefx = 0;
	relolist = relo2list = 0;
	last_d = 0;
	iflist = if2list = if2list_end = 0;
	varglist = varg2list = varg2list_end = 0;
	imap = 0;
	amax1 = imap_len = lastj = nocopy = 0;
	com11 = lasta = lasta0 = lasta00 = n_com1 = nderp = 0;
	last_e = 0;
	co_first = 1;
	}

 static void
fscream(EdRead *R, const char *name, int nargs, const char *kind)
{
	scream(R, ASL_readerr_argerr,
		"line %ld: attempt to call %s with %d %sargs\n",
		R->Line, name, nargs, kind);
	}

 static void
sorry_CLP(EdRead *R, const char *what)
{
	fprintf(Stderr,
		"Sorry, %s cannot handle %s.\n",
		progname ? progname : "", what);
	exit_ASL(R,ASL_readerr_CLP);
	}

#ifdef Double_Align
#define memadj(x) x
#else
#define memadj(x) (((x) + (sizeof(long)-1)) & ~(sizeof(long)-1))
#endif

 extern efunc f_OPPLTERM, f_OPHOL, f_OPVARVAL, f_OPFUNCALL;
 extern sfunc f_OPIFSYM;

 static void
new_derp(int a, int b, real *c)
{
	derp *d;
	if (a == nv1)
		return;
	nderp++;
	d = (derp *)mem(sizeof(derp));
	d->next = last_d;
	last_d = d;
	d->a.i = a;
	d->b.i = b;
	d->c.rp = c;
	}

 static derp *
new_relo(expr *e, derp *Dnext, int *ap)
{
	relo *r;
	derp *d;

	if (last_d != Dnext) {
		*ap = e->a;
		for(d = last_d; d->next != Dnext; d = d->next);
		d->next = 0;
		}
	else {
		last_d = 0;
		new_derp(e->a, *ap = lasta++, &one);
		}
	if (!last_d)
		return 0;
	r = (relo *)mem(sizeof(relo));
	r->next = relolist;
	r->next2 = relo2list;
	relo2list = relolist = r;
	r->D = r->Dcond = last_d;
	r->Dnext = Dnext;
	return r->D;
	}

 static relo *
new_relo1(derp *Dnext)
{
	relo *r;

	r = (relo *)mem(sizeof(relo));
	r->next = relolist;
	relolist = r;
	r->D = 0;
	r->Dnext = Dnext;
	return r;
	}

 static expr *
new_expr(int opcode, expr *L, expr *R, int deriv)
{
	expr *rv;
	extern efunc f_OP1POW, f_OP2POW, f_OPPOW;
	efunc *o;
	int L1, R1, i;

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
	rv = (expr *)mem(sizeof(expr));
	rv->op = o;
	rv->L.e = L;
	rv->R.e = R;
	if (deriv) {
		L1 = L && L->op != f_OPNUM;
		R1 = R && R->op != f_OPNUM;
		if (L1 | R1) {
			rv->a = lasta++;
			if (L1)
				new_derp(L->a, rv->a, &rv->dL);
			if (R1)
				new_derp(R->a, rv->a, &rv->dR);
			rv->bak = last_e;
			last_e = rv;
			if (R)
				rv->dL2 = rv->dLR = rv->dR2 = 0;
			else if (o == f_OP2POW)
				rv->dL2 = 2;
			else
				rv->dL2 = 0;
			if (L1)
			    if (R1)
				switch(opcode) {
					case PLUS:	i = Hv_plusLR;	break;
					case MINUS:	i = Hv_minusLR;	break;
					case MULT:	i = Hv_timesLR;	break;
					default:	i = Hv_binaryLR;
					}
			    else switch(opcode) {
					case PLUS:
					case MINUS:	i = Hv_plusL;	break;
					case MULT:	i = Hv_timesL;	break;
					case UMINUS:	i = Hv_negate;	break;
					default:	i = Hv_unary;
					}
			else
				switch(opcode) {
					case PLUS:	i = Hv_plusR;	break;
					case MINUS:	i = Hv_minusR;	break;
					case MULT:	i = Hv_timesR;	break;
					default:	i = Hv_binaryR;
					}
			rv->dO.i = i;
			}
		}
	return rv;
	}

 static char op_type[] = {
#include "op_type.hd"
	};

 static expr *
eread(EdRead *R, int deriv)
{
	arglist *al;
	argpair *ap, *da, *sap;
	char **sa;
	char *dig;
	de *d;
	derp *dsave;
	efunc *op;
	expr *L, *arg, **args, **args1, **argse, *rv;
	expr_f *rvf;
	expr_if *rvif;
	expr_n *rvn;
	expr_va *rva;
	func_info *fi;
	int *at, *nn, *nn0;
	int a0, a1, i, i1, j, j1, k, kd, kd2, ks, numargs, symargs;
	int (*Xscanf)(EdRead*, const char*, ...);
	long L1;
	plterm *p;
	real *b, **fh, *hes, r, *ra;
	unsigned int Ls;
	static real dvalue[] = {
#include "dvalue.hd"
		};

	Xscanf = xscanf;
	switch(edag_peek(R)) {
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
			rvf->dO.i = Hv_func;
			args = rvf->args;
			argse = args + j;
			k = ks = symargs = numargs = 0;
			while(args < argse) {
				arg = *args++ = eread(R, deriv);
				if ((op = arg->op) == f_OPHOL)
					symargs++;
				else if (op == (efunc*)f_OPIFSYM)
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
			if (deriv) {
				kd2 = (kd = k) ? numargs*(numargs+1) >> 1 : 0;
				rvf->a = lasta++;
				}
			else {
				kd = kd2 = 0;
				rvf->a = nv1;
				}
			Ls = sizeof(arglist)
					+ (k + kd + ks)*sizeof(argpair)
					+ kd*kd*sizeof(real *)
					+ (numargs+kd+kd2)*sizeof(real)
					+ symargs*sizeof(char *)
					+ j*sizeof(int);
			if (kd)
				Ls += numargs*sizeof(real);
			ra = (real *)mem(Ls);
			dig = 0;
			if (kd < numargs && kd)
				dig = (char*)mem(numargs);
			b = kd ? ra + numargs : ra;
			hes = b + numargs;
			al = rvf->al = (arglist *)(hes + kd2);
			al->n = numargs + symargs;
			al->nr = numargs;
			al->ra = ra;
			if (kd) {
				al->derivs = b;
				al->hes = hes;
				memset(b, 0, (numargs + kd2)*sizeof(real));
				}
			else
				al->derivs = al->hes = 0;
			al->dig = dig;
			al->funcinfo = fi->funcinfo;
			al->AE = asl->i.ae;
			al->sa = (Const char**)(sa = (char **)(al + 1));
			ap = rvf->ap = (argpair *)(sa + symargs);
			rvf->da = da = ap + k;
			sap = rvf->sap = da + kd;
			rvf->fh = fh = (real **)(sap + ks);
			at = al->at = (int *)(fh + kd*kd);
			symargs = numargs = i = 0;
			nn = nn0 = (int *)b;
			for(args = rvf->args; args < argse; at++) {
				arg = *args++;
				if ((op = arg->op) == f_OPHOL) {
					*at = --symargs;
					*sa++ = ((expr_h *)arg)->sym;
					}
				else if (op == (efunc*)f_OPIFSYM) {
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
							new_derp(arg->a,
								rvf->a, b);
							*b = 0;
							da->e = arg;
							(da++)->u.v = b;
							*nn++ = i;
							}
						}
					if (dig)
						*dig++ = j;
					b++;
					ra++;
					i++;
					}
				}
			rvf->ape = ap;
			rvf->sape = sap;
			rvf->dae = da;
			if (kd) {
				rvf->bak = last_e;
				last_e = (expr  *)rvf;
				rvf->dO.i = Hv_func;
				for(i1 = 0; i1 < kd; i1++) {
					i = nn0[i1];
					for(j1 = 0; j1 < kd; j1++) {
						j = nn0[j1];
						*fh++ = &hes[i >= j ? (i*(i+1)>>1)+j
								    : (j*(j+1)>>1)+i];
						}
					}
				}
			return (expr *)rvf;

		case 'h':
			return holread(R);

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
			rvn = (expr_n *)mem(sizeof(expr_n));
			rvn->op = f_OPNUM_ASL;
			rvn->v = r;
			return (expr *)rvn;

		case 'o':
			break;

		case 'v':
			if (Xscanf(R,"%d",&k) != 1 || k < 0)
				badline(R);
			if (k >= nvar0)
				k += nvinc;
			if (k > max_var)
				badline(R);
			if (k < nv01 && deriv && !zc[k]++)
				zci[nzc++] = k;
			return (expr *)(var_e + k);

		default:
			badline(R);
		}

	if (Xscanf(R, asl->i.opfmt, &k) != 1 || k < 0 || k >= sizeof(op_type))
		badline(R);
	switch(op_type[k]) {

		case 1:	/* unary */
			rv = new_expr(k, eread(R, deriv), 0, deriv);
			rv->dL = dvalue[k];	/* for UMINUS, FLOOR, CEIL */
			return rv;

		case 2:	/* binary */
			if (dvalue[k] == 11)
				deriv = 0;
			L = eread(R, deriv);
			rv = new_expr(k, L, eread(R, deriv), deriv);
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
			rva->L.d = d =
				(de *)mem(i*sizeof(de) + sizeof(expr *));
			rva->next = varglist;
			varglist = varg2list = rva;
			if (!last_d) {
				new_derp(lasta, lasta, &one);
				lasta++;
				}
			rva->d0 = dsave = last_d;
			rva->bak = last_e;
			a0 = a1 = lasta;
			for(j = 0; i > 0; i--, d++) {
				last_d = dsave;
				last_e = 0;
				d->e = L = eread(R, deriv);
				d->ee = last_e;
				if (L->op == f_OPNUM || L->a == nv1) {
					d->d = dsave;
					d->dv.i = nv1;
					}
				else {
					if (deriv)
						d->d = new_relo(L, dsave,
								&d->dv.i);
					if (a1 < lasta)
						a1 = lasta;
					lasta = a0;
					}
				}
			d->e = 0;	/* sentinnel expr * */
			rva->a = lasta = a1;
			last_d = dsave;
			if (deriv) {
				new_derp(0, lasta++, &one);
				/* f_MINLIST or f_MAXLIST will replace the 0 */
				rva->R.D = last_d;
				nocopy = 1;
				last_e = (expr *)rva;
				rva->dO.i = Hv_vararg;
				}
			else {
				rva->R.D = 0;
				last_e = rva->bak;
				}
			return (expr *)rva;

		case 4: /* piece-wise linear */
			i = -1;
			Xscanf(R, "%d", &i);
			if (i <= 1)
				badline(R);
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
			if (edag_peek(R) != 'v'
			 || Xscanf(R, "%d", &k) != 1
			 || k < 0 || k >= max_var)
				badline(R);
			if (k >= nvar0)
				k += nvinc;
			rv = (expr *)mem(sizeof(expr));
			rv->op = f_OPPLTERM;
			rv->L.p = p;
			rv->R.e = (expr *)(var_e + k);
			if (deriv) {
				new_derp(k, rv->a = lasta++, &rv->dL);
				rv->bak = last_e;
				last_e = rv;
				rv->dO.i = Hv_plterm;
				}
			return rv;

		case 5: /* if */
			rvif = (expr_if *)mem(sizeof(expr_if));
			rvif->op = r_ops[k];
			rvif->next = iflist;
			iflist = if2list = rvif;
			if (!last_d) {
				new_derp(lasta, lasta, &one);
				lasta++;
				}
			rvif->d0 = dsave = last_d;
			rvif->bak = last_e;
			rvif->e = eread(R, 0);
			last_e = 0;
			a0 = lasta;
			rvif->T = L = eread(R, deriv);
			rvif->Te = last_e;
			j = 0;
			if (L->op == f_OPNUM) {
				rvif->dT = dsave;
				rvif->Tv.i = nv1;
				}
			else if ((j = deriv))
				rvif->dT = new_relo(L, dsave, &rvif->Tv.i);
			a1 = lasta;
			lasta = a0;
			last_d = dsave;
			last_e = 0;
			rvif->F = L = eread(R, deriv);
			rvif->Fe = last_e;
			if (L->op == f_OPNUM) {
				rvif->dF = dsave;
				rvif->Fv.i = nv1;
				}
			else if (j)
				rvif->dF = new_relo(L, dsave, &rvif->Fv.i);
			if (lasta < a1)
				lasta = a1;
			last_d = dsave;
			if (j) {
				new_derp(0, rvif->a = lasta++, &one);
				rvif->D = last_d;
				nocopy = 1;
				last_e = (expr *)rvif;
				rvif->dO.i = Hv_if;
				}
			else {
				rvif->a = nv1;
				rvif->D = 0;
				last_e = rvif->bak;
				}
			return (expr *)rvif;

		case 11: /* OPCOUNT */
			deriv = 0;
			/* no break */
		case 6: /* sumlist */
			i = 0;
			Xscanf(R, "%d", &i);
			if (i <= 2 && (op_type[k] == 6 || i < 1))
				badline(R);
			rv = (expr *)mem(sizeof(expr) - sizeof(real)
					+ (deriv ? i+i+1 : i)*sizeof(expr *));
			rv->op = r_ops[k];
			rv->a = deriv ? lasta++ : nv1;
			rv->L.ep = args = (expr **)&rv->dR;
			if (deriv) {
				j = 0;
				do {
					*args++ = L = eread(R, deriv);
					if (L->op != f_OPNUM) {
						new_derp(L->a, rv->a, &one);
						j++;
						}
					}
					while(--i > 0);
				if (j) {
					rv->bak = last_e;
					last_e = rv;
					rv->dO.i = Hv_sumlist;
					argse = args;
					args1 = rv->L.ep;
					while(args1 < args) {
						arg = *args1++;
						if (arg->op != f_OPNUM)
							*argse++ = arg;
						}
					*argse = 0;
					}
				}
			else do
				*args++ = eread(R, deriv);
				while(--i > 0);
			rv->R.ep = args;
			return rv;
			}
	badline(R);
	return 0;
	}

 static list *
new_list(list *nxt)
{
	list *rv = (list *)mem(sizeof(list));
	rv->next = nxt;
	return rv;
	}

 static list *
crefs(void)
{
	int i;
	list *rv = 0;

	while(nzc > 0) {
		if ((i = zci[--nzc]) >= nv0) {
			rv = new_list(rv);
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
derpadjust(derp *d0, int a, derp *dnext)
{
	derp *d, *d1;
	int *r, *re;
	relo *rl;
	expr_if *il, *ile;
	expr_va *vl, *vle;
	de *de1;

	if (!(d = d0))
		return dnext;
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
derpcopy(cexp *ce, derp *dnext)
{
	derp	*d, *dprev;
	int	*map;
	derp		d00;

	if (!(d = ce->d))
		return dnext;
	map = imap;
	for(dprev = &d00; d; d = d->next) {
		new_derp(map[d->a.i], map[d->b.i], d->c.rp);
		dprev = dprev->next = last_d;
		}
	dprev->next = dnext;
	return d00.next;
	}

 static void
imap_alloc(void)
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
comsubs(int alen, cde *d)
{
	list *L;
	int a, i, j, k;
	int *r, *re;
	cexp *ce;
	derp *D, *dnext;
	relo *R;

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
		if (nzc > 0) {
			R = new_relo1(dnext);
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
				dnext = R->D = derpcopy(ce, R->D);
				}
				while(++i < nzc);
			nzc = 0;
			}
		}
	if (D || R) {
		if (!R)
			R = new_relo1(dnext);
		D = R->D = derpadjust(D, a, R->D);
		if (d->e->op != f_OPVARVAL)
			d->e->a = imap[d->e->a];
		}
	d->d = D;
	a += alen;
	d->zaplen = (a > lasta00 ? a - nv1 : 0)*sizeof(real);
	if (amax < a)
		amax = a;
	}

 static void
co_read(EdRead *R, cde *d, int wd)
{
	int alen;

	d->com11 = com11;
	d->n_com1 = n_com1;
	com11 += n_com1;
	n_com1 = 0;

	if (amax1 < lasta)
		amax1 = lasta;
	if (co_first) {
		co_first = 0;
		if (imap_len < lasta)
			imap_alloc();
		f_b = funnelfix(f_b);
		f_c = funnelfix(f_c);
		f_o = funnelfix(f_o);
		}
	if (!lastj) {
		lasta = lasta0;
		last_d = 0;
		}
	lastj = 0;
	last_e = 0;
	d->e = eread(R, wd);
	d->ee = last_e;
	alen = lasta - lasta0;
	if (imap_len < lasta)
		imap_alloc();
	comsubs(alen, d);
	}

 static linpart *
linpt_read(EdRead *R, int nlin)
{
	int (*Xscanf)(EdRead*, const char*, ...);
	linpart *L, *rv;

	if (nlin <= 0)
		return 0;
	Xscanf = xscanf;
	L = rv = (linpart *)mem(nlin*sizeof(linpart));
	do {
		if (Xscanf(R, "%d %lf", &L->v.i, &L->fac) != 2)
			badline(R);
		L++;
		}
		while(--nlin > 0);
	return rv;
	}

 static int
funnelkind(cexp *ce, int *ip)
{
	int i, j, k, nzc0, rv;
	int *vr, *vre;

	ce->vref = 0;
	if (!(nzc0 = nzc))
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
	funnel *f, **fp;
	linpart *L, *Le;
	expr *e;
	cplist *cl, *cl0;
	cexp *ce;

	ce = cexps + k - nv0;
	L = ce->L = linpt_read(R, ce->nlin = nlin);
	nocopy = 0;
	last_d = 0;
	last_e = 0;
	ce->z.i = la0 = lasta;
	nderps += nderp;
	nderp = 0;
	e = ce->e = eread(R, 1);
	if (la0 == lasta) {
		a = lasta++;
		if (e->op != f_OPNUM)
			new_derp(e->a, a, &edagread_one);
		}
	else
		a = e->a;
	var_e[k].a = a;
	ce->zlen = lasta - la0;
	for(Le = L + nlin; L < Le; L++) {
		new_derp(i = L->v.i, a, &L->fac);
		if (!zc[i]++)
			zci[nzc++] = i;
		}
	i = 0; /* only needed to shut up an erroneous warning */
	if ((fk = funnelkind(ce, &i))) {
		/* arrange to funnel */
		fp = k < nv0b ? &f_b : k < nv0c ? &f_c : &f_o;
		ce->funneled = f = (funnel *)mem(sizeof(funnel));
		f->next = *fp;
		*fp = f;
		f->ce = ce;
		if (imap_len < lasta)
			imap_alloc();
		if (fk == 1) {
			f->fulld = last_d;
			a = lasta00;
			for(i = nzc; --i >= 0; )
				if ((j = zci[i]) >= nv0)
					imap[var_e[j].a] = a++;
			if ((na = ce->zlen) || a > lasta00)
				na += a - nv1;
			f->fcde.zaplen = na*sizeof(real);
			i = nzc;
			derpadjust(last_d, a, 0);
			}
		else {
			f->fulld = 0;
			f->fcde.e = e;
			comsubs(ce->zlen, &f->fcde);
			memcpy(zci, vrefx, i*sizeof(int));
			}
		last_d = 0;
		cl0 = 0;
		while(i > 0)
			if ((a = var_e[zci[--i]].a) != nv1) {
				new_derp(a, lasta0, 0);
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
	ce->ee = last_e;
	ce->cref = crefs();
	--ncom_togo;
	}

 static void
cexp1_read(EdRead *R, int j, int k, int nlin)
{
	linpart *L, *Le;
	cexp1 *ce = cexps1 + (k - nv01);
	expr *e;
	expr_v *v;
	int la0;

	n_com1++;
	L = ce->L = linpt_read(R, ce->nlin = nlin);

	if (!lastj) {
		last_d = 0;
		if (amax1 < lasta)
			amax1 = lasta;
		lasta = lasta0;
		lastj = j;
		}
	last_e = 0;
	la0 = lasta;
	e = ce->e = eread(R, 1);
	ce->ee = last_e;
	v = var_e + k;
	if (la0 == lasta) {
		j = lasta++;
		if (e->op != f_OPNUM)
			new_derp(e->a, j, &edagread_one);
		}
	else
		j = e->a;
	v->a = j;
	v->dO.r = 0;
	for(Le = L + nlin; L < Le; L++)
		new_derp(L->v.i, j, &L->fac);
	}

 static expr *
hvadjust(expr *e)
{
	expr *e0;

	for(e0 = 0; e; e = e->bak) {
		e->fwd = e0;
		e0 = e;
		e->a = e->dO.i;
		}
	return e0;
	}

 static void
co_adjust(cde *d, int n)
{
	cde *de1;

	for(de1 = d + n; d < de1; d++)
		d->ef = hvadjust(d->ee);
	}

 static void
ifadjust(expr_if *e)
{
	for(; e; e = e->next) {
		e->Tv.rp = &adjoints[e->Tv.i];
		e->Fv.rp = &adjoints[e->Fv.i];
		e->Tf = hvadjust(e->Te);
		e->Ff = hvadjust(e->Fe);
		}
	}

 static void
vargadjust(expr_va *e)
{
	de *d;

	for(; e; e = e->next) {
		for(d = e->L.d; d->e; d++) {
			d->dv.rp = &adjoints[d->dv.i];
			d->ef = hvadjust(d->ee);
			}
		}
	}

 static void
funneladj1(funnel *f)
{
	real	*a	= adjoints;
	derp	*d;
	cplist	*cl;

	for(a = adjoints; f; f = f->next) {
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
funneladjust(VOID)
{
	cexp *c, *ce;
	linpart *L, *Le;
	c = cexps;
	for(ce = c + ncom0; c < ce; c++) {
		if ((L = c->L))
			for(Le = L + c->nlin; L < Le; L++)
				L->v.vp = (Char*)&var_e[L->v.i];
		c->ef = hvadjust(c->ee);
		}

	funneladj1(f_b);
	funneladj1(f_c);
	funneladj1(f_o);
	}

 static void
com1adjust(VOID)
{
	cexp1 *c, *ce;
	linpart *L, *Le;

	for(c = cexps1, ce = c + ncom1; c < ce; c++) {
		for(L = c->L, Le = L + c->nlin; L < Le; L++)
			L->v.vp = (Char*)&var_e[L->v.i];
		c->ef = hvadjust(c->ee);
		}
	}

 static void
adjust_compl_rhs(VOID)
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
		if (Cvar[i] && (e = C[i].e) && e->op == f_OPNUM
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
adjust(ASL *a0, int flags)
{
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
	ifadjust(iflist);
	vargadjust(varglist);
	if (ncom0)
		funneladjust();
	com1adjust();
	co_adjust(con_de, n_con);
	co_adjust(obj_de, n_obj);
	if (k_seen) {
		if (!A_vals)
			goff_comp_ASL(a0);
		else if (Fortran)
			colstart_inc_ASL(a0);
		}
	if (n_cc > nlcc && nlc < n_con
	 && !(flags & ASL_no_linear_cc_rhs_adjust))
		adjust_compl_rhs();
	}

 static void
br_read(EdRead *R, int nc, real *L, real *U, int *Cvar, int nv)
{
	ASL *asl;
	int i, inc, j, k;
	int (*Xscanf)(EdRead*, const char*, ...);

	asl = R->asl;

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

 static expr *
aholread(EdRead *R)
{
	int i, k;
	expr_h *rvh;
	char *s1;
	FILE *nl = R->nl;

	k = getc(nl);
	if (k < '1' || k > '9')
		badline(R);
	i = k - '0';
	while((k = getc(nl)) != ':') {
		if (k < '0' || k > '9')
			badline(R);
		i = 10*i + k - '0';
		}
	rvh = (expr_h *)mem(memadj(sizeof(expr_h) + i));
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
hv2init_ignore(ASL *asl, int hid_limit, int nobj, real *ow, real *y) {}

 int
fgh_read_ASL(ASL *a, FILE *nl, int flags)
{
	EdRead ER, *R;
	Jmp_buf JB;
	cgrad *cg, **cgp;
	char fname[128];
	expr_v *e;
	func_info *fi;
	int i, j, k, *ka, maxfwd1, nc, nc0, nc1, nco, ncom, nlcon, nlin;
	int no, nv, nvc, nvo, nvr, nxv, readall;
	int (*Xscanf)(EdRead*, const char*, ...);
	ograd *og, **ogp;
	real t;
	size_t *kaz, nz;
	unsigned x;

	ASL_CHECK(a, ASL_read_fgh, "fgh_read");
	flagsave_ASL(a, flags); /* includes allocation of LUv, LUrhs, A_vals or Cgrad, etc. */
	asl = (ASL_fgh*)a;
#define asl ((ASL_fgh*)a)

	ed_reset(a);
	R = EdReadInit_ASL(&ER, a, nl, 0);
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
	if (nfunc)
		func_add(a);
	if (binary_nl)
		holread = bholread;
	else
		holread = aholread;

	ncom = comb + comc + como + comc1 + como1;
	nc0 = n_con;
	nc = nc1 = nc0 + a->i.nsufext[ASL_Sufkind_con];
	no = n_obj;
	nvc = c_vars;
	nvo = o_vars;
	nco = nc + no + nlcon;
	if (no < 0 || nco <= 0)
		scream(R, ASL_readerr_corrupt,
			"ed2read: nc = %d, no = %d, nlcon = %d\n",
			nc0, no, nlcon);
	nxv = a->i.nsufext[ASL_Sufkind_var];
	nvr = n_var; /* nv for reading */
	nv0 = nv1 = nvr + nxv;
	max_var = nv = nv1 + ncom;
	combc = comb + comc;
	ncom0 = ncom_togo = combc + como;
	nzclim = ncom0 >> 3;
	ncom1 = comc1 + como1;
	nv0b = nv1 + comb;
	nv0c = nv0b + comc;
	nv01 = nv1 + ncom0;
	amax = lasta = lasta0 = lasta00 = nv1 + 1;
	ka = 0;
	kaz = 0;
	if ((maxfwd1 = maxfwd + 1) > 1)
		nvref = maxfwd1*((ncom0 < vrefGulp ? ncom0 : vrefGulp) + 1);
	x = nco*sizeof(cde) + no*sizeof(ograd *)
		+ nv*(sizeof(expr_v) + 2*sizeof(int))
		+ ncom0*sizeof(cexp)
		+ ncom1*sizeof(cexp1)
		+ nfunc*sizeof(func_info *)
		+ nvref*sizeof(int)
		+ no;
	nvar0 = a->i.n_var0;
	if (!(nvinc = a->i.n_var_ - nvar0 + nxv))
		nvar0 += ncom0 + ncom1;
	if (pi0) {
		memset(pi0, 0, nc*sizeof(real));
		if (havepi0)
			memset(havepi0, 0, nc);
		}
	if (X0)
		memset(X0, 0, nv0*sizeof(real));
	if (havex0)
		memset(havex0, 0, nv0);
	e = var_e = (expr_v *)M1zapalloc(x);
	var_ex = e + nv0;
	var_ex1 = var_ex + ncom0;
	con_de = (cde *)(e + nv);
	lcon_de = con_de + nc;
	for(k = 0; k < nv; e++) {
		e->op = f_OPVARVAL;
		e->a = k++;
		}
	obj_de = lcon_de + nlcon;
	Ograd = (ograd **)(obj_de + no);
	cexps = (cexp *)(Ograd + no);
	cexpsc = cexps + comb;
	cexpso = cexpsc + comc;
	cexps1 = (cexp1 *)(cexpse = cexps + ncom0);
	funcs = (func_info **)(cexps1 + ncom1);
	zc = (int *)(funcs + nfunc);
	zci = zc + nv;
	vrefx = zci + nv;
	objtype = (char *)(vrefx + nvref);
	if (nvref) {
		vrefnext = vrefx + maxfwd1;
		nvref -= maxfwd1;
		}
	if (n_cc && !cvar)
		cvar = (int*)M1alloc(nc*sizeof(int));
	if (cvar)
		memset(cvar, 0, nc*sizeof(int));
	last_d = 0;
	nz = 0;
	Xscanf = xscanf;
	for(;;) {
		ER.can_end = 1;
		i = edag_peek(R);
		if (i == EOF) {
			free(imap);
			adjoints = (real *)M1zapalloc(amax*Sizeof(real));
			adjoints_nv1 = &adjoints[nv1];
			adjust(a, flags);
			nzjac = nz;
			if (!Lastx)
				Lastx = (real *)M1alloc(nv0*sizeof(real));
			fclose(nl);
			nderps += nderp;
			a->p.Objval  = a->p.Objval_nomap  = obj2val_ASL;
			a->p.Objgrd  = a->p.Objgrd_nomap  = obj2grd_ASL;
			a->p.Conval  = con2val_ASL;
			a->p.Jacval  = jac2val_ASL;
			a->p.Hvcomp  = a->p.Hvcomp_nomap  = hv2comp_ASL;
			a->p.Hvcompd = hv2compd_ASL;
			a->p.Hvcomps = hv2comps_ASL;
			a->p.Hvinit  = a->p.Hvinit_nomap  = hv2init_ignore;
			a->p.Conival = a->p.Conival_nomap = con2ival_ASL;
			a->p.Congrd  = a->p.Congrd_nomap  = con2grd_ASL;
			a->p.Lconval = lcon2val_ASL;
			a->p.Xknown  = x2known_ASL;
			return prob_adj_ASL(a);
			}
		ER.can_end = 0;
		k = -1;
		switch(i) {
			case 'C':
				Xscanf(R, "%d", &k);
				if (k < 0 || k >= nc0)
					badline(R);
				co_read(R, con_de + k, 1);
				break;
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
				if (!fi->funcp && !(fi->funcp = dynlink(fname)))
					scream(R, ASL_readerr_unavail,
						"function %s not available\n",
						fname);
				funcs[i] = fi;
				break;
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
							i = ka[i]++;
							A_vals[i] = t;
							A_rownos[i] = j;
							}
						}
					else {
						while(k--) {
							if (Xscanf(R, "%d %lf",
								&i, &t) != 2)
								badline(R);
							i = kaz[i]++;
							A_vals[i] = t;
							A_rownos[i] = j;
							}
						}
					break;
					}
				cgp = Cgrad + j;
				j = 0;
				while(k--) {
					*cgp = cg = (cgrad *)mem(sizeof(cgrad));
					cgp = &cg->next;
					if (k_seen) {
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
			case 'L':
				Xscanf(R, "%d", &k);
				if (k < 0 || k >= nlcon)
					badline(R);
				co_read(R, lcon_de + k, 0);
				break;
			case 'O':
				if (Xscanf(R, "%d %d", &k, &j) != 2
				 || k < 0 || k >= no)
					badline(R);
				objtype[k] = j;
				co_read(R, obj_de + k, 1);
				break;
			case 'V':
				if (Xscanf(R, "%d %d %d", &k, &nlin, &j) != 3)
					badline(R);
				if (k >= nvar0)
					k += nvinc;
				if (k < nv0 || k >= nv)
					badline(R);
				if (j)
					cexp1_read(R, j, k, nlin);
				else
					cexp_read(R, k, nlin);
				break;
			case 'S':
				Suf_read_ASL(R, readall);
				break;
			case 'r':
				br_read(R, asl->i.n_con0, LUrhs, Urhsx, cvar, nv0);
				break;
			case 'b':
				br_read(R, asl->i.n_var0, LUv, Uvx, 0, 0);
				break;
			case 'K':
			case 'k':
				k_seen++;
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
					if (Xscanf(R, "%d %lf", &j, &t) != 2)
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
#ifdef __cplusplus
	}
#endif
