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

#ifdef DEBUG
#include "assert.h"
#else
#define assert(x) /*nothing*/
#endif

#ifdef PSHVREAD
#include "jacpdim.h"
#define efunc efunc2
#include "opnos.hd"
#define PSHV(x) x
#define asltype ASL_read_pfgh
#define who "pfgh_read"
#define ASLTYPE ASL_pfgh
#else
#define PSHV(x) /*nothing*/
#include "asl_pfg.h"
#define asltype ASL_read_pfg
#define who "pfg_read"
#define ASLTYPE ASL_pfg
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define GAP_MAX 10

#ifdef PSHVREAD
 static void hes_setup ANSI((ASL*, int, int, int, int, int));
/* #define Static to avoid confusion in pi. */
#define Static Static_pshv
#else
#define Static Static_psed
#endif
typedef struct Static Static;
 static ograd *cotermwalk ANSI((Static*, expr **ep, ps_func *f, int wantg, int omitdv));
 static efunc *OPNUM, *OPVARVAL;

#define NFHASH 23

#undef f_OPNUM
#include "r_opn0.hd"

#undef nzc

#if defined(IEEE_MC68k) || defined(IBM)
  enum { W0 = 0, W1 = 1 };
#else
  enum { W0 = 1, W1 = 0 };
#endif

/*******/

 typedef struct
expr_vx {
	expr_v v;
	linarg *la;
	int a0;
	int a1;
	} expr_vx;

 typedef struct
Elemtemp {
	unsigned int esize;
	int nmax;
	int k;
	void **mp;
	} Elemtemp;

 typedef struct
PSfind {
	ps_func *f;
	Elemtemp *b, *g;
	} PSfind;

 typedef struct EU EU;
 struct EU { efunc_n *op; union {real v; EU *p; } u; };

 struct
Static {
	ASLTYPE *asl;
	ASL *a;
	Elemtemp *_last_psb_elem;
	derp *_last_d;
	expr *(*_holread) ANSI((EdRead*));
	expr *_expr_free, **_slscratch;	/* used in awalk(sumlist) */
	expr *_last_e;
	expr_if *_iflist, *_if2list, *_if2list_end;
	EU *_expr_n_free;
	expr_v **_larvlist;
	expr_v **_varp;
	expr_va *_varglist, *_varg2list, *_varg2list_end;
	int *_imap, *_vrefnext, *_vrefx;
	int *_zc, *zc1, *_zci, *zci1, *_zl;
	size_t _lthashmask, _nrange, _rangehashmask;
	int klthash, krangehash;
	int _allJ;
	int _amax1;
	int _cexp_k;
	int _cexp_n;
	int _conno;
	int _groupno;
	int _imap_len;
	int _k_Elemtemp;
	int _k_seen;
	int _kimap;
	int _kzc;
	int _lasta;
	int _lasta0;
	int _lasta00;
	int _max_var;
	int _ncom;
	int _ncom_togo;
	int _nderp;
	int _noa;
	int _nocopy;
	int _nndv;
	int _nsce;
	int _nv0b;
	int _nv0c;
	int _nv0r;
	int _nv0x;
	int _nv1;
	int _nvref;
	int _nzc;
	int _nzclim;
	int _slmax;
	int _slscratchlev;
	int _termno;		/* current term number */
	int _wantCgroups;
	int _wantOgroups;
	int _zc_lim;
	int nvar0, nvinc;
	int size_exprn;
	la_ref *_laref_free;
	linarg **_lthash;		/* hash table */
	linarg *_ltfree;		/* free list */
	linarg *_tlist;		/* linargs in this term */
	ograd *_freeog;		/* free list */
	range **_rangehash;
	real *_rnz;
	real _lt_scale;
	relo *_relolist, *_relo2list;
	EdRead *R;
	};

#define allJ		S->_allJ
#define amax1		S->_amax1
#define cexp_k		S->_cexp_k
#define cexp_n		S->_cexp_n
#define Conno		S->_conno
#define expr_free	S->_expr_free
#define expr_n_free	S->_expr_n_free
#define freeog		S->_freeog
#define Groupno		S->_groupno
#define holread		S->_holread
#define if2list		S->_if2list
#define if2list_end	S->_if2list_end
#define iflist		S->_iflist
#define imap		S->_imap
#define imap_len	S->_imap_len
#define k_Elemtemp	S->_k_Elemtemp
#define k_seen		S->_k_seen
#define kimap		S->_kimap
#define kzc		S->_kzc
#define laref_free	S->_laref_free
#define larvlist	S->_larvlist
#define last_d		S->_last_d
#define last_e		S->_last_e
#define last_psb_elem	S->_last_psb_elem
#define lasta		S->_lasta
#define lasta0		S->_lasta0
#define lasta00		S->_lasta00
#define lt_scale	S->_lt_scale
#define ltfree		S->_ltfree
#define lthash		S->_lthash
#define lthashmask	S->_lthashmask
#define max_var		S->_max_var
#define max_var1	asl->P.max_var1_
#define Ncom		S->_ncom
#define ncom_togo	S->_ncom_togo
#define nderp		S->_nderp
#define noa		S->_noa
#define nocopy		S->_nocopy
#define nndv		S->_nndv
#define nrange		S->_nrange
#define nsce		S->_nsce
#define nv0b		S->_nv0b
#define nv0c		S->_nv0c
#define nv0r		S->_nv0r
#define nv0x		S->_nv0x
#define nv1		S->_nv1
#define nvref		S->_nvref
#define nzc		S->_nzc
#define nzclim		S->_nzclim
#define rangehash	S->_rangehash
#define rangehashmask	S->_rangehashmask
#define relo2list	S->_relo2list
#define relolist	S->_relolist
#define rnz		S->_rnz
#define slmax		S->_slmax
#define slscratch	S->_slscratch
#define slscratchlev	S->_slscratchlev
#define Termno		S->_termno
#define tlist		S->_tlist
#define varg2list	S->_varg2list
#define varg2list_end	S->_varg2list_end
#define varglist	S->_varglist
#define varp		S->_varp
#define vrefnext	S->_vrefnext
#define vrefx		S->_vrefx
#define wantCgroups	S->_wantCgroups
#define wantOgroups	S->_wantOgroups
#define zc		S->_zc
#define zc_lim		S->_zc_lim
#define zci		S->_zci
#define zl		S->_zl

static list *crefs ANSI((Static*));

 static Static *
S_init(Static *S, ASLTYPE *asl)
{
	memset(S, 0, sizeof(Static));
	S->asl = asl;
	S->a = (ASL*)asl;
	return S;
	}

 static void *
new_mblkzap(ASLTYPE *asl, int k)
{
	void *rv = new_mblk_ASL((ASL*)asl,k);
	memset(rv, 0, sizeof(void*)<<k);
	return rv;
	}

 static void
sorry_CLP(EdRead *R, const char *what)
{
	fprintf(Stderr,
		"Sorry, %s cannot handle %s.\n",
		progname ? progname : "", what);
	exit_ASL(R,ASL_readerr_CLP);
	}

 static void
ed_reset(ASLTYPE *asl)
{
	asl->i.memLast = asl->i.memNext = 0;
	memset(asl->mblk_free, 0, MBLK_KMAX*sizeof(char*));
	memset(&asl->P.merge + 1, 0, sizeof(ps_info) - sizeof(long));
#ifdef PSHVREAD
	asl->P.hop_free = 0;
	asl->P.hes_setup_called = 0;
#endif
	}

#ifdef Double_Align
#define memadj(x) x
#else
#define memadj(x) (((x) + (sizeof(long)-1)) & ~(sizeof(long)-1))
#endif
#ifdef DEBUG
 static expr *opzork;
 static derp *dzork;
 static ograd *ogzork;
 static expr *ezork;
 static EU *enzork;
 static int dzork1, dzork2, izork = -1;
#endif

 static Elemtemp *
new_Elemtemp(Static *S, unsigned int esize, void **mp)
{
	Elemtemp *e;
	int k;
	ASL *asl = S->a;

	e = (Elemtemp *)new_mblk(k_Elemtemp);
	e->esize = esize;
	e->mp = mp;
	e->k = k = htcl(8*esize);
	*mp = new_mblk(k);
	e->nmax = (sizeof(void*) << k) / esize;
	return e;
	}

 static void
del_Elemtemp(Static *S, Elemtemp *e)
{
	ASL *asl = S->a;
	del_mblk(e->k, *e->mp);
	del_mblk(k_Elemtemp, e);
	}

 static void
upgrade_Elemtemp(Static *S, Elemtemp *e)
{
	void *m, *m0;
	int k;
	ASL *asl = S->a;

	k = e->k++;
	memcpy(m = new_mblk(e->k), m0 = *e->mp, e->esize * e->nmax);
	del_mblk(k++, m0);
	*e->mp = m;
	e->nmax = (sizeof(void*) << k) / e->esize;
	}

 static void
free_laref(Static *S, la_ref **L)
{
	la_ref *L1, *L2;

	if ((L1 = *L)) {
		while((L2 = L1->next))
			L1 = L2;
		L1->next = laref_free;
		laref_free = *L;
		*L = 0;
		}
	}

 static la_ref *
new_laref(Static *S, la_ref *nxt)
{
	la_ref *rv;

	if ((rv = laref_free))
		laref_free = rv->next;
	else
		rv = (la_ref *)mem_ASL(S->a, sizeof(la_ref));
	rv->next = nxt;
	return rv;
	}

 static list *
new_list(ASL *asl, list *nxt)
{
	list *rv = (list *)mem(sizeof(list));
	rv->next = nxt;
	return rv;
	}

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
	nderp++;
	d = (derp *)mem_ASL(S->a, sizeof(derp));
	d->next = last_d;
	last_d = d;
	d->a.i = a;
	d->b.i = b;
	d->c.rp = c;
#ifdef DEBUG
	if (d == dzork)
		printf("");
	if (++dzork1 == dzork2)
		printf("");
#endif
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
	if (!last_d)
		return 0;
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

 static void
efree(Static *S, expr *e)
{
	EU *en;
	expr **args, **argse, *e1;
 top:
	switch(optypeb[Intcast e->op]) {

	 case 2: /* binary */
		efree(S, e->R.e);
		/* no break */

	 case 1: /* unary */
		e1 = e->L.e;
		e->L.e = expr_free;
		expr_free = e;
		e = e1;
		goto top;

	 case 6: /* sumlist */
		args = e->L.ep;
		argse = e->R.ep;
		while(args < argse)
			efree(S, *args++);
		e->L.e = expr_free;
		expr_free = e;
		break;

	 case 9:
		en = (EU*)e;
		en->u.p = expr_n_free;
		expr_n_free = en;
	 }
	}

 static expr *
new_expr(Static *S, int o, expr *L, expr *R)
{
	expr *rv;

	if ((rv = expr_free))
		expr_free = rv->L.e;
	else
		rv = (expr *)mem_ASL(S->a, sizeof(expr));
	PSHV(rv->dL2 = 0);
	if (o == f_OPPOW) {
		if (Intcast R->op == f_OPNUM) {
			if (((expr_n *)R)->v == 2.) {
				o = f_OP2POW;
				R = 0;
				PSHV(rv->dL2 = 2);
				}
			else
				o = f_OP1POW;
			}
		else if (Intcast L->op == f_OPNUM)
			o = f_OPCPOW;
		}
	rv->op = (efunc *)(size_t)o;
	rv->L.e = L;
	rv->R.e = R;
#ifdef DEBUG
	if (rv == ezork)
		printf("");
#endif
	return rv;
	}

 static void
dexpr(Static *S, expr *e, expr *L, expr *R)
{
	int L1, R1;
#ifdef PSHVREAD
	int i, opcode;
#endif

	e->a = noa;
	L1 = L && L->op != OPNUM && L->a != noa;
	R1 = R && R->op != OPNUM && R->a != noa;
	if (L1 | R1) {
		if (L1)
			new_derp(S, L->a, lasta, &e->dL);
		if (R1)
			new_derp(S, R->a, lasta, &e->dR);
		e->a = lasta++;
#ifdef PSHVREAD
		e->bak = last_e;
		last_e = e;
		opcode = Intcast e->op;
		if (R)
			e->dLR = e->dR2 = 0;
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
				case DIV:
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
		e->dO.i = i;
#endif
		}
	}

 static real dvalue[] = {
#include "dvalue.hd"
		};

 static expr_n *
new_expr_n(Static *S, real t)
{
	EU *en;

	if ((en = expr_n_free))
		expr_n_free = en->u.p;
	else
		en = (EU *)mem_ASL(S->a, S->size_exprn);
#ifdef DEBUG
	if (en == enzork)
		printf("");
#endif
	en->u.v = t;
	en->op = (efunc_n *)f_OPNUM;
	return (expr_n*)en;
	}

 static expr *
eread(EdRead *R)
{
	ASLTYPE *asl;
	Static *S;
	arglist *al;
	argpair *ap, *sap;
	char **sa;
	de *d;
	expr *L, *arg, **args, **argse, *rv;
	expr_f *rvf;
	expr_if *rvif;
	expr_va *rva;
	fint L1;
	func_info *fi;
	int *at, i, j, k, ks, numargs, op, symargs;
	int (*Xscanf)(EdRead*, const char*, ...);
	plterm *p;
	real *b, *ra;
	real r;

	S = (Static *)R->S;
	asl = S->asl;
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
			rvf->op = (efunc *)f_OPFUNCALL;
			rvf->fi = fi;
			PSHV(rvf->dO.i = Hv_func);
			args = rvf->args;
			argse = args + j;
			k = ks = symargs = numargs = 0;
			while(args < argse) {
				arg = *args++ = eread(R);
				if ((op = Intcast arg->op) == f_OPHOL)
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
			ra = (real *)mem(sizeof(arglist)
					+ (k+ks)*sizeof(argpair)
					+ numargs*sizeof(real)
					+ symargs*sizeof(char *)
					+ j*sizeof(int));
			al = rvf->al = (arglist *)(ra + numargs);
			al->n = numargs + symargs;
			al->nr = numargs;
			al->ra = ra;
			al->derivs = al->hes = 0;
			al->funcinfo = fi->funcinfo;
			al->AE = asl->i.ae;
			al->sa = (Const char**)(sa = (char **)(al + 1));
			ap = rvf->ap = (argpair *)(sa + symargs);
			sap = rvf->sap = ap + k;
			at = al->at = (int *)(sap + ks);
			symargs = numargs = 0;
			for(args = rvf->args; args < argse; at++) {
				arg = *args++;
				if ((op = Intcast arg->op) == f_OPHOL) {
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
					if (op == f_OPNUM)
						*ra = ((expr_n *)arg)->v;
					else  {
						ap->e = arg;
						(ap++)->u.v = ra;
						}
					ra++;
					}
				}
			rvf->ape = ap;
			rvf->sape = sap;
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
			return (expr *)new_expr_n(S, r);

		case 'o':
			break;

		case 'v':
			if (Xscanf(R,"%d",&k) != 1 || k < 0)
				badline(R);
			if (k >= S->nvar0)
				k += S->nvinc;
			if (k > max_var)
				badline(R);
			return (expr *)(var_e + k);

		default:
			badline(R);
		}

	if (Xscanf(R, asl->i.opfmt, &k) != 1 || k < 0 || k >= N_OPS)
		badline(R);
	switch(optype[k]) {

		case 1:	/* unary */
			rv = new_expr(S, k, eread(R), 0);
			return rv;

		case 2:	/* binary */
			L = eread(R);
			rv = new_expr(S, k, L, eread(R));
			return rv;

		case 3:	/* vararg (min, max) */
			i = -1;
			Xscanf(R, "%d", &i);
			if (i <= 0)
				badline(R);
			rva = (expr_va *)mem(sizeof(expr_va));
			PSHV(rva->val = 0;)
			rva->op = (efunc *)(size_t)k;
			rva->L.d = d = (de *)mem(i*sizeof(de) + sizeof(expr *));
			rva->next = varglist;
			varglist = rva;
			for(j = 0; i > 0; i--, d++)
				d->e = eread(R);
			d->e = 0;	/* sentinnel expr * */
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
			rv->op = (efunc *)(size_t)k;
			rv->L.p = p;
			rv->R.e = eread(R);
			return rv;

		case 5: /* if */
			rvif = (expr_if *)mem(sizeof(expr_if));
			PSHV(rvif->val = 0;)
			rvif->op = (efunc *)(size_t)k;
			rvif->next = iflist;
			iflist = rvif;
			rvif->e = eread(R);
			rvif->T = L = eread(R);
			rvif->F = L = eread(R);
			return (expr *)rvif;

		case 11: /* OPCOUNT */
		case 6: /* sumlist */
			i = 0;
			Xscanf(R, "%d", &i);
			if (i <= 2 && (optype[k] == 6 || i < 1))
				badline(R);
			if (slmax < i)
				slmax = i;
			args = (expr**)new_mblk(htcl(i*sizeof(expr*)));
			rv = new_expr(S, k, (expr*)args, (expr*)(args+i));
			do
				*args++ = eread(R);
				while(--i > 0);
			return rv;
			}
	badline(R);
	return 0;
	}

 static void
ogfree(Static *S, ograd *og)
{
	ograd *og1;

	if ((og1 = og)) {
#ifdef DEBUG
		if (ogzork) {
			do if (og == ogzork)
				printf("");
			   while((og = og->next));
			og = og1;
			}
#endif
		while(og1->next)
			og1 = og1->next;
		og1->next = freeog;
		freeog = og;
		}
	}

 static real
ogfree1(Static *S, ograd **ogp)
{
	ograd *og = *ogp;
	*ogp = og->next;
	og->next = freeog;
	freeog = og;
	return og->coef;
	}

 static ograd *
new_ograd(Static *S, ograd *next, int varno, real coef)
{
	ograd *rv;

	if ((rv = freeog))
		freeog = rv->next;
	else
		rv = (ograd *)mem_ASL(S->a, sizeof(ograd));
	rv->next = next;
	rv->varno = varno;
	rv->coef = coef;
#ifdef DEBUG
	if (rv == ogzork)
		printf("");
#endif
	return rv;
	}

 static linarg *
lahash(Static *S, linarg *la)
{
	ASLTYPE *asl;
	Ulong x;
	int k, nnz;
	linarg *la1, *la2, **lap, **lap1, **q, **q0, **qe;
	ograd *og, *og1;
	size_t mask;
	union {real r; Ulong L[2];} u;

	asl = S->asl;
	x = nnz = la->nnz;
	for(og = la->nz; og; og = og->next) {
		u.r = og->coef;
		x = (x << 1 | (x & 0x80000000) >> 31) ^
			  (og->varno*101 + u.L[W0] + 257*u.L[W1]);
		}
	for(lap = &lthash[x & lthashmask]; (la1 = *lap); lap = &la1->hnext)
		if (la1->nnz == nnz) {
			og = la->nz;
			og1 = la1->nz;
			for(;;) {
				if (!og) {
					if (!og1)
						return la1;
					break;
					}
				if (!og1)
					break;
				if (og->varno != og1->varno
				 || og->coef != og1->coef)
					break;
				og = og->next;
				og1 = og1->next;
				}
			}
	*lap = la;
	if (++S->asl->P.nlttot > lthashmask) {
		mask = lthashmask;
		k = S->klthash;
		q = q0 = lthash;
		qe = q + mask;
		lthashmask = mask = (mask << 1) | 1;
		lap = lthash = (linarg**)new_mblkzap(asl, S->klthash = k + 1);
		while(q <= qe) {
			for(la2 = *q++; la2; la2 = la1) {
				la1 = la2->hnext;
				x = la2->nnz;
				for(og = la2->nz; og; og = og->next) {
					u.r = og->coef;
					x = (x << 1 | (x & 0x80000000) >> 31) ^
						  (og->varno*101 + u.L[W0] + 257*u.L[W1]);
					}
				lap1 = lap + (x & mask);
				la2->hnext = *lap1;
				*lap1 = la2;
				}
			}
		}
	return la;
	}

 static range *
new_range(Static *S, range *r, range **rp)
{
	ASLTYPE *asl;
	int k, len, uilen;
	range *r1, *r2, *r3, **rp1, **rq, **rq0, **rqe;
	size_t mask;

	asl = S->asl;
	uilen = r->nv*sizeof(int);
	len = sizeof(range) + uilen;
	r1 = (range*)mem(len);
	r1->nintv = 0;
	r1->n = r->n;
	r1->nv = r->nv;
	r1->refs = 0;
	r1->lasttermno = -1;
	r1->hnext = r1->hunext = 0;
	if (uilen)
		memcpy(r1->ui = (int*)(r1+1), r->ui, uilen);
	r1->lap = (linarg**)new_mblk(htcl(len = r->n*sizeof(linarg*)));
	memcpy(r1->lap, r->lap, len);
	r2 = r1->rlist.next = asl->P.rlist.next;
	r1->rlist.prev = r2->rlist.prev;
	r2->rlist.prev = asl->P.rlist.next = r1;
	*rp = r1;
	if (++nrange > rangehashmask) {
		mask = rangehashmask;
		k = S->krangehash;
		rq = rq0 = rangehash;
		rqe = rq + mask;
		rangehashmask = mask = (mask << 1) | 1;
		rp = rangehash = (range**)new_mblkzap(asl, S->krangehash = k + 1);
		while(rq <= rqe) {
			for(r2 = *rq++; r2; r2 = r3) {
				r3 = r2->hunext;
				rp1 = rp + (r2->chksum & mask);
				r2->hunext = *rp1;
				*rp1 = r2;
				}
			}
		del_mblk(k, rq0);
		}
	return r1;
	}

 static int
lacompar(const void *a, const void *b, void *v)
{
	ograd *oa, *ob;
	int i;
	real t;

	if (a == b)
		return 0;
	Not_Used(v);
	oa = (*(linarg **)a)->nz;
	ob = (*(linarg **)b)->nz;
	for(;;) {
		if (!oa) {
			if (ob)
				return -1;
			break;
			}
		if (!ob)
			return 1;
		if ((i = oa->varno - ob->varno))
			return i;
		if ((t = oa->coef - ob->coef))
			return t > 0. ? 1 : -1;
		oa = oa->next;
		ob = ob->next;
		}
	return 0;
	}

 static int
ndiff(range *r, range *r1)
{
	/* Return 0 if r and r1 have the same linargs; else return 1. */

	linarg **la, **la1, **la1e, **lae;

	la = r->lap;
	lae = la + r->n;
	la1 = r1->lap;
	la1e = la1 + r1->n;
	for(;;++la, ++la1) {
		if (la >= lae)
			return la1 < la1e;
		if (la1 >= la1e)
			return 1;
		if (lacompar((char*)la, (char *)la1, NULL))
			return 1;
		}
	return 0;
	}

 static range *
uhash(Static *S, range *r)
{
	ASLTYPE *asl;
	int len, n, nv, *ui, *uie;
	range *r1, **rp;
	size_t L;

	asl = S->asl;
	L = 0;
	nv = r->nv;
	ui = r->ui;
	uie = ui + nv;
	len = sizeof(int)*nv;
	while (ui < uie)
		L = 37*L + *ui++;
	r->chksum = L;
	ui = r->ui;
	rp = rangehash + (L & rangehashmask);
	n = r->n;
	if (asl->P.merge)
	    while((r1 = *rp)) {
		if (r1->nv == nv && r1->n == n && !memcmp(ui, r1->ui, len)) {
			if (!memcmp(r->lap, r1->lap, n*sizeof(linarg*)))
				return *rp;
			if (!ndiff(r, r1))
				return *rp;
			}
		rp = &r1->hunext;
		}
	return new_range(S, r, rp);
	}

 static range *
rhash(Static *S, range *r, int addnew)
{
	ASLTYPE *asl;
	int len, n;
	linarg **lae, **lap;
	ograd *og;
	range *r1, **rp;
	unsigned long L;

	asl = S->asl;
	L = 0;
	lap = r->lap;
	n = r->n;
	lae = lap + n;
	len = n * sizeof(linarg*);
	while(lap < lae) {
		L *= 37;
		for(og = (*lap++)->nz; og; og = og->next)
			L = 101*L + og->varno;
		}
	r->chksum = L;
	lap = r->lap;
	rp = rangehash + (L & rangehashmask);
	if (asl->P.merge)
	    while((r1 = *rp)) {
		if (r1->n == n && !memcmp(lap, r1->lap, len))
			return r1;
		rp = &r1->hnext;
		}
	if (addnew)
		return new_range(S, r, rp);
	return *rp;
	}

 static int
compar(const void *a0, const void *b0, void *v)
{
	int a, b, c;
	Static *S = (Static *)v;

	if ((a = *(int*)a0) >= max_var) {
		a = zl[a-nv0x];
		if ((b = *(int*)b0) >= max_var) {
			if ((c = a - zl[b-nv0x]))
				return c;
			return *(int*)a0 - b;
			}
		if (a == b)
			return -1;
		}
	else if ((b = *(int*)b0) >= max_var) {
		b = zl[b-nv0x];
		if (a == b)
			return 1;
		}
	return a - b;
	}

 static int
hscompar(const void *a0, const void *b0, void *v)
{
	int a, b, c;
	Static *S = (Static *)v;

	if ((a = *(int*)a0) >= Ncom) {
		a = zl[a];
		if ((b = *(int*)b0) >= Ncom) {
			if ((c = a - zl[b]))
				return c;
			return *(int*)a0 - b;
			}
		a -= nv0x;
		if (a == b)
			return -1;
		}
	else if ((b = *(int*)b0) >= Ncom) {
		b = zl[b] - nv0x;
		if (a == b)
			return 1;
		}
	return a - b;
	}

 static void
zcsort(Static *S, int *c, int *ci, int i, int n, int p)
{
	int j;

	if (n < nzclim || p < 0)
		qsortv(ci, n, sizeof(int), compar, S);
	else for(j = 0; i < p; i++)
		if (c[i])
			ci[j++] = i;
	}

 static ograd*
compress(Static *S, ograd *og, real *cp, int *comvar)
{
	ograd *og1;
	int i, ix, j, j1, k, nzc1, *zc1, *zci1;
	real c, t;

	c = og->varno < 0 ? ogfree1(S,&og) : 0.;
	ix = nzc1 = 0;
	zc1 = S->zc1;
	zci1 = S->zci1;
	for(og1 = og; og1; og1 = og1->next) {
		zc1[i = og1->varno] = 1;
		zci1[nzc1++] = i;
		rnz[i] = og1->coef;
		if (ix < i)
			ix = i;
		}
	if (ix < nv0x) {
		*cp = c;
		*comvar = 0;
		for(og1 = og; og1; og1 = og1->next)
			zc1[og1->varno] = 0;
		return og;
		}
	*comvar = 1;
	for(i = 0; i < nzc1; )
		if ((j = zci1[i]) < nv0x)
			i++;
		else {
			if (!zc[j]++)
				zci[nzc++] = j;
			t = rnz[j];
			j1 = j - nv0x;
			if ((og1 = (S->asl->P.dv + j1)->ll)) {
				if (og1->varno < 0) {
					c += t*og1->coef;
					og1 = og1->next;
					}
				for(; og1; og1 = og1->next)
					if (!zc1[k = og1->varno]++) {
						zci1[nzc1++] = k;
						rnz[k] = t*og1->coef;
						}
					else
						rnz[k] += t*og1->coef;
				}
			zc1[j] = 0;
			zci1[i] = zci1[--nzc1];
			}
	*cp = c;
	ogfree(S, og);
	og = 0;
	if (nzc1 > 0) {
		zcsort(S, zc1, zci1, 0, nzc1, max_var);
		while(nzc1 > 0) {
			i = zci1[--nzc1];
			zc1[i] = 0;
			if ((t = rnz[i])) {
				og = new_ograd(S, og, i, t);
				if (!zc[i]++)
					zci[nzc++] = i;
				}
			}
		}
	return og;
	}

 static linarg *
afree(Static *S, ograd *og, expr **ep)
{
	linarg *la, *la1, *rv;
	real c, s, t;
	ograd *og1, *ogx;
	int comvar, nnz;
	la_ref *refs;
	ASLTYPE *asl = S->asl;

	rv = 0;
	if (!og || !(og = compress(S,og,&c,&comvar)))
		goto done;

	if ((la = ltfree))
		ltfree = la->hnext;
	else {
		la = (linarg *)mem(sizeof(linarg));
		la->refs = 0;
		}

	s = og->coef;
	if (s < 0)
		s = -s;
	la->nz = ogx = og1 = og;
	nnz = 1;
	while((og1 = og1->next)) {
		nnz++;
		t = og1->coef;
		if (t < 0)
			t = -t;
		if (s < t) {
			s = t;
			ogx = og1;
			}
		}

	la->nnz = nnz;
	if ((s = ogx->coef) != 1.)
		for(og1 = og; og1; og1 = og1->next)
			og1->coef /= s;
	lt_scale = s;
	la1 = lahash(S, la);
	if (la1 == la) {
		la->refs = 0;
		la->v = 0;
		la->termno = Termno;
		la->tnext = tlist;
		tlist = la;
		la->lnext = asl->P.lalist;
		asl->P.lalist = la;
		la->hnext = 0;
		}
	else {
		/* !!? Eventually create or adjust la1->v to save cycles. */
		/* This requires possible adjustment for differing */
		/* constants and scale factors. */
		if (la1->termno == Termno)
			asl->P.ndupst++;	/* research statistic */
		else {
			free_laref(S, &la->refs);
			la1->termno = Termno;
			la1->tnext = tlist;
			tlist = la1;
			asl->P.ndupdt++;	/* research statistic */
			}
		ogfree(S, og);
		la->hnext = ltfree;
		ltfree = la;
		la = la1;
		}
	if (ep && (nnz > 1 || comvar)) {
		la->refs = refs = new_laref(S, la->refs);
		refs->ep = ep;
		refs->c = c;
		refs->scale = s;
		}
	if (nnz > 1)
		rv = la;
 done:
	return rv;
	}

 static ograd *
af_sum(Static *S, ograd *Log, ograd *Rog)
{
	ograd *oL, *oR, *oR1, *og, **ogp;

	oL = Log;
	oR = Rog;
	ogp = &og;
	for(;;) {
		if (!oL) {
			*ogp = oR;
			break;
			}
		if (!oR) {
			*ogp = oL;
			break;
			}
		if (oL->varno > oR->varno) {
			*ogp = oR;
			ogp = &oR->next;
			oR = *ogp;
			}
		else {
			if (oL->varno == oR->varno) {
				oL->coef += oR->coef;
				oR1 = oR->next;
				oR->next = 0;
				ogfree(S, oR);
				oR = oR1;
				if (oL->coef == 0.) {
					oR1 = oL->next;
					oL->next = 0;
					ogfree(S, oL);
					oL = oR1;
					continue;
					}
				}
			*ogp = oL;
			ogp = &oL->next;
			oL = *ogp;
			}
		}
	return og;
	}

 static cgrad *
cf_sum(Static *S, cgrad *Lcg, ograd *Rog)
{
	cgrad *cg, **cgp, *oL;
	ograd *oR, *oR1;

	oL = Lcg;
	oR = Rog;
	cgp = &cg;
	cg = 0;
	for(;;) {
		if (!oL) {
			assert(!oR);
			break;
			}
		if (!oR) {
			*cgp = oL;
			break;
			}
		assert(oL->varno <= oR->varno);
		if (oL->varno == oR->varno) {
			oL->coef += oR->coef;
			oR1 = oR->next;
			oR->next = 0;
			ogfree(S, oR);
			oR = oR1;
			}
		*cgp = oL;
		cgp = &oL->next;
		oL = *cgp;
		}
	return cg;
	}

 static void
sumlist_afree(Static *S, ograd *Laf, expr *e, expr **argso, expr **sls, expr **slscr)
{
	int nlin, nnl;
	expr *e1, *e2, **ep, **ep1;
	ASL *asl = S->a;

	nlin = sls - slscr;
	ep = e->L.ep;
	nnl = argso - ep;
	switch(nlin) {
	  case 1:
		e1 = slscr[0];
		break;
	  case 2:
		e1 = new_expr(S, f_OPPLUS, slscr[0], slscr[1]);
		e1->dL = e1->dR = 1.;
		break;
	  default:
		ep1 = (expr**)new_mblk(htcl(nlin*sizeof(expr*)));
		e1 = new_expr(S, f_OPSUMLIST, (expr*)ep1, (expr*)(ep1+nlin));
		memcpy(ep1, slscr, nlin*sizeof(expr*));
		}
	switch(nnl) {
	  case 1:
		e2 = *ep;
		break;
	  case 2:
		e2 = new_expr(S, f_OPPLUS, ep[0], ep[1]);
		e2->dL = e2->dR = 1.;
		break;
	  default:
		ep1 = (expr**)new_mblk(htcl(nnl*sizeof(expr*)));
		e2 = new_expr(S, f_OPSUMLIST, (expr*)ep1, (expr*)(ep1+nnl));
		memcpy(ep1, ep, nnl*sizeof(expr*));
		}
	del_mblk(htcl((e->R.ep - ep)*sizeof(expr*)), ep);
	e->op = (efunc*)f_OPPLUS;
	e->dL = e->dR = 1.;
	e->L.e = e1;
	e->R.e = e2;
	afree(S, Laf, &e->L.e);
	}

 static void
sdvcite(Static *S, int k)
{
	range *r;
	split_ce *cs;
	linarg *la, **lap, **lape;

	cs = S->asl->P.Split_ce + (k - max_var);
	r = cs->r;
	lap = r->lap;
	lape = lap + r->n;
	while(lap < lape) {
		la = *lap++;
		if (la->termno != Termno) {
			free_laref(S, &la->refs);
			la->termno = Termno;
			la->tnext = tlist;
			tlist = la;
			}
		}
	}

 static ograd *
awalk(Static *S, expr *e)		/* return 0 if e is not linear */
{
	ASLTYPE *asl;
	int k, k1, kscr;
	expr *L, **args, **argse, **argso, **sls, **slscr;
	expr_va *rva;
	de *d;
	expr_if *rvif;
	expr_f *rvf;
	ograd *Laf, *Raf, *rv, *taf;
	linarg *la, **nl;
	argpair *ap, *ape;
	real t;

	switch(optypeb[k = Intcast e->op]) {

		case 1:	/* unary */
			Laf = awalk(S, e->L.e);
			if (k == f_OPUMINUS) {
				if ((taf = Laf)) {
					do taf->coef = -taf->coef;
						while((taf = taf->next));
					return Laf;
					}
				}
			if (Laf)
				afree(S, Laf, &e->L.e);
			break;

		case 2:	/* binary */
			Laf = awalk(S, e->L.e);
			Raf = awalk(S, e->R.e);
			if (Laf && Raf)
			 switch(k) {
			  case f_OPMINUS:
				taf = Raf;
				do taf->coef = -taf->coef;
					while((taf = taf->next));
				/* no break; */
			  case f_OPPLUS:
				return af_sum(S, Laf, Raf);

			  case f_OPMULT:
				if (Raf->varno < 0 && !Raf->next) {
 swap:
					taf = Laf;
					Laf = Raf;
					Raf = taf;
					goto scale;
					}
				if (Laf->varno < 0 && !Laf->next) {
					taf = Raf;
 scale:					t = Laf->coef;
					if (t == 0.) {
						ogfree(S, Raf);
						return Laf;
						}
					do taf->coef *= t;
						while((taf = taf->next));
					ogfree(S, Laf);
					return Raf;
					}
				break;

			  case f_OPDIV:
				if (Raf->varno < 0 && !Raf->next) {
					Raf->coef = 1. / Raf->coef;
					goto swap;
					}
			  }
			afree(S, Laf, &e->L.e);
			afree(S, Raf, &e->R.e);
			break;

		case 3:	/* vararg (min, max) */
			rva = (expr_va *)e;
			for(d = rva->L.d; (L = d->e); d++)
				afree(S, awalk(S,L), &d->e);
			break;

		case 4: /* piece-wise linear */
			afree(S, awalk(S,e->R.e), &e->R.e);
			break;

		case 5: /* if */
			rvif = (expr_if *)e;
			afree(S, awalk(S,rvif->e), &rvif->e);
			/* The above may introduce "range" rows that */
			/* are not involved in any derivatives... */
			afree(S, awalk(S,rvif->T), &rvif->T);
			afree(S, awalk(S,rvif->F), &rvif->F);
			break;

		case 6: /* sumlist */
			args = e->L.ep;
			argse = e->R.ep;
			k1 = argse - args;
			while(!(Laf = awalk(S,*args++)))
				if (args >= argse)
					return 0;
			kscr = -1;
			asl = S->asl;
			if (!slscratchlev++)
				slscr = slscratch;
			else {
				kscr = htcl(k1*sizeof(expr*));
				slscr = (expr**)new_mblk(kscr);
				}
			sls = slscr;
			argso = args - 1;
			*sls++ = *argso;
			while(args < argse)
				if ((Raf = awalk(S,*args))) {
					*sls++ = *args++;
					Laf = af_sum(S, Laf, Raf);
					}
				else
					*argso++ = *args++;
			rv = 0;
			if (argso == e->L.ep) {
				rv = Laf;
				goto delscratch;
				}
			sumlist_afree(S, Laf, e, argso, sls, slscr);
 delscratch:
			--slscratchlev;
			if (kscr >= 0)
				del_mblk(kscr, slscr);
			return rv;

		case 7: /* function call */
			rvf = (expr_f *)e;
			ap = rvf->ap;
			for(ape = rvf->ape; ap < ape; ap++)
				afree(S, awalk(S,ap->e), &ap->e);

		case 8: /* OPHOL (Hollerith) */
			break;

		case 9: /* OPNUM */
			return new_ograd(S, 0, -1, ((expr_n*)e)->v);

		case 10:/* OPVARVAL */
			asl = S->asl;
			k = (expr_v *)e - var_e;
			if ((k < 0 || k >= max_var)
			 && (k = ((expr_vx*)e)->a0) < 0) {
				if ((la = ((expr_vx*)e)->la)
				  && la->termno != Termno) {
					la->termno = Termno;
					la->tnext = tlist;
					tlist = la;
					}
				return 0;
				}
			if ((k1 = (k - nv0x)) < 0) {
				if (!zc[k]++)
					zci[nzc++] = k;
				goto ogret;
				}
			if (k1 >= Ncom) {
				if (!zc[k = ((expr_vx*)e)->a1]++)
					zci[nzc++] = k;
				sdvcite(S,k);
				return 0;
				}
			if ((nl = (asl->P.dv + k1)->nl)) {
				if (!zc[k]++)
					zci[nzc++] = k;
				while((la = *nl++))
					if (la->termno != Termno) {
						la->termno = Termno;
						la->tnext = tlist;
						tlist = la;
						}
				return 0;
				}
			if (!zc[k]++)
				zci[nzc++] = k;
		ogret:
			return new_ograd(S, 0, k, 1.);

		case 11: /* OPCOUNT */
			args = e->L.ep;
			for(argse = e->R.ep; args < argse; args++)
				afree(S, awalk(S, *args), args);
			break;

		default:
			scream(S->R, ASL_readerr_bug,
				"awalk: unexpected optype[%d] = %d\n",
				k, optype[k]);
			}
	return 0;
	}

 static void
cexp_upgrade(Static *S, int t)
{
	int k, n1, n2;
	cexp *ce;
	int *z;
	expr_v **vp;
	split_ce *cs;
	ASLTYPE *asl = S->asl;

	n2 = t - Ncom;
	k = htcl(t*(sizeof(cexp) + sizeof(int) + sizeof(expr_v*))
		+ n2*sizeof(split_ce));
	memset(ce = (cexp *)new_mblk(k), 0, n1 = sizeof(char*) << k);
	n1 = (n1 + Ncom*sizeof(split_ce))
		/ (sizeof(cexp) + sizeof(int) + sizeof(split_ce)
			+ sizeof(expr_v*));
	n2 = n1 - Ncom;
	cs = (split_ce *)(ce + n1);
	vp = (expr_v **)(cs + n2);
	z = (int *)(vp + n1);
	if (cexps) {
		if (nsce)
			memcpy(cs, asl->P.Split_ce, nsce*sizeof(split_ce));
		memcpy(ce, cexps, cexp_n*sizeof(cexp));
		memcpy(z, zl, cexp_n*sizeof(int));
		memcpy(vp, varp, cexp_n*sizeof(expr_v*));
		del_mblk(cexp_k, cexps);
		}
	nsce = n2;
	asl->P.Split_ce = cs;
	cexps = ce;
	zl = z;
	cexp_k = k;
	cexp_n = n1;
	varp = vp;
	}

 static void
zc_upgrade(Static *S)
{
	int k, n, n0;
	int *z;
	ASLTYPE *asl = S->asl;

	k = htcl(sizeof(int)*(max_var1 + 1)) + 1;
	z = (int*)new_mblk(k);
	n = (sizeof(char*)/sizeof(int)) << (k-1);
	memset(z + n, 0, n*sizeof(int));
	if (zci) {
		n0 = (sizeof(char*)/sizeof(int)) << (kzc - 1);
		memcpy(z, zci, n0*sizeof(int));
		memcpy(z+n, zci+n0, n0*sizeof(int));
		del_mblk(kzc, zci);
		}
	kzc = k;
	zci = z;
	zc = z + n + 1;
	zc_lim = n;
	}

 static void
la_replace(Static *S, linarg *la)
{
	la_ref *r;
	expr_v *v;
	expr_vx *vx;
	expr *eout;
	ASLTYPE *asl = S->asl;

	if (la->nnz > 1) {
		if (!(v = la->v)) {
			vx = (expr_vx *)mem(sizeof(expr_vx));
			vx->la = la;
			vx->a0 = vx->a1 = -1;
			la->v = v = (expr_v *)vx;
			v->a = max_var + nndv++;
			lasta00++;
			v->op = (efunc *)f_OPVARVAL;
			if (larvlist) {
				*larvlist = v;
				larvlist = (expr_v**)&v->v;
				}
			}
		}
	else
		v = var_e + la->nz->varno;
	for(r = la->refs; r; r = r->next) {
		efree(S, *r->ep);
		eout = (expr *)v;
		if (r->scale != 1.) {
			if (r->scale == -1.) {
				eout = new_expr(S, f_OPUMINUS, eout, 0);
				eout->dL = -1.;
				}
			else
				eout = new_expr(S, f_OPMULT, eout,
					(expr*)new_expr_n(S, r->scale));
			}
		if (r->c) {
			eout = new_expr(S, f_OPPLUS, eout,
					(expr*)new_expr_n(S, r->c));
			eout->dL = 1.;
			}
		*r->ep = eout;
		}
	free_laref(S, &la->refs);
	}

 static void
tlistgen(Static *S, ps_func *f)
{
	int *ci, *cie, i, t;
	linarg *la, **lap, **lape, *tl;
	ograd *og;
	range *r;
	psb_elem *b, *be;

	t = f->nb;
	b = f->b;
	be = b + t;
	t = ++Termno;
	tl = 0;
	for(; b < be; b++) {
		if ((ci = b->ce)) {
			cie = ci + *ci;
			do {
				if (!zc[i = nv0x + *++ci])
					zci[nzc++] = i;
				}
				while(ci < cie);
			}
		r = b->U;
		lap = r->lap;
		lape = lap + r->n;
		while(lap < lape) {
			la = *lap++;
			if (la->termno != t) {
				la->termno = t;
				la->tnext = tl;
				tl = la;
				for(og = la->nz; og; og = og->next)
					if (!zc[og->varno]++)
						zci[nzc++] = og->varno;
				}
			}
		}
	tlist = tl;
	}

 static int
might_expand(Static *S, expr *e)
{
 top:
	switch(Intcast e->op) {
		case f_OPPLUS:
		case f_OPMINUS:
		case f_OPUMINUS:
		case f_OPSUMLIST:
			return 1;
		case f_OPMULT:
			if (Intcast e->R.e->op == f_OPNUM) {
				e = e->L.e;
				goto top;
				}
			if (Intcast e->L.e->op == f_OPNUM) {
				e = e->R.e;
				goto top;
				}
			break;
		case f_OPVARVAL:
			if (((expr_v*)e)->a >= nv0x)
				return 1;
		}
	return 0;
	}

 static void
ce_split(Static *S, int i, ps_func *f)
{
	int j, j1, je, k, n;
	cexp *c, *ce;
	expr **ep;
	expr_v **vp, **vp0;
	expr_vx *vx;
	split_ce *cs;
	psb_elem *b;
	ASLTYPE *asl = S->asl;

	asl->P.ndvsplit++;
	n = f->nb;
	j1 = asl->P.ndvspout;
	j = k = j1 + Ncom;
	asl->P.ndvspout += n;
	je = j + n;
	if (je > cexp_n)
		cexp_upgrade(S, je);
	c = cexps + j;
	b = f->b;
	cs = asl->P.Split_ce + (j-Ncom);
	for(ce = c + n; c < ce; c++, b++, cs++) {
		c->e = b->D.e;
		cs->r = b->U;
		cs->ce = b->ce;
		}
	c = cexps + i;
	vp0 = vp = varp + j;
	j = max_var + nndv;
	je = j + n;
	nndv += n;
	j1 += max_var;
	while(j < je) {
		vx = (expr_vx *)mem(sizeof(expr_vx));
		vx->la = 0;
		*vp++ = (expr_v*)vx;
		vx->v.a = vx->a0 = j++;
		vx->a1 = j1++;
		vx->v.op = (efunc *)f_OPVARVAL;
		}
	if (n == 2)
		c->e = new_expr(S, f_OPPLUS, (expr*)vp0[0], (expr*)vp0[1]);
	else {
		ep = (expr**)new_mblk(htcl(n*sizeof(expr*)));
		memcpy(ep, vp0, n*sizeof(expr*));
		c->e = new_expr(S, f_OPSUMLIST, (expr*)ep, (expr*)(ep+n));
		}
	if ((max_var1 += n) >= zc_lim)
		zc_upgrade(S);

	/* adjust zl for use in compar() */
	i += nv0x;
	j = k + nv0x;
	n += k;
	do {
		zl[k++] = i;
		zci[nzc++] = j++;	/* for crefs */
		}
		while(k < n);
	}

 static void
linpart_augment(Static *S, cexp *c, ograd *og, ps_func *f)
{
	linpart *L, *Le;
	int n;
	ograd *og1;
	real t;
	ASL *asl = S->a;

	if (og->varno < 0) {
		if ((t = ogfree1(S, &og)))
			f->b->D.e = new_expr(S, f_OPPLUS, f->b->D.e,
					(expr*)new_expr_n(S, t));
		if (!og)
			return;
		}
	if ((n = c->nlin)) {
		L = c->L;
		Le = L + n;
		og1 = 0;
		do {
			--Le;
			og1 = new_ograd(S, og1, Le->v.i, Le->fac);
			}
			while(Le > L);
		del_mblk(htcl(htcl(n*sizeof(linpart))), c->L);
		og = af_sum(S, og, og1);
		}
	og1 = og;
	for(n = 0; og; og = og->next)
		n++;
	c->L = L = (linpart*)new_mblk(htcl((c->nlin = n)*sizeof(linpart)));
	for(og = og1; og; og = og->next, L++) {
		L->v.i = og->varno;
		L->fac = og->coef;
		}
	ogfree(S, og1);
	}

 static ograd *
linterms(Static *S, cexp *c, real scale)
{
	linpart *L = c->L;
	linpart *Le = L + c->nlin;
	ograd *og = 0;
	while(Le > L) {
		--Le;
		og = new_ograd(S, og, Le->v.i, scale*Le->fac);
		}
	return og;
	}

 static void
dvwalk(Static *S, int i)
{
	int n;
	ograd *og, *og0;
	linarg *tl, **tp;
	ps_func f;
	dv_info *dvi;
	cexp *c;
	int split;
	ASLTYPE *asl = S->asl;

	Termno++;
	tlist = 0;
	c = cexps + i;
	split = zl[i] & 2;
	zl[i] = 0;
	if (split && might_expand(S, c->e)) {
		asl->P.ndvspcand++;
		if ((og = cotermwalk(S, &c->e, &f, 0, 0)) && f.nb >= 1)
			linpart_augment(S, c, og, &f);
		if (f.nb > 1) {
			zl[i] = f.nb;
			ce_split(S, i, &f);
			c = cexps + i;
			og = 0;
			}
		else if (f.nb == 1) {
			c->e = f.b->D.e;
			og = 0;
			}
		tlistgen(S, &f);
		del_Elemtemp(S, last_psb_elem);
		}
	else
		og = awalk(S, c->e);
	og0 = og;
	if (c->nlin)
		og = af_sum(S, og, linterms(S,c,1.));
	asl->P.dvsp0[i+1] = asl->P.ndvspout + Ncom;
	dvi = asl->P.dv + i;
	if (og0) {
		c->cref = crefs(S);
		dvi->ll = og;
		dvi->nl = 0;
		return;
		}
	if ((dvi->lt = afree(S, og, 0)))
		dvi->scale = lt_scale;
	for(n = 1, tl = tlist; tl; tl = tl->tnext)
		n++;
	dvi->ll = 0;
	dvi->nl = tp = (linarg **)mem(n*sizeof(linarg*));
	for(tl = tlist; tl; tl = tl->tnext)
		la_replace(S, *tp++ = tl);
	*tp = 0;
	c->cref = crefs(S);
	}

 static int
colindvref(Static *S, expr *e, int ndv)
{
	expr **a, **ae;
	int j, k, rv = 0;

 top:
	switch(Intcast e->op) {
		case f_OPPLUS:
		case f_OPMINUS:
			rv |= colindvref(S, e->R.e, ndv);
			/* no break */
		case f_OPUMINUS:
			e = e->L.e;
			goto top;
		case f_OPSUMLIST:
			a = e->L.ep;
			ae = e->R.ep;
			while(a < ae)
				rv |= colindvref(S, *a++, ndv);
			break;
		case f_OPVARVAL:
			k = ((expr_v *)e)->a;
			if ((k -= nv0x) < 0)
				break;
			if (zl[k]) {
				rv |= zl[k];
				break;
				}
			zl[k] = 1;
			if ((j = colindvref(S, (S->cexps + k)->e, k))) {
				rv |= j;
				zl[k] |= j;
				}
			break;
		case f_OPMULT:
			if (Intcast e->R.e->op == f_OPNUM) {
				e = e->L.e;
				goto top;
				}
			if (Intcast e->L.e->op == f_OPNUM) {
				e = e->R.e;
				goto top;
				}
			/* no break */
		default:
			if (ndv >= 0)
				rv = zl[ndv] |= 2;
		}
	return rv;
	}

 static void
dv_walk(Static *S)
{
	int i, m, n;
	expr_v *v;
	ASLTYPE *asl = S->asl;

	m = asl->i.n_con0;
	if ((n = Ncom)) {
		for(i = 0; i < n_obj; i++)
			colindvref(S, obj_de[i].e, -1);
		for(i = 0; i < m; i++)
			colindvref(S, con_de[i].e, -1);
		larvlist = &v;
		for(i = 0; i < n; i++)
			dvwalk(S, i);
		*larvlist = 0;
		if ((i = asl->P.ndvspout))
			nv0b = nv0x + Ncom + i;
		larvlist = 0;
		}
	}

 static int	/* adjust zci, return k s.t. zci[i] < nv0 for i < k */
nzcperm(Static *S)
{
	int i, j, k;

	k = nzc;
	for(i = 0; i < k; )
		if (zci[i] >= nv0x) {
			j = zci[--k];
			zci[k] = zci[i];
			zci[i] = j;
			}
		else
			i++;
	return k;
	}

 static void
sumlist_adj(ASLTYPE *asl, expr *e, expr *e1)
{
	int k, n;
	expr **ep, **ep0, **ep1;

	ep = e->R.ep;
	k = htcl((n = ep - e->L.ep)*sizeof(expr*));
	if (n == (1 << k)) {
		ep1 = (expr**)new_mblk(k+1);
		memcpy(ep1, ep0 = e->L.ep, n*sizeof(expr*));
		del_mblk(k, ep0);
		e->L.ep = ep1;
		ep = ep1 + n;
		}
	*ep++ = e1;
	e->R.ep = ep;
	}

 static void
termwalk(Static *S, expr **ep, PSfind *p)
{
	ograd *og;
	linarg **lap, **lap1;
	linarg *tl;
	int *cp, *cp0, *ce, *cee, i, j, k, kl, n, ncp, nzc2, *ui;
	size_t len;
	range *r;
	list *L;
	expr *e, *e1, *e2;
	split_ce *cs;
	psb_elem *b;
	Elemtemp *bt;
	ps_func *f;
	ASLTYPE *asl = S->asl;

	Termno++;
	tlist = 0;
	afree(S, awalk(S, *ep), ep);

	i = k = nzcperm(S);
	f = p->f;
	if (!larvlist)
		for(; i < nzc; i++)
			if ((j = zci[i]) < max_var) {
				for(L = cexps[j-nv0x].cref; L; L = L->next)
					if (!zc[L->item.i]++)
						zci[nzc++] = L->item.i;
				}
			else {
				cs = asl->P.Split_ce + (j-max_var);
				if ((ce = cs->ce)) {
					cee = ce + *ce;
					while(ce++ < cee)
						if (!zc[j = *ce + nv0x]++)
							zci[nzc++] = j;
					}
				}

	r = (range *)rnz;	/* scratch */
	if ((ncp = nzc - k)) {
		zcsort(S, zc, zci+k, nv0x, ncp, -1 /*max_var*/);
		cp = cp0 = (int*)(r+1);
		i = k;
		*cp = ncp;
		do {
			zc[j = zci[i++]] = 0;
			*++cp = j - nv0x;
			} while(i < nzc);
		}
	else
		cp0 = 0;

	for(n = 0, tl = tlist; tl; tl = tl->tnext) {
		n++;
		for(og = tl->nz; og; og = og->next)
			if (!zc[og->varno]++)
				zci[k++] = og->varno;
		}
	nzc2 = k;
	if (zc[-1])
		--nzc2;	/* ignore constant */
	if (n <= 0 && Intcast (*ep)->op == f_OPNUM)
		goto done;

	r->n = n;
	r->nv = nzc2;
	r->lap = lap1 = lap = (linarg**)
		new_mblk(kl = htcl(n*sizeof(linarg*)));
	for(tl = tlist; tl; tl = tl->tnext)
		la_replace(S, *lap1++ = tl);
	if (n > 1)
		qsortv(lap, n, sizeof(linarg*), lacompar, NULL);
	r->ui = ui = cp0 ? cp+1 : (int*)(r+1);
	zcsort(S, zc, zci, 0, nzc2, nv0x);
	for(j = 0; j < nzc2; j++)
		*ui++ = zci[j];
	r = n >= nzc2 ? uhash(S,r) : rhash(S,r,1);
	del_mblk(kl, lap);

	e1 = *ep;
	if (!r || (i = r->lasttermno) == -1 || r->lastgroupno != Groupno) {
		bt = p->b;
		if ((i = f->nb++) >= bt->nmax)
			upgrade_Elemtemp(S, bt);
		b = f->b + i;
		b->conno = Conno;
		b->termno = i;
		b->groupno = Groupno;
		if ((b->U = r)) {
			r->lasttermno = i;
			r->lastgroupno = Groupno;
			}
		b->D.e = e1;
		if (cp0) {
			cp = (int*)mem(len = (ncp+1)*sizeof(int));
			memcpy(cp, cp0, len);
			cp0 = cp;
			}
		b->ce = cp0;
		while(k > 0)
			zc[zci[--k]]= 0;
		}
	else {
		b = f->b + i;
		while(k > 0)
			zc[zci[--k]] = 0;
		if ((cp = b->ce)
		 && (*cp != ncp || memcmp(cp+1, cp0+1, ncp*sizeof(int)))) {
			/* fix up f->ce */
			n = *cp;
			while(k < n)
				zc[zci[k++] = *++cp] = 1;
			for(n = 0; n < ncp; )
				if (!zc[j = cp0[++n]]++)
					zci[k++] = j;
			qsortv(zci, k, sizeof(int), hscompar, S);
			b->ce = cp = (int*)mem((k+1)*sizeof(int));
			*cp++ = k;
			n = 0;
			while(n < k)
				zc[*cp++ = zci[n++]] = 0;
			}
		e = b->D.e;
		switch(Intcast e->op) {
		  case f_OPPLUS:
			ep = (expr**)new_mblk(2);
			ep[0] = e->L.e;
			ep[1] = e->R.e;
			ep[2] = e1;
			e->L.ep = ep;
			e->R.ep = ep + 3;
			e->op = (efunc *)f_OPSUMLIST;
			break;

		  case f_OPSUMLIST:
			sumlist_adj(asl, e, e1);
			break;
		  default:
			b->D.e = e2 = new_expr(S, 0, e, e1);
			e2->op = (efunc *)f_OPPLUS;
		  }
		}
 done:
	nzc = 0;
	}


 static void
co_finish(ps_func *f)
{
	range *r;
	psb_elem *b, *be;
	psg_elem *g, *ge;

	b = f->b;
	for(be = b + f->nb; b < be; b++)
		if ((r = b->U))
			r->lasttermno = -1;
	g = f->g;
	for(ge = g + f->ng; g < ge; g++)
		for(b = g->E, be = b + g->ns; b < be; b++)
			if ((r = b->U))
				r->lasttermno = -1;
	}

 static expr *
ecopy(Static *S, expr *e)
{
	expr **a, **a1, **ae;
	int n, op;

	switch(op = Intcast e->op) {
		case f_OPPLUS:
		case f_OPMINUS:
			e = new_expr(S, op, ecopy(S, e->L.e), ecopy(S, e->R.e));
			break;
		case f_OPMULT:
			if (Intcast e->L.e->op == f_OPNUM) {
				e = new_expr(S, op, ecopy(S, e->R.e), (expr*)
					new_expr_n(S, ((expr_n*)e->L.e)->v));
				break;
				}
			assert(Intcast e->R.e->op == f_OPNUM);
			e = new_expr(S, op, ecopy(S, e->L.e),
				(expr*)new_expr_n(S, ((expr_n*)e->R.e)->v));
			break;
		case f_OPUMINUS:
			e = new_expr(S, op, ecopy(S, e->L.e), 0);
			break;
		case f_OPSUMLIST:
			a = e->L.ep;
			ae = e->R.ep;
			a1 = (expr**)new_mblk_ASL(S->a, htcl((n = ae-a)*sizeof(expr*)));
			e = new_expr(S, op, (expr*)a1, (expr*)(a1+n));
			while(a < ae)
				*a1++ = ecopy(S, *a++);
			break;
		case f_OPVARVAL:
			break;
#ifdef DEBUG
		default:
			fprintf(Stderr, "Impossible case in ecopy!\n");
			exit(1);
#endif
		}
	return e;
	}

 static int getgroup ANSI((Static *S, real scale, expr *e, PSfind *p));

 static ograd *
ltermwalk(Static *S, real scale, expr **ep, PSfind *p)
{
	ASLTYPE *asl;
	EU *en;
	dv_info *dvi;
	expr **a, **ae, *e;
	int k, k0;
	ograd *og, *rv;
	cexp *ce;

	asl = S->asl;
	rv = 0;
 top:
	e = *ep;
	switch(Intcast e->op) {
		case f_OPPLUS:
			rv = af_sum(S, rv, ltermwalk(S, scale, &e->L.e, p));
			ep = &e->R.e;
			goto top;
		case f_OPMINUS:
			rv = af_sum(S, rv, ltermwalk(S, scale, &e->L.e, p));
			ep = &e->R.e;
			scale = -scale;
			goto top;
		case f_OPUMINUS:
			ep = &e->L.e;
			scale = -scale;
			goto top;
		case f_OPSUMLIST:
			a = e->L.ep;
			ae = e->R.ep;
			while(a < ae)
				rv = af_sum(S, rv, ltermwalk(S, scale, a++, p));
			break;
		case f_OPVARVAL:
			k = e->a;
			if (k < nv0x) {
				rv = af_sum(S, rv, new_ograd(S,0,k,scale));
				break;
				}
			k0 = k;
			if ((k -= nv0x) >= Ncom)
				goto dflt;
			if (zl[k]) {
				ce = cexps + k;
				*ep = ecopy(S, ce->e);
				if (ce->nlin)
					rv = af_sum(S, rv, linterms(S,ce,scale));
				goto top;
				}
			if (!zc[e->a]++)
				zci[nzc++] = e->a;
			dvi = asl->P.dv + k;
			if (dvi->nl)
				goto dflt;
			og = new_ograd(S, 0, k0, scale);
			rv = af_sum(S, rv, og);
			break;
		case f_OPNUM:
			rv = af_sum(S, rv,
				new_ograd(S, 0, -1, scale*((expr_n *)e)->v));
			break;
		case f_OPMULT:
			en = (EU *)e->R.e;
			if (Intcast en->op == f_OPNUM) {
				*ep = e->L.e;
		case_opnum:
				if (en->u.v == 0.) {
					efree(S,*ep);
					*ep = (expr*)en;
					}
				else {
					rv = af_sum(S, rv, ltermwalk(S,
							scale*en->u.v, ep, p));
					en->u.p = expr_n_free;
					expr_n_free = en;
					}
				e->L.e = expr_free;
				expr_free = e;
				break;
				}
			en = (EU *)e->L.e;
			if (Intcast en->op == f_OPNUM) {
				*ep = e->R.e;
				goto case_opnum;
				}
		default:
		dflt:
			if (p->g && getgroup(S, scale, e, p))
				break;
			if (scale != 1.)
				*ep = scale == -1.
					? new_expr(S, f_OPUMINUS, e, 0)
					: new_expr(S, f_OPMULT, e,
						(expr *)new_expr_n(S, scale));
			termwalk(S, ep, p);
			asl->P.ns0++;
		}
	return rv;
	}

 static void
PSfind_init(Static *S, ps_func *f, PSfind *psf, int wantg)
{
	f->nxval = f->nb = f->ng = 0;
	psf->f = f;
	psf->b = new_Elemtemp(S, sizeof(psb_elem), (void**)&f->b);
	if (wantg)
		psf->g = new_Elemtemp(S, sizeof(psg_elem), (void**)&f->g);
	else {
		psf->g = 0;
		f->g = 0;
		last_psb_elem = psf->b;
		}
	}

 static void
psb_copy(psb_elem *b, psb_elem *b0, int n)
{
	range *r;
	psb_elem *be;

	memcpy(b, b0, n*sizeof(psb_elem));
	for(be = b + n; b < be; b++) {
		/* *b = *b0++; */ /* DEC Alpha gcc got this wrong */
		if (b->conno != -1 && (r = b->U)) {
			b->next = r->refs;
			r->refs = b;
			}
		}
	}

 static int
getgroup(Static *S, real scale, expr *e, PSfind *p)
{
	ps_func *f, f1;
	PSfind p1;
	expr *e0, *e1, *e2;
	ograd *og, *og1;
	psb_elem *b, *be;
	psg_elem *g;
	Elemtemp *gt;
	linarg **lap, **lape;
	linpart *L;
	range *U;
	int i, nzc1, *zc1, *zci1;
	ASL *asl = S->a;

	for(e1 = e; optype[Intcast e1->op] == 1; e1 = e1->L.e)
		e0 = e1;
	if (e1 == e || !might_expand(S, e1))
		return 0;
	PSfind_init(S, &f1, &p1, 0);
	f = p->f;
	gt = p->g;
	if ((i = f->ng++) >= gt->nmax)
		upgrade_Elemtemp(S, gt);
	Groupno = f->ng;
	memset(g = f->g + i, 0, sizeof(psg_elem));
	g->g = e;
	g->ge = e0;
	g->scale = scale;
	if ((og = ltermwalk(S, 1., &e0->L.e, &p1)))
		og = compress(S, og, &g->g0, &i);
	for(e1 = e; e1 != e0; e1 = e2) {
		e2 = e1->L.e;
		e2->R.e = e1;	/* back pointer, used in psderprop */
		}
	zc1 = S->zc1;
	zci1 = S->zci1;
	Groupno = nzc1 = 0;
	if ((og1 = og)) {
		for(i = 1; (og = og->next); i++);
		g->nlin = i;
		g->L = L = (linpart*)mem(i*sizeof(linpart));
		for(og = og1;; L++) {
			zc1[zci1[nzc1++] = L->v.i = og->varno]= 1;
			L->fac = og->coef;
			if (!(og = og->next))
				break;
			}
		ogfree(S, og1);
		}
	g->esum.op = (efunc_n *)f_OPNUM;
	g->ns = f1.nb;
	g->E = b = (psb_elem*)mem(f1.nb * sizeof(psb_elem));
	psb_copy(b, f1.b, f1.nb);
	del_Elemtemp(S, p1.b);
	for(be = b + f1.nb; b < be; b++) {
		if (!(U = b->U))
			continue;
		lap = U->lap;
		lape = lap + U->n;
		while(lap < lape)
			for(og = (*lap++)->nz; og; og = og->next)
				if (!zc1[og->varno]++)
					zci1[nzc1++] = og->varno;
		}
	zcsort(S, zc1, zci1, 0, nzc1, nv0x);
	og = 0;
	while(nzc1 > 0) {
		og = new_ograd(S, og, i = zci1[--nzc1], 0.);
		zc1[i] = 0;
		}
	g->og = og;
	return 1;
	}

 static ograd *
cotermwalk(Static *S, expr **ep, ps_func *f, int wantg, int omitdv)
{
	int comvar, n, x;
	ograd *og;
	real t;
	PSfind psf;
	psb_elem *b;
	psg_elem *g, *g1, *ge;

	PSfind_init(S, f, &psf, wantg);
	t = 0.;
	og = ltermwalk(S, 1., ep, &psf);
	if (omitdv && og)
		og = compress(S, og, &t, &comvar);
	b = f->b;
	if (f->nb + f->ng == 0)
		*ep = (expr*)new_expr_n(S, t);
	else if (t) {
		if (f->nb)
			b->D.e = new_expr(S, f_OPPLUS, b->D.e,
					(expr*)new_expr_n(S, t));
		else {
			f->nb = 1;
			memset(b, 0, sizeof(psb_elem));
			b->D.e = (expr *)new_expr_n(S, t);
			}
		}
	co_finish(f);
	if (!larvlist) {
		n = f->nb;
		if ((x = n * sizeof(psb_elem) + f->ng * sizeof(psg_elem))) {
			ASL *asl = S->a;
			g = (psg_elem*)(x >= 256 ? M1alloc(x) : mem(x));
			b = (psb_elem*)(g + f->ng);
			if (n)
				psb_copy(b, f->b, n);
			else
				b = 0;
			if ((n = f->ng)) {
				memcpy(g1 = g, f->g, n*sizeof(psg_elem));
				for(ge = g + n; g1 < ge; g1++)
					g1->ge->L.e = (expr*)&g1->esum;
				}
			else
				g = 0;
			del_Elemtemp(S, psf.b);
			if (wantg)
				del_Elemtemp(S, psf.g);
			f->b = b;
			f->g = g;
			}
		}
	return og;
	}

 static void
psfind(Static *S)
{
	ASLTYPE *asl;
	fint x;
	int i, j, k, m, mx, nv;
	linarg *la, **lap;
	ps_func *f, *f1;
	range *r, *r0;
	size_t L;
#ifdef PSHVREAD
	Ihinfo *ihi, *ihi0, *ihi1;
	int ihdlim, ihdlim1, ihdmax, ihdmin, n;
#endif
	asl = S->asl;
	m = asl->i.n_con0;
	mx = m + asl->i.nsufext[ASL_Sufkind_con];
	x = (n_obj+mx)*sizeof(ps_func)
		+ slmax*sizeof(expr*)
		+ (Ncom+1)*sizeof(int);
	asl->P.ops = f = (ps_func *)M1alloc(x);
	asl->P.cps = f1 = f + n_obj;
	slscratch = (expr**)(f1 + mx);
	asl->P.dvsp0 = (int*)(slscratch + slmax);
	*asl->P.dvsp0 = Ncom;

	k = 3;
	if ((nv = n_var) > 7) {
		nv >>= 3;
		do ++k; while(nv >>= 1);
		}
	S->klthash = S->krangehash = k;
	lthash = (linarg**)new_mblkzap(asl, k);
	rangehash = (range**)new_mblkzap(asl, k);
	L = 1;
	lthashmask = rangehashmask = (L << k) - 1;
	nrange = 0;

	/* Separate allocation of zc1, zci1, rnz to allow freeing them later. */
	/* Allow room in rnz for use as scratch in termwalk. */
	i = max_var + (sizeof(range)+sizeof(int)+sizeof(real)-1)/sizeof(real);
	rnz = (real *)M1alloc(i*(sizeof(real)+2*sizeof(int)));
	S->zc1 = (int *)(rnz + i);
	S->zci1 = S->zc1 + i;
	memset(S->zc1, 0, max_var*sizeof(int));

	Conno = -1;	/* do not record refs for split defined vars */
	dv_walk(S);
	asl->P.ndvspin = asl->P.ns0;
	asl->P.ns0 = 0;
	PSHV(n = 0);
	for(i = 0; i < n_obj; i++, f++) {
		Conno = -2 - i;
		Ograd[i] = af_sum(S, Ograd[i], cotermwalk(S, &obj_de[i].e, f,
							wantOgroups, 1));
		PSHV(n += f->ng);
		}
	PSHV(asl->P.nobjgroups = n);
	PSHV(n = 0);
	for(i = 0; i < m; i++, f1++) {
		Conno = i;
		Cgrad[i] = cf_sum(S, Cgrad[i], cotermwalk(S, &con_de[i].e, f1,
							wantCgroups, 1));
		PSHV(n += f1->ng);
		}
	PSHV(asl->P.ncongroups = n);

	/* Compute nintv value for ranges */
	/* and set lasttermno to first var. */

	r0 = (range*)&asl->P.rlist;
	for(r = asl->P.rlist.next; r != r0; r = r->rlist.next) {
		i = 0;
		if ((j = r->n) > 0) {
			lap = r->lap;
			r->lasttermno = (*lap)->nz->varno;
			while(i < j) {
				la = lap[i];
				if (la->v) {
					i++;
					la->termno = 0; /* for psgcomp */
					}
				else {
					lap[i] = lap[--j];
					lap[j] = la;
					}
				}
			}
		r->nintv = i;
		}
#ifdef PSHVREAD
	if ((ihdlim = ihd_limit) > 0) {
		ihdmin = ihdlim1 = ihdlim + 1;
		n = ihdlim1*sizeof(Ihinfo);
		asl->P.ihi = ihi0 = (Ihinfo*)M1zapalloc(n);
		ihdmax = 0;
		for(r = asl->P.rlist.next; r != r0; r = r->rlist.next)
			if ((n = r->n) > 0) {
				if (n > r->nv)
					n = r->nv;
				if (n > ihdlim)
					n = ihdlim1;
				else {
					if (ihdmax < n)
						ihdmax = n;
					if (ihdmin > n)
						ihdmin = n;
					}
				ihi = ihi0 + n - 1;
				r->rlist.prev = ihi->r;
				ihi->r = r;
				ihi->nr++;
				}
		asl->P.ihdmax = ihdmax;
		asl->P.ihdmin = asl->P.ndhmax = ihdmin;
		ihi1 = ihi = ihi0 + ihdlim;
		ihi->ihd = ihdlim1;	/* sentinel */
		for(i = ihdlim; i > 0; --i)
			if ((n = (--ihi)->nr)) {
				ihi->next = ihi1;
				ihi1 = ihi;
				ihi->ihd = i;
				ihi->k = htcl(((i*(i+1))>>1)*n*sizeof(real));
				}
		asl->P.ihi1 = ihi1;
		}
#endif
	del_mblk(S->krangehash, rangehash);
	del_mblk(S->klthash, lthash);
	lthash = 0;
	rangehash = 0;
	}

/*******/

 static void
ewalk(Static *S, expr *e, int deriv)
{
	int a0, a1, i, j, k, kf, numargs, op;
	real *b, *b0, *ra;
	unsigned int len;
#ifdef PSHVREAD
	ASL *asl;
	argpair *da;
	int i1, j0, j1, k0, k1, kn, *nn, *nn0;
	real **fh;
	expr **args1;
#endif
	expr *L, *arg, **args, **argse;
	expr_va *rva;
	de *d;
	derp *dsave;
	expr_if *rvif;
	expr_f *rvf;
	arglist *al;
	argpair *ap, *ap0, *ape;
	char *dig;

	switch(optypeb[k = Intcast e->op]) {

		case 1:	/* unary */
			e->dL = dvalue[k];	/* for UMINUS, FLOOR, CEIL */
			ewalk(S, e->L.e, deriv);
			if (deriv)
				dexpr(S, e, e->L.e, 0);
			break;

		case 2:	/* binary */
			e->dL = 1.;
			e->dR = dvalue[k];	/* for PLUS, MINUS, REM */
			if (dvalue[k] == 11)
				return;
			ewalk(S, e->L.e, deriv);
			ewalk(S, e->R.e, deriv);
			if (deriv)
				dexpr(S, e, e->L.e, e->R.e);
			break;

		case 3:	/* vararg (min, max) */
			rva = (expr_va *)e;
			rva->next2 = varg2list;
			varg2list = rva;
			if (!last_d) {
				new_derp(S, lasta, lasta, &edagread_one);
				lasta++;
				}
			rva->d0 = dsave = last_d;
#ifdef PSHVREAD
			rva->bak = last_e;
#endif
			a0 = a1 = lasta;
			j = 0;
			for(d = rva->L.d; (L = d->e); d++) {
				last_d = dsave;
				op = Intcast L->op;
				PSHV(last_e = 0;)
				ewalk(S, L, deriv);
				PSHV(d->ee = last_e;)
				if (op == f_OPNUM || L->a == noa) {
					d->d = dsave;
					d->dv.i = noa;
					}
				else if (deriv) {
					d->d = new_relo(S, L, dsave, &d->dv.i);
					j++;
					if (a1 < lasta)
						a1 = lasta;
					lasta = a0;
					}
				}
			last_d = dsave;
			if (j) {
				rva->a = lasta = a1;
				new_derp(S, 0, lasta++, &edagread_one);
				/* f_MINLIST or f_MAXLIST will replace the 0 */
				rva->R.D = last_d;
				nocopy = 1;
#ifdef PSHVREAD
				last_e = (expr *)rva;
				rva->dO.i = Hv_vararg;
#endif
				}
			else {
				rva->a = noa;
				rva->R.D = 0;
				PSHV(last_e = rva->bak);
				}
			break;

		case 4: /* piece-wise linear */
			ewalk(S, e->R.e, deriv);	/* should be a no-op */
			if (deriv) {
				new_derp(S, e->R.e->a, e->a = lasta++, &e->dL);
#ifdef PSHVREAD
				e->bak = last_e;
				last_e = e;
				e->dO.i = Hv_plterm;
#endif
				}
			break;

		case 5: /* if */
			rvif = (expr_if *)e;
			rvif->next2 = if2list;
			if2list = rvif;
			PSHV(rvif->bak = last_e;)
			ewalk(S, rvif->e, 0);
			if (deriv && !last_d) {
				new_derp(S, lasta, lasta, &edagread_one);
				lasta++;
				}
			rvif->d0 = dsave = last_d;
			a0 = lasta;
			L = rvif->T;
			op = Intcast L->op;
			PSHV(last_e = 0;)
			ewalk(S, L, deriv);
			PSHV(rvif->Te = last_e;)
			j = 0;
			if (op == f_OPNUM || L->a == noa) {
				rvif->dT = dsave;
				rvif->Tv.i = noa;
				}
			else if ((j = deriv))
				rvif->dT = new_relo(S, L, dsave, &rvif->Tv.i);
			a1 = lasta;
			lasta = a0;
			last_d = dsave;
			L = rvif->F;
			op = Intcast L->op;
			PSHV(last_e = 0;)
			ewalk(S, L, deriv);
			PSHV(rvif->Fe = last_e;)
			if (op == f_OPNUM || L->a == noa) {
				rvif->dF = dsave;
				rvif->Fv.i = noa;
				}
			else if ((j = deriv))
				rvif->dF = new_relo(S, L, dsave, &rvif->Fv.i);
			if (lasta < a1)
				lasta = a1;
			last_d = dsave;
			if (j) {
				new_derp(S, 0, rvif->a = lasta++, &edagread_one);
				rvif->D = last_d;
				nocopy = 1;
#ifdef PSHVREAD
				last_e = (expr *)rvif;
				rvif->dO.i = Hv_if;
#endif
				}
			else {
				rvif->a = noa;
				rvif->D = 0;
				PSHV(last_e = rvif->bak);
				}
			break;

		case 6: /* sumlist */
			a0 = lasta;
			j = 0;
			args = e->L.ep;
			argse = e->R.ep;
			e->a = lasta++;
			while(args < argse) {
				L = *args++;
				op = Intcast L->op;
				ewalk(S, L, deriv);
				if (op != f_OPNUM && L->a != noa && deriv) {
					new_derp(S, L->a, e->a, &edagread_one);
					j++;
					}
				}
			if (!j) {
				e->a = noa;
				lasta = a0;
				}
#ifdef PSHVREAD
			else {
				e->bak = last_e;
				last_e = e;
				e->dO.i = Hv_sumlist;
				j0 = e->R.ep - e->L.ep;
				k0 = htcl(j0*sizeof(expr*));
				j1 = j0 + j + 1;
				k1 = htcl(j1*sizeof(expr*));
				if (k1 > k0) {
					ASL *asl = S->a;
					args1 = (expr**)new_mblk(k1);
					memcpy(args1, e->L.ep,
						j0*sizeof(expr*));
					del_mblk(k0, e->L.ep);
					e->L.ep = args1;
					args = e->R.ep = args1 + j0;
					}
				argse = e->R.ep;
				args1 = e->L.ep;
				while(args1 < args) {
					arg = *args1++;
					if (arg->op != OPNUM && arg->a != noa)
						*argse++ = arg;
					}
				*argse = 0;
				}
#endif
			break;

		case 7: /* function call */
			rvf = (expr_f *)e;
			ap = rvf->sap;
			for(ape = rvf->sape; ap < ape; ap++)
				ewalk(S, ap->e,0);
			ap = rvf->ap;
			ape = rvf->ape;
			kf = ape - ap;
			al = rvf->al;
			ra = al->ra;
			numargs = al->nr;
			for(i = 0; ap < ape; ap++) {
				arg = ap->e;
				op = Intcast arg->op;
				ewalk(S, arg, deriv);
				arg->op = (efunc *)(size_t)op;
				if (arg->a != noa)
					i += deriv;
				}
			rvf->a = i ? lasta++ : noa;
			if (!numargs) {
#ifdef PSHVREAD
				rvf->da = rvf->dae = 0;
				rvf->fh= 0;	/* for debugging */
#endif
				break;
				}
			dig = 0;
			if (i) {
				b = (real *)mem_ASL(S->a, len =
#ifdef PSHVREAD
					(numargs*(numargs+3) >> 1)
#else
					numargs
#endif
						*sizeof(real));
				memset(b, 0, len);
#ifdef PSHVREAD
				al->hes = b + numargs;
#endif
				if (kf < numargs) {
					ASLTYPE *asl = S->asl;
					dig = (char*)mem(numargs);
					}
				}
			else
				b = 0;
			al->derivs = b;
			al->dig = dig;
			b0 = b;
#ifdef PSHVREAD
			asl = S->a;
			if (i) {
				da = (argpair*)mem(kf*sizeof(argpair)
					 + kf*kf*sizeof(real*));
				fh = (real**)(da + kf);
				kn = htcl(2*sizeof(int)*numargs);
				nn = (int*)new_mblk(kn);
				}
			else {
				da = 0;
				fh = 0;
				kn = 12345; /* silence false "used uninitialized" warning */
				nn = 0;
				}
			rvf->da = da;
			rvf->fh = fh;
			nn0 = nn;
#endif

			for(ap = ap0 = rvf->ap; ap < ape; ap++) {
				j = 1;
				arg = ap->e;
				arg->op = r_ops[op = Intcast arg->op];
				switch(op) {
					case f_OPNUM:
						goto loopend;
					case f_OPVARVAL:
						arg->op = (efunc *)(size_t)op;
					}
				b = b0 + (j = ap->u.v - ra);
				*b = 0;
				if (arg->a == noa)
					goto loopend;
				new_derp(S, arg->a, rvf->a, b);
#ifdef PSHVREAD
				da->e = arg;
				(da++)->u.v = b;
				*nn++ = j;
				j = 0;
#endif
			loopend:
				if (dig)
					*dig++ = j;
				}
#ifdef PSHVREAD
			rvf->dae = da;
			if (i) {
				rvf->bak = last_e;
				last_e = (expr  *)rvf;
				rvf->dO.i = Hv_func;
				b = al->hes;
				for(i1 = 0; i1 < kf; i1++) {
					i = nn0[i1];
					for(j1 = 0; j1 < kf; j1++) {
						j = nn0[j1];
						*fh++ = &b[i >= j
							? (i*(i+1)>>1)+j
							: (j*(j+1)>>1)+i];
						}
					}
				}
			if (nn0)
				del_mblk(kn, nn0);
#endif
			/* no break */

		case 8: /* OPHOL (Hollerith) */
		case 9: /* OPNUM */
			break;

		case 10:/* OPVARVAL */
			k = (expr_v *)e - S->var_e;
			if (k < 0 || k >= max_var) {
				if ((k = ((expr_vx*)e)->a1) < 0)
					return;
				}
			if (k >= 0 && deriv && !zc[k]++)
				zci[nzc++] = k;
			return;

		case 11: /* OPCOUNT */
			args = e->L.ep;
			argse = e->R.ep;
			while(args < argse)
				ewalk(S, *args++, 0);
			break;

		default:
			scream(S->R, ASL_readerr_bug,
				"unexpected optype[%d] = %d\n", k, optype[k]);
			}
#ifdef DEBUG
	if (e == opzork)
		printf("");
#endif
	e->op = r_ops[k];
	}

 static list *
crefs(Static *S)
{
	int Nv0, Nzc, i, maxvar1 = S->max_var1;
	list *rv = 0;

	Nv0 = nv0x;
	for(Nzc = nzc; Nzc > 0; ) {
		if ((i = zci[--Nzc]) >= Nv0 && i < maxvar1) {
			rv = new_list(S->a, rv);
			rv->item.i = i;
			}
		zc[i] = 0;
		}
	nzc = Nzc;
	return rv;
	}

 static funnel *
funnelfix(funnel *f)
{
	cexp *ce;
	funnel *fnext, *fprev;
	derp *d;

	for(fprev = 0; f; f = fnext) {
		fnext = f->next;
		f->next = fprev;
		fprev = f;
		ce = f->ce;
		if ((d = ce->d))
			ce->z.i = d->b.i;
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
	ASL *asl;

	if (!(d = d0))
		return dnext;
	asl = S->a;
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
			while((il = il->next2) != ile);
		}
	if (varg2list != varg2list_end) {
		vle = varg2list_end;
		varg2list_end = vl = varg2list;
		do {
			for(de1 = vl->L.d; de1->e; de1++)
				de1->dv.i = r[de1->dv.i];
			}
			while((vl = vl->next2) != vle);
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

 static derp *
derp_ogcopy(Static *S, ograd *og, derp *dnext, int k)
{
	last_d = dnext;
	if (og->varno < 0)
		og = og->next;
	while(og) {
		new_derp(S, og->varno, k, &og->coef);
		og = og->next;
		}
	return last_d;
	}

 static void
imap_alloc(Static *S)
{
	expr_v *v;
	int i, *r;
	linarg *tl;
	ASLTYPE *asl = S->asl;

	if (imap) {
		r = (int*)new_mblk(i = htcl(lasta*sizeof(int)));
		memcpy(r, imap, imap_len*sizeof(int));
		del_mblk(kimap, imap);
		imap = r;
		kimap = i;
		imap_len = (sizeof(char*)/sizeof(int)) << i;
		return;
		}
	i = amax1 > lasta ? amax1 : lasta;
	r = imap = (int*)new_mblk(kimap = htcl((i+100)*sizeof(int)));
	imap_len = (sizeof(char*)/sizeof(int)) << kimap;
	r += i = max_var1;
	while(i > 0)
		*--r = --i;
	i = nv0x;
	for(tl = asl->P.lalist; tl; tl = tl->lnext)
		if ((v = tl->v))
			r[v->a] = i++;
	r[noa] = i;
	}

 static void
comsubs(Static *S, int alen, cde *d)
{
	list *L;
	int a, i, j, jx, k, nzc1;
	int *r, *re;
	cexp *ce;
	derp *D, *dnext;
	dv_info *dvi;
	relo *R;
	expr *e;
	split_ce *cs;
	ograd *og;
	ASLTYPE *asl = S->asl;

	D = last_d;
	a = lasta00;
	dnext = 0;
	R = 0;
	nzc1 = nzc;
	nzc = 0;
	for(i = j = 0; i < nzc1; i++)
		if ((k = zci[i]) >= nv0x && k < max_var1)
			zci[j++] = k;
		else
			zc[k] = 0;
	if ((nzc1 = j)) {
		for(i = 0; i < nzc1; i++) {
			if ((j = zci[i] - nv0x) >= Ncom) {
				cs = asl->P.Split_ce + (j-Ncom);
				if ((r = cs->ce)) {
					re = r + *r;
					while(r++ < re)
						if (!zc[j = *r + nv0x]++)
							zci[nzc1++] = j;
					}
				}
			else if (j >= 0 && (L = cexps[j].cref)) {
				if (cexps[j].funneled) do {
					if (zc[j = L->item.i]
					 || asl->P.dv[j-nv0x].ll)
						continue;
					zc[j] = 1;
					zci[nzc1++] = j;
					}
					while((L = L->next));
				else do {
					if (!zc[L->item.i]++)
						zci[nzc1++] = L->item.i;
					}
					while((L = L->next));
				}
			}
		if (nzc1 > 1)
			zcsort(S, zc, zci, 0, nzc1, -1);
		R = new_relo1(S, dnext);
		i = 0;
		do {
			j = zci[i];
			jx = j - nv0x;
			zc[j] = 0;
			ce = &cexps[jx];
			e = ce->e;
			k = varp[jx]->a;
			if (ce->funneled) {
				if (j >= max_var)
					imap[((expr_vx*)varp[jx])->a0] = a;
				imap[k] = a++;
				}
			else {
				r = imap + ce->z.i;
				re = r + ce->zlen;
				while(r < re)
					*r++ = a++;
				/* zlen == 1 if e->op == OPNUM */
				imap[k] = e->op == OPNUM ? a-1 : imap[e->a];
				if (!ce->d
				 && jx < Ncom
				 && (og = (dvi = &asl->P.dv[jx])->ll)
				 && (dvi->nl || dvi->lt)) {
					dnext = R->D =
						derp_ogcopy(S,og,R->D,imap[j]);
					continue;
					}
				}
			dnext = R->D = derpcopy(S, ce, R->D);
			}
			while(++i < nzc1);
		}
	if (D || R) {
		if (!R)
			R = new_relo1(S, dnext);
		D = R->D = derpadjust(S, D, a, R->D);
		if (Intcast d->e->op != f_OPVARVAL)
			d->e->a = imap[d->e->a];
		}
	d->d = D;
	a += alen;
	d->zaplen = (a > lasta00 ? a - nv1 : 0)*sizeof(real);
	if (amax < a)
		amax = a;
	}

 static void
co_read(EdRead *R, cde *d, int k)
{
	d = (cde *)((char *)d + k*sizeof(cde));
	d->e = eread(R);
	}

 static void
co_walk(Static *S, cde *d)
{
	int alen;

	if (amax1 < lasta)
		amax1 = lasta;
	lasta = lasta0;
	last_d = 0;
	PSHV(last_e = 0;)
	ewalk(S, d->e, 1);
	PSHV(d->ee = last_e;)
	alen = lasta - lasta0;
	if (imap_len < lasta)
		imap_alloc(S);
	comsubs(S, alen, d);
	}

 static int
lpcompar(const void *a, const void *b, void *v)
{
	Not_Used(v);
	return ((linpart *)a)->v.i - ((linpart *)b)->v.i;
	}

 static linpart *
linpt_read(EdRead *R, int nlin)
{
	ASL *asl;
	int i0, needsort;
	int (*Xscanf)(EdRead*, const char*, ...);
	linpart *L, *rv;

	asl = R->asl;
	if (nlin <= 0)
		return 0;
	Xscanf = xscanf;
	L = rv = (linpart *)new_mblk(htcl(nlin*sizeof(linpart)));
	i0 = needsort = 0;
	do {
		if (Xscanf(R, "%d %lf", &L->v.i, &L->fac) != 2)
			badline(R);
		if (i0 > L->v.i)
			needsort++;
		i0 = L->v.i;
		L++;
		}
		while(--nlin > 0);
	if (needsort)
		qsortv(rv, L-rv, sizeof(linpart), lpcompar, NULL);
	return rv;
	}

 static int
funnelkind(Static *S, int i0, int *ip)
{
	cexp *ce, *cej;
	dv_info *dvi;
	int i, j, k, k1, k2, nzc0, rv;
	int *vr, *vre;
	linarg *la, **lap, **lape;
	list *tl;
	range *r;
	ASLTYPE *asl = S->asl;

	ce = cexps + i0;
	ce->vref = 0;
	if (!(nzc0 = nzc))
		return nocopy;
	k1 = k2 = nzcperm(S);
	dvi = asl->P.dv;
	if (i0 >= Ncom || !dvi[i0].nl)
		dvi = 0;
	for(i = k = 0; i < nzc; i++)
		if ((j = zci[i]) < nv0x) {
			if (k >= maxfwd) {
				k = 0;
				break;
				}
			vrefx[k++] = j;
			}
		else  {
			cej = cexps + (j -= nv0x);
			if (!(vr = cej->vref)) {
				if (dvi && j < Ncom && dvi[j].ll)
					for(tl = cej->cref; tl; tl = tl->next)
					    if (!zc[tl->item.i]++)
						zci[nzc++] = tl->item.i;
				continue;
				}
			if (!*++vr) {
				k = 0;
				break;
				}
			vre = vr + *vr;
			while(++vr <= vre)
				if (!zc[*vr]++)
					zci[nzc++] = *vr;
			}
	if (k2 < k)
		k2 = k;
	k2 += 2;
	if (k2 > nvref) {
		nvref = (maxfwd + 1)*(ncom_togo < vrefGulp
					? ncom_togo : vrefGulp);
		if (nvref < k2)
			nvref = 2*k2;
		vrefnext = (int *)M1alloc(nvref*Sizeof(int));
		}
	vr = ce->vref = vrefnext;
	vrefnext += k2;
	nvref -= k2;
	vr[0] = k1;
	vr[1] = k;
	vr += 2;
	i = 0;
	if (k)
		while(i < k)
			*vr++ = vrefx[i++];
	else
		while(i < k1)
			*vr++ = zci[i++];
	if (!nocopy) {
		k1 = 0;
		if (i0 < Ncom) {
			if ((lap = (asl->P.dv + i0)->nl))
				while((la = *lap++))
					if (la->nnz > 1)
						k1++;
			}
		else {
			r = asl->P.Split_ce[i0-Ncom].r;
			lap = r->lap;
			lape = lap + r->nintv;
			while(lap < lape)
				if ((*lap++)->nnz > 1)
					k1++;
			}
		if (k) {
			if (nderp > 3*(k+k1)) {
				*ip = k;
				return 2;
				}
			}
		}
	rv = 0;
	if (nocopy || nderp > 3*(nzc0+k1))
		rv = 1;
	while(nzc > nzc0)
		zc[zci[--nzc]] = 0;
	return rv;
	}

 static void
cexp_read(EdRead *R, int k, int nlin)
{
	cexp *ce;
	expr *e, *zero;
	Static *S = (Static *)R->S;
	ASLTYPE *asl = S->asl;

	ce = cexps + k - nv0x;
	ce->L = linpt_read(R, ce->nlin = nlin);
	e = eread(R);
	if (e->op == (efunc *)f_OPVARVAL) {
		/* avoid confusion with very simple defined variables */
		zero = (expr*)new_expr_n(S, 0.);
		zero->op = (efunc*)f_OPNUM;
		e = new_expr(S, 0, e, zero);
		}
	ce->e = e;
	}

 static cplist *
funnelderp(Static *S, int a, cplist *cl0)
{
	cplist *cl;
	ASL *asl = S->a;

	new_derp(S, a, lasta0, 0);
	cl = (cplist *)mem(sizeof(cplist));
	cl->next = cl0;
	cl->ca.i = imap[a];
	cl->cfa = last_d->c.rp = (real *)mem(sizeof(real));
	return cl;
	}

 static void
cexp_walk(Static *S, int k)
{
	ASLTYPE *asl = S->asl;
	cexp *ce;
	cplist *cl;
	dv_info *dvi;
	expr *e;
	funnel *f, **fp;
	int a, *ap, fk, i, j, ka, la0, na, nlin;
	linarg *la, **lap, **lape;
	linpart *L, *Le;
	ograd *og;
	range *r;

	/* for now, assume all common exprs need to be differentiated */

	ce = cexps + k;
	nlin = ce->nlin;
	L = ce->L;
	nocopy = 0;
	last_d = 0;
	la0 = lasta;
	nderps += nderp;
	nderp = 0;
	e = ce->e;
	switch(Intcast e->op) {
		case f_OPVARVAL:
		case f_OPNUM:
			ap = 0;
			break;
		default:
			ap = &e->a;
		}
	PSHV(last_e = 0;)
	ewalk(S, e, 1);
	PSHV(ce->ee = last_e;)
	ka = k + nv0x;
	if ((ce->zlen = lasta - la0))
		ce->z.i = la0;
	else {
		ce->z.i = k < Ncom ? ka : ((expr_vx*)varp[k-Ncom])->a0;
		ce->zlen = 1;
		}
	j = ap ? *ap : ka;
	if (nlin) {
		if (nlin == 1)
			new_derp(S, L->v.i, j, &L->fac);
		else if (k < Ncom) {
			dvi = asl->P.dv + k;
			if (dvi->lt)
				new_derp(S, dvi->lt->v->a, j, &dvi->scale);
			}
		for(Le = L + nlin; L < Le; L++)
			if (!zc[i = L->v.i]++)
				zci[nzc++] = i;
		}
	i = 0; /* only needed to shut up an erroneous warning */
	if ((fk = funnelkind(S, k, &i))) {
		/* arrange to funnel */
		fp = ka < nv0b ? &f_b : ka < nv0c ? &f_c : &f_o;
		ce->funneled = f = (funnel *)mem(sizeof(funnel));
		f->next = *fp;
		*fp = f;
		f->ce = ce;
		if (imap_len < lasta)
			imap_alloc(S);
		if (fk == 1) {
			f->fulld = last_d;
			a = lasta00;
			for(i = nzc; --i >= 0; )
				if ((j = zci[i]) >= nv0x && j < max_var1)
					imap[varp[j-nv0x]->a] = a++;
			if ((na = ce->zlen) || a > lasta00)
				na += a - nv1;
			f->fcde.zaplen = na*sizeof(real);
			i = nzc;
			derpadjust(S, last_d, a, 0);
			}
		else {
			f->fulld = 0;
			f->fcde.e = e;
			comsubs(S, ce->zlen, &f->fcde);
			memcpy(zci, vrefx, i*sizeof(int));
			}
		last_d = 0;
		cl = 0;
		while(i > 0) {
			if ((a = zci[--i]) >= nv0x && a < max_var1)
				a = varp[a-nv0x]->a;
			if (a != noa)
				cl = funnelderp(S, a, cl);
			}
		if (k < Ncom) {
			if ((lap = (asl->P.dv + k)->nl))
				while((la = *lap++))
					if (la->v)
						cl = funnelderp(S, la->v->a, cl);
			}
		else {
			r = asl->P.Split_ce[k-Ncom].r;
			lap = r->lap;
			lape = lap + r->nintv;
			while(lap < lape)
				if ((la = *lap++)->nnz > 1)
					cl = funnelderp(S, la->v->a, cl);
			}
		f->cl = cl;
		varp[k]->a = lasta0++;
		lasta = lasta0;
		}
	lasta0 = lasta;
	if (!(ce->d = last_d)
	 && e->op == OPNUM
	 && (og = asl->P.dv[k].ll)
	 && og->varno < 0)
		((expr_n*)e)->v = 0;
	while(nzc > 0)
		zc[zci[--nzc]] = 0;
	--ncom_togo;
	}

#ifdef PSHVREAD
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
co_adjust(ps_func *p, int n)
{
	ps_func *pe;
	psb_elem *b, *be;
	psg_elem *g, *ge;

	for(pe = p + n; p < pe; p++) {
		for(b = p->b, be = b + p->nb; b < be; b++)
			b->D.ef = hvadjust(b->D.ee);
		for(g = p->g, ge = g + p->ng; g < ge; g++) {
			g->esum.op = (efunc_n*)OPNUM;
			for(b = g->E, be = b + g->ns; b < be; b++)
				b->D.ef = hvadjust(b->D.ee);
			}
		}
	}
#else	/*PSHVREAD*/
#define co_adjust(a,b) /*nothing*/
#endif	/*PSHVREAD*/

 static void
ifadjust(ASLTYPE *asl, expr_if *e)
{
	real *Adjoints = asl->i.adjoints_;
	for(; e; e = e->next) {
		e->Tv.rp = &Adjoints[e->Tv.i];
		e->Fv.rp = &Adjoints[e->Fv.i];
#ifdef PSHVREAD
		e->Tf = hvadjust(e->Te);
		e->Ff = hvadjust(e->Fe);
#endif
		}
	}

 static void
vargadjust(ASLTYPE *asl, expr_va *e)
{
	de *d;
	real *Adjoints = asl->i.adjoints_;

	for(; e; e = e->next) {
		for(d = e->L.d; d->e; d++) {
			d->dv.rp = &Adjoints[d->dv.i];
#ifdef PSHVREAD
			d->ef = hvadjust(d->ee);
#endif
			}
		}
	}

 static void
funneladj1(ASLTYPE *asl, funnel *f)
{
	real	*a	= adjoints;
	derp	*d;
	cplist	*cl;
	funnel	*fnext;

	for(; f; f = fnext) {
		fnext = f->next;
		f->next = 0;
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
funneladjust(Static *S)
{
	cexp *c, *ce;
	linpart *L, *Le;
	ASLTYPE *asl = S->asl;

	c = cexps;
	for(ce = c + Ncom; c < ce; c++) {
		if ((L = c->L))
			for(Le = L + c->nlin; L < Le; L++)
				L->v.vp = (void*)&var_e[L->v.i];
#ifdef PSHVREAD
		c->ef = hvadjust(c->ee);
#endif
		}
#ifdef PSHVREAD
	for(ce += asl->P.ndvspout; c < ce; c++)
		c->ef = hvadjust(c->ee);
#endif

	funneladj1(asl, f_b);
	funneladj1(asl, f_c);
	funneladj1(asl, f_o);
	}

 static void
cg_zzap(ASLTYPE *asl)
{
	cgrad *cg, **cgp,**cgp1, **cgpe;

	cgp1 = Cgrad;
	cgpe = cgp1 + asl->i.n_con0;
	while(cgp1 < cgpe) {
		cgp = cgp1++;
		while((cg = *cgp))
			if (cg->coef)
				cgp = &cg->next;
			else
				*cgp = cg->next;
		}
	}

 static void
adjust_compl_rhs(ASLTYPE *asl)
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
		if (Cvar[i] && (e = C[i].e) && Intcast e->op == f_OPNUM
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
	derp *d, **dp;
	ASLTYPE *asl = S->asl;
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
	ifadjust(asl, iflist);
	vargadjust(asl, varglist);
	if (Ncom)
		funneladjust(S);
	co_adjust(asl->P.cps, asl->i.n_con0);
	co_adjust(asl->P.ops, n_obj);
	if (asl->i.n_con0 && !allJ)
		cg_zzap(asl);
	if (k_seen) {
		if (!A_vals)
			goff_comp_ASL((ASL*)asl);
		else if (Fortran)
			colstart_inc_ASL((ASL*)asl);
		}
	if (n_cc > nlcc && nlc < n_con
	 && !(flags & ASL_no_linear_cc_rhs_adjust))
		adjust_compl_rhs(asl);
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
	rvh->op = (efunc *)f_OPHOL;
	rvh->a = noa;
	return (expr *)rvh;
	}

 static expr *
bholread(EdRead *R)
{
	int i;
	expr_h *rvh;
	char *s;
	FILE *nl = R->nl;
	ASL *asl = R->asl;
	Static *S = (Static *)R->S;

	if (xscanf(R, "%d", &i) != 1)
		badline(R);
	rvh = (expr_h *)mem(memadj(sizeof(expr_h) + i));
	s = rvh->sym;
	if (fread(s, i, 1, nl) != 1)
		badline(R);
	s[i] = 0;
	rvh->op = (efunc *)f_OPHOL;
	rvh->a = noa;
	for(;;) switch(*s++) {
			case 0: goto break2; /* so we return at end of fcn */
			case '\n': R->Line++;
			}
 break2:
	return (expr *)rvh;
	}

 static int
qwalk(Static *S, expr *e)
{
	int i, j, k;
	expr **args, **argse;
	expr_f *rvf;
	argpair *ap, *ape;

 top:
	switch(optype[k = Intcast e->op]) {

		case 1:	/* unary */
			switch(k) {
			  case f_OPCPOW:
				if (qwalk(S, e->R.e))
					return 3;
				return 0;
			  case f_OPUMINUS:
				e = e->L.e;
				goto top;
			  case f_OP2POW:
				i = qwalk(S, e->L.e);
				if (i < 2)
					return i << 1;
			  }
			return 3;

		case 2:	/* binary */
			switch(k) {
			 case f_OPPLUS:
			 case f_OPMINUS:
				i = qwalk(S, e->L.e);
				if (i == 3)
					break;
				j = qwalk(S, e->R.e);
				if (i < j)
					i = j;
				return i;
			 case f_OPDIV:
				if (qwalk(S, e->R.e))
					return 3;
				e = e->L.e;
				goto top;
			 case f_OPMULT:
				i = qwalk(S, e->L.e);
				if (i < 3) {
					i += qwalk(S, e->R.e);
					if (i < 3)
						return i;
					}
				}
			return 3;

		case 6: /* sumlist */
			j = 0;
			args = e->L.ep;
			argse = e->R.ep;
			while(args < argse) {
				i = qwalk(S, *args++);
				if (j < i) {
					j = i;
					if (j == 3)
						break;
					}
				}
			return j;

		case 7: /* function call */
			rvf = (expr_f *)e;
			ap = rvf->ap;
			ape = rvf->ape;
			for(i = 0; ap < ape; ap++) {
				if (qwalk(S, ap->e))
					return 3;
				}

		case 9: /* OPNUM */
			return 0;
		case 10:/* OPVARVAL */
			k = (expr_v *)e - S->var_e;
			if (k < 0)
				goto k_adj;
			if (k < nv0x)
				return 1;
			if (k >= max_var) {
		k_adj:
				k = ((expr_vx*)e)->a1;
				return k < 0 ? 1 : S->asl->I.v_class[k-nv0x];
				}
			return S->asl->I.v_class[k - nv0x];
		}
	return 3;
	}

 static int
co_walkloop(Static *S, ps_func *f, int n, char *c, ograd **o)
{
	ps_func *fe;
	psb_elem *b, *be;
	psg_elem *g, *ge;
	int k, k1, kx;

	kx = 0;
	for(fe = f + n; f < fe; f++) {
		if (c) {
			k = *o++ != 0;
			for(g = f->g, ge = g + f->ng; g < ge; g++) {
				if (Intcast g->g->op != f_OP2POW) {
					k = 3;
					goto have_c;
					}
				if (g->nlin)
					k = 2;
				for(b = g->E, be = b + g->ns; b < be; b++) {
					if ((k1 = qwalk(S, b->D.e)) > 1) {
						k = 3;
						goto have_c;
						}
					k = 2;
					}
				}
			for(b = f->b, be = b + f->nb; b < be; b++) {
				if ((k1 = qwalk(S, b->D.e)) > k) {
					k = k1;
					if (k == 3)
						goto have_c;
					}
				}
 have_c:		*c++ = (char)k;
			if (kx < k)
				kx = k;
			}
		for(b = f->b, be = b + f->nb; b < be; b++)
			co_walk(S, &b->D);
		for(g = f->g, ge = g + f->ng; g < ge; g++) {
			ewalk(S, g->g, 1);
			for(b = g->E, be = b + g->ns; b < be; b++)
				co_walk(S, &b->D);
			}
		}
	return kx;
	}

 static void
do_ewalk(Static *S)
{
	int i, *map, n, nc;
	char *v, *ve;
	expr_v *e, *ee, **vp, **vpe;
	efunc *vv;
	linarg *la, **lap;
	ograd *og;
	cexp *c;
	ASLTYPE *asl = S->asl;

	psfind(S);
	nc = max_var1 - nv0x;
	noa = lasta00 + nc;
	nv1 = lasta00++;
	amax = lasta = lasta0 = max_var + nndv + 1;

	if (Ncom && asl->I.o_class) {
		i = max_var1 - nv0x;
		v = asl->I.v_class = (char*)M1alloc(i);
		ve = v + i;
		c = cexps;
		while(v < ve)
			*v++ = (char)((i = qwalk(S, (c++)->e)) ? i : 1);
		}

	n = Ncom + asl->P.ndvspout;
	while(nzc > 0)
		zc[zci[--nzc]] = 0;
	for(i = 0; i < n; i++) {
#ifdef DEBUG
		if (i == izork)
			printf("");
#endif
		cexp_walk(S, i);
		}

	if (imap_len < lasta)
		imap_alloc(S);
	f_b = funnelfix(f_b);
	f_c = funnelfix(f_c);
	f_o = funnelfix(f_o);

	asl->I.c_class_max = co_walkloop(S, asl->P.cps, asl->i.n_con0,
					asl->I.c_class, (ograd**)Cgrad);
	asl->I.o_class_max = co_walkloop(S, asl->P.ops, n_obj,
					asl->I.o_class, Ograd);
	del_mblk(kzc, zci);
	zci = zc = 0;	/* for debugging */
	e = var_e;
	vv = OPVARVAL;
	for(ee = e + nv0x + Ncom; e < ee; e++)
		e->op = vv;
	lap = &asl->P.lalist;
	map = imap;
	while((la = *lap))
		if ((e = la->v)) {
			e->op = vv;
			e->a = map[e->a];
			lap = &la->lnext;
			}
		else {
			og = la->nz;
			assert(og != 0);
			assert(og->next == 0);
			assert(og->varno < nv0x);
			la->v = var_e + og->varno;
			*lap = la->lnext;
			}
	if (n) {
		asl->P.vp = vp = varp;
#ifdef PSHVREAD
		asl->P.zlsave = zl;
#endif
		e = var_ex;
		for(ee = e + Ncom; e < ee; e++)
			*vp++ = e;
		vpe = vp + asl->P.ndvspout;
		while(vp < vpe)
			(*vp++)->op = vv;
		}
	}

 int
pfg_read_ASL(ASL *a, FILE *nl, int flags)
{
	ASLTYPE *asl;
	EdRead ER, *R;
	Jmp_buf JB;
	Static SS, *S;
	cgrad *cg, **cgp;
	char fname[128];
	expr_v *e;
	func_info *fi;
	int allG, i, i1, j, k, *ka, kseen, maxfwd1, nc, nc0, nco, nlin, nlogc, no;
	int nv, nvc, nvo, nvr, nvextra, nvx, readall;
	int (*Xscanf)(EdRead*, const char*, ...);
	ograd *og, **ogp;
	real t;
	size_t *kaz, nz;
	unsigned x;

	ASL_CHECK(a, asltype, who);
	flagsave_ASL(a, flags); /* includes allocation of LUv, LUrhs, A_vals or Cgrad, etc. */
	asl = (ASLTYPE*)a;
	S = S_init(&SS, asl);
	ed_reset(asl);
	SS.R = R = EdReadInit_ASL(&ER, a, nl, S);
	if (flags & ASL_return_read_err) {
		a->i.err_jmp_ = &JB;
		i = setjmp(JB.jb);
		if (i) {
			a->i.err_jmp_ = 0;
			return i;
			}
		}
	if ((nlogc = a->i.n_lcon_) && !(flags & ASL_allow_CLP)) {
		if (a->i.err_jmp_)
			return ASL_readerr_CLP;
		sorry_CLP(R, "logical constraints");
		}
	if (!(flags & ASL_find_default_no_groups))
		flags |= ASL_findgroups;
	Xscanf = xscanf;
	readall = flags & ASL_keep_all_suffixes;
	PSHV(asl->P.pshv_g1 = 1;)
	k_Elemtemp = htcl(sizeof(Elemtemp));
	allJ = (flags & ASL_J_zerodrop) == 0;
	allG = (flags & ASL_G_zerodrop) == 0;
	wantCgroups = flags & ASL_findCgroups;
	wantOgroups = flags & ASL_findOgroups;
	OPNUM = r_ops[f_OPNUM];
	OPVARVAL = r_ops[f_OPVARVAL];
	if (!size_expr_n)
		size_expr_n = sizeof(expr_n);
	S->size_exprn = size_expr_n;
	asl->P.rlist.next = asl->P.rlist.prev = (range*)&asl->P.rlist;
	if (nfunc)
		func_add(a);
	if (binary_nl)
		holread = bholread;
	else
		holread = aholread;

	asl->P.ncom = Ncom = comb + comc + como + comc1 + como1;
	nc0 = n_con;
	nc = nc0 + a->i.nsufext[ASL_Sufkind_con];
	no = n_obj;
	nvc = c_vars;
	nvo = o_vars;
	nco = nc + no + nlogc;
	if (no < 0 || nco <= 0)
		scream(R, ASL_readerr_corrupt,
			"pshvread: nc = %d, no = %d, nlogc = %d\n",
			nc0, no, nlogc);
	if (pi0) {
		memset(pi0, 0, nc*sizeof(real));
		if (havepi0)
			memset(havepi0, 0, nc);
		}
	nvextra = a->i.nsufext[ASL_Sufkind_var];
	nv0r = nvr = n_var;
	asl->P.nv0_ = lasta00 = nv0x = nvx = n_var + nvextra;
	max_var1 = max_var = nv = nvr + Ncom + nvextra;
	combc = comb + comc;
	ncom0 = ncom_togo = combc + como;
	nzclim = max_var >> 3;
	ncom1 = comc1 + como1;
	nv0b = nvr + comb;
	nv0c = nv0b + comc;
	if ((maxfwd1 = maxfwd + 1) > 1)
		nvref = maxfwd1*((ncom0 < vrefGulp ? ncom0 : vrefGulp) + 1);
	x = nco*sizeof(cde) + no*sizeof(ograd *)
		+ nv*sizeof(expr_v)
		+ nfunc*sizeof(func_info *)
		+ nvref*sizeof(int)
		+ no;
	SS.nvar0 = a->i.n_var0;
	if (!(SS.nvinc = a->i.n_var_ - SS.nvar0 + nvextra))
		SS.nvar0 += ncom0 + ncom1;
	if (flags & ASL_find_co_class)
		x += nco;
	if (X0)
		memset(X0, 0, nvx*sizeof(real));
	if (havex0)
		memset(havex0, 0, nvx);
	e = var_e = (expr_v *)M1zapalloc(x);
	con_de = (cde *)(e + nv);
	lcon_de = con_de + nc;
	obj_de = lcon_de + nlogc;
	Ograd = (ograd **)(obj_de + no);
	var_ex = e + nvx;
	var_ex1 = var_ex + ncom0;
	for(k = 0; k < nv; e++) {
		e->op = (efunc *)f_OPVARVAL;
		e->a = k++;
		}
	funcs = (func_info **)(Ograd + no);
	zci = 0;
	zc_upgrade(S); /* initially allow nv+1 entries in zci and zc[-1] */
	vrefx = (int *)(funcs + nfunc);
	if (Ncom) {
		cexps = 0;
		asl->P.ndvspout = 0;
		memset(asl->P.dv = (dv_info*)mem(Ncom*sizeof(dv_info)),
			0, Ncom*sizeof(dv_info));
		cexp_upgrade(S, Ncom);
		for(k = 0; k < Ncom; k++)
			varp[k] = &var_ex[k];
		}
	objtype = (char *)(vrefx + nvref);
	if (flags & ASL_find_co_class) {
		asl->I.o_class = objtype + no;
		asl->I.c_class = asl->I.o_class + no;
		}
	if (nvref) {
		vrefnext = vrefx + maxfwd1;
		nvref -= maxfwd1;
		}
	if (n_cc && !cvar)
		cvar = (int*)M1alloc(nc*sizeof(int));
	if (cvar)
		memset(cvar, 0, nc*sizeof(int));
	last_d = 0;
	ka = 0;
	kaz = 0;
	kseen = 0;
	nz = 0;
	for(;;) {
		ER.can_end = 1;
		i = edag_peek(R);
		if (i == EOF) {
			fclose(nl);
			do_ewalk(S);
			if (imap)
				del_mblk(kimap, imap);
			/* Make amax long enough for nlc to handle */
			/* var_e[i].a for common variables i. */
			if (ncom0) {
				i = comb + como;
				if (i < combc)
					i = combc;
				if ((i += noa + 1) > amax)
					amax = i;
				}
			adjoints = (real *)M1zapalloc(amax*Sizeof(real));
			adjoints_nv1 = &adjoints[nv1];
			nderps += nderp;
			adjust(S, flags);
			nzjac = nz;
			if (!Lastx)
				Lastx = (real *)M1alloc(nvx*sizeof(real));
#ifdef PSHVREAD
			a->i.ncxval = (int*)M1zapalloc(nco*sizeof(int));
			a->i.noxval = a->i.ncxval + nc;
			a->p.Conival = a->p.Conival_nomap = conpival_ASL;
			a->p.Congrd  = a->p.Congrd_nomap  = conpgrd_ASL;
			a->p.Objval = a->p.Objval_nomap = objpval_ASL;
			a->p.Objgrd = a->p.Objgrd_nomap = objpgrd_ASL;
			a->p.Conval = conpval_ASL;
			a->p.Jacval = jacpval_ASL;
			a->p.Lconval= lconpval_ASL;
			a->p.Hvcomp = a->p.Hvcomp_nomap = hvpcomp_ASL;
			a->p.Hvcompd = hvpcompd_ASL;
			a->p.Hvcomps = hvpcomps_ASL;
			a->p.Hvinit = a->p.Hvinit_nomap = hvpinit_ASL;
			a->p.Hesset = hes_setup;
			a->p.Xknown = xp2known_ASL;
			a->p.Duthes = a->p.Duthes_nomap = duthes_ASL;
			a->p.Fulhes = a->p.Fulhes_nomap = fullhes_ASL;
			a->p.Sphes  = a->p.Sphes_nomap  = sphes_ASL;
			a->p.Sphset = a->p.Sphset_nomap = sphes_setup_ASL;
#else
			a->p.Xknown = xp1known_ASL;
#endif
			return prob_adj_ASL(a);
			}
		ER.can_end = 0;
		k = -1;
		switch(i) {
			case 'C':
				Xscanf(R, "%d", &k);
				if (k < 0 || k >= nc0)
					badline(R);
				co_read(R, con_de, k);
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
			case 'L':
				Xscanf(R, "%d", &k);
				if (k < 0 || k >= nlogc)
					badline(R);
				co_read(R, lcon_de, k);
				ewalk(S, lcon_de[k].e, 0);
				break;
			case 'V':
				if (Xscanf(R, "%d %d %d", &k, &nlin, &j) != 3)
					badline(R);
				if (k >= SS.nvar0)
					k += SS.nvinc;
				if (k < nvr || k >= nv)
					badline(R);
				cexp_read(R, k, nlin);
				break;
			case 'G':
				if (Xscanf(R, "%d %d", &j, &k) != 2
				|| j < 0 || j >= no || k <= 0 || k > nvo)
					badline(R);
				ogp = Ograd + j;
				while(k--) {
					if (Xscanf(R, "%d %lf", &i, &t) != 2)
						badline(R);
					if (allG || t) {
						*ogp = og = (ograd *)
							mem(sizeof(ograd));
						ogp = &og->next;
						og->varno = i;
						og->coef = t;
						}
					}
				*ogp = 0;
				break;
			case 'J':
				if (Xscanf(R, "%d %d", &j, &k) != 2
				|| j < 0 || j >= nc || k <= 0 || k > nvc)
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
				j = 0;
				while(k--) {
					if (kseen) {
						if (Xscanf(R, "%d %lf", &i, &t) != 2)
							badline(R);
						}
					else
						if (Xscanf(R, "%d %d %lf", &i, &j, &t) != 3)
							badline(R);
					*cgp = cg = (cgrad *)
						mem(sizeof(cgrad));
					cgp = &cg->next;
					cg->varno = i;
					cg->goff = j;
					cg->coef = t;
					}
				*cgp = 0;
				break;
			case 'O':
				if (Xscanf(R, "%d %d", &k, &j) != 2
				 || k < 0 || k >= no)
					badline(R);
				objtype[k] = j;
				co_read(R, obj_de, k);
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
					x = nvx*sizeof(real);
					if (want_xpi0 & 4)
						x += nvx;
					X0 = (real *)M1zapalloc(x);
					if (want_xpi0 & 4)
						havex0 = (char*)(X0 + nvx);
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
#undef asl

#ifdef PSHVREAD

 static int
heswork(expr *e)
{
	argpair *da, *dae;
	int i, j, n;
	expr *e1, **ep;
	de *d;

	n = 0;
	if (e && e->op != OPVARVAL &&  e->op != OPNUM) do {
	    switch(e->a) {

		case Hv_timesR:
		case Hv_timesL:
		case Hv_negate:
		case Hv_plterm:
			n += 4;
			break;

		case Hv_binaryR:
		case Hv_unary:
			n += 6;
			break;

		case Hv_timesLR:
			n += 10;
			break;

		case Hv_binaryLR:
			n += 14;
			break;


		case Hv_vararg:
			d = ((expr_va*)e)->L.d;
			i = heswork(d->e);
			while((e1 = (++d)->e)) {
				j = heswork(e1);
				if (i < j)
					i = j;
				}
			n = i + 2;
			break;

		case Hv_if:
			i = heswork(((expr_if *)e)->Tf);
			j = heswork(((expr_if *)e)->Ff);
			if (i < j)
				i = j;
			n += i + 2;
			break;


		case Hv_sumlist:
			for(ep = e->R.ep; (e1 = *ep); ep++)
				n++;
			break;

		case Hv_func:
			da = ((expr_f *)e)->da;
			dae = ((expr_f *)e)->dae;
			i = dae - da;
			n += i*(i+4);
			break;

		case Hv_plusR:
		case Hv_plusL:
		case Hv_minusR:
			n += 2;
			break;

		case Hv_plusLR:
		case Hv_minusLR:
			n += 3;
			break;

		default:/*DEBUG*/
			fprintf(Stderr, "bad e->a = %d in heswork\n", e->a);
			exit(1);
		}
	    } while((e = e->fwd));
	return n;
	}

 static cexp *
hesfunnel(Static *S, int *ip, int nrefs, ograd **ogp, linarg ***lapp, linarg ***lapep)
{
	cexp *c;
	expr *e;
	int hw, i, k, m, n;
	linarg *la, **lap, **lape;
	ograd *og;
	dv_info *dvi;
	list *L;
	range *r;
	ASLTYPE *asl = S->asl;

	i = *ip;
	c = cexps + i;
	n = 0;
	if (i >= Ncom) {
		r = asl->P.Split_ce[i-Ncom].r;
		lap = r->lap;
		lape = lap + (n = r->n);
		}
	else if (!(lap = (dvi = asl->P.dv + i)->nl)) {
		asl->P.linmultr++;
		og = dvi->ll;
		if (og->varno < 0)
			og = og->next;
		for(*ogp = og; og; og = og->next)
			n++;
		if (n > 1 && asl->p.hffactor > 0) {
			asl->P.linhesfun++;
			*ip = n;
			return c;
			}
		return 0;
		}
	else {
		lape = lap;
		while(*lape)
			lape++;
		n = lape - lap;
		}
	if (!(e = c->ef))
		return 0;
	*lapp = lap;
	*lapep = lape;
	*ogp = 0;
	asl->P.nlmultr++;
	while(lap < lape) {
		la = *lap++;
		for(og = la->nz; og; og = og->next)
			if (!zc[og->varno]++)
				zci[nzc++] = og->varno;
		}
	m = nzc;
	while(nzc > 0)
		zc[zci[--nzc]] = 0;
	for(L = c->cref; L; L = L->next)
		n++;
	if (m > n)
		m = n;
	if ((k = m*nrefs - n) <= 0)
		return 0;
	hw = heswork(e);
	if (k*hw <= nrefs*n*(n+3)*asl->p.hffactor)
		return 0;
	*ip = n;
	asl->P.nlhesfun++;
	return c;
	}

 ASL_pfgh *
pscheck_ASL(ASL *a, const char *who1)
{
	ASL_CHECK(a, ASL_read_pfgh, who1);
	return (ASL_pfgh*)a;
	}

 static void
hes_setup(ASL *a, int flags, int obj, int nobj, int con, int ncon)
{
	range *r, *r0;
	int *c, *ce, *cei, i, k, n, n0, nmax, nvx, *z;
	hes_fun *hf, *hfth;
	cexp *C;
	linarg *la, **lap, **lape;
	ograd *og;
	list *L;
	expr_v **vp;
	psb_elem *b;
	ASL_pfgh *asl;
	Static SS, *S;

	asl = pscheck_ASL(a, "hes_setup");

	S = S_init(&SS, asl);
	Ncom = asl->P.ncom;
	nv0x = nvx = asl->P.nv0_;
	max_var = nvx + Ncom;
	max_var1 = max_var + asl->P.ndvspout;
	nzclim = max_var1 >> 3;
	varp = asl->P.vp;
	zl = asl->P.zlsave;

	zc_upgrade(S);	/* make zc and zci available */
	n = nmax = 0;
	r0 = (range*)&asl->P.rlist;
	for(r = asl->P.rlist.next; r != r0; r = r->rlist.next) {
		if (r->n <= 0) {
			r->cei = 0;
			continue;
			}
		if (nmax < r->n)
			nmax = r->n;
		n0 = 0;
		cei = 0;
		for(b = r->refs; b; b = b->next) {
			i = b->conno;
			if (i < 0) {
				i = -i - 2;
				if (i < obj || i-obj >= nobj)
					continue;
				}
			else {
				if (i < con || i-con >= ncon)
					continue;
				}
			if ((c = b->ce)) {
				if (!n0) {
					cei = c;
					n0 = *c;
					}
				ce = c + *c;
				while(c < ce)
					if (!zc[i = *++c]++)
						zci[n++] = i;
				}
			}
		if (n0 < n)
			cei = 0;
		if ((r->cei = cei)) {
			while(n > 0)
				zc[zci[--n]] = 0;
			continue;
			}
		if (!n)
			continue;
		r->cei = cei = (int*)mem((n+1)*sizeof(int));
		*cei++ = n;
		if (n > 1)
			qsortv(zci, n, sizeof(int), hscompar, S);
		cei += n;
		do zc[*--cei = zci[--n]] = 0;
			while(n > 0);
		}
	z = zc + nvx;
	k = -1;
	for(r = asl->P.rlist.next; r != r0; r = r->rlist.next) {
		if (!(cei = r->cei))
			continue;
		c = cei + 1;
		ce = c + *cei;
		do {
			if (z[i = *c++]++ && k < i)
				k = i;
			} while(c < ce);
		}
	hfth = 0;
	asl->P.linmultr = asl->P.linhesfun = 0;
	asl->P.nlmultr = asl->P.nlhesfun = 0;
	lap = lape = 0; /* silence false "used uninitialized" warning */
	while (k >= 0)
		if (z[i = k--] > 1
		 && (C = hesfunnel(S,&i,z[i],&og,&lap,&lape))) {
			hf = (hes_fun*)mem(sizeof(hes_fun));
			hf->hfthread = hfth;
			C->hfun = hfth = hf;
			hf->c = C;
			hf->n = i;
			if ((hf->og = og)) {
				hf->vp = 0;
				hf->grdhes = 0;
				}
			else {
				hf->vp = vp = (expr_v **)
						mem(i*sizeof(expr_v*));
				while(lap < lape) {
					la = *lap++;
					*vp++ = la->v;
					}
				for(L = C->cref; L; L = L->next)
					*vp++ = varp[L->item.i - nvx];
				i += i*i;
				hf->grdhes = (real*)mem(i*sizeof(real));
				}
			}
	if ((asl->I.hesthread = hfth)) {
		asl->I.x0kind_init = ASL_need_funnel;
		if (x0kind != ASL_first_x)
			x0kind |= ASL_need_funnel;
		}
	else
		asl->I.x0kind_init = 0;
	del_mblk(kzc, zci);
	asl->P.hes_setup_called = 1 | (flags << 2);
	asl->P.khesoprod = 5;
	asl->P.hop_free = 0;
	asl->P.nmax = nmax;
	if (!(flags & 1))
		return;
	i = htcl(n = nvx*sizeof(range*));
	k = htcl(3*n);
	if (k == i + 1) {
		asl->P.rtodo = (range**)new_mblkzap(asl, k);
		asl->P.utodo = (uHeswork **)(asl->P.rtodo + nvx);
		asl->P.otodo = (Hesoprod**)(asl->P.utodo + nvx);
		}
	else {
		asl->P.rtodo = (range**)new_mblkzap(asl, i);
		asl->P.utodo = (uHeswork**)new_mblkzap(asl, i);
		asl->P.otodo = (Hesoprod**)new_mblkzap(asl, i);
		}
	k = htcl(nmax*(sizeof(int)+sizeof(real)));
	asl->P.dOscratch = (real *)new_mblkzap(asl, k);
	asl->P.iOscratch = (int *)(asl->P.dOscratch + nmax);
	}

#endif	/* PSHVREAD */

#ifdef __cplusplus
}
#endif
