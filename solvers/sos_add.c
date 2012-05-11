/****************************************************************
Copyright (C) 2002 AMPL Optimization LLC
All Rights Reserved
Based largely on suf_sos.c, which bears the following Copyright
notice and disclaimer.  AMPL Optimization LLC similarly disclaims
all warranties with regard to this software and grants permission
to use, copy, modify, and distribute it and its documentation. */

/****************************************************************
Copyright (C) 1999-2001 Lucent Technologies
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

/* This is a replacement for suf_sos() for solvers that know about
   SOS1 sets but do not supply the requisite convexity constraints
   (or do not know about SOS2 sets), or that know nothing of SOS
   sets but do handle integer variables (in which case sos_finish()
   should be called with null values for the arguments in which
   SOS information is returned, such as sosref_p).

   Before calling the .nl reader, one invokes

	void *SI = sos_add(nl,flags);

   (where nl is the open FILE* that will be passed to the .nl reader)
   to prepare for sos_finish to insert convexity constraints and adjust
   n_con, n_var, nzc, and nbv, etc.  After calling the .nl reader,
   one invokes

	if (SI)
		sos_finish(SI,...);

   to finish modifying the problem -- adding the convexity constraints,
   etc., for SOS sets implied by suffixes .sosno and .ref, and retaining
   SOS1 and SOS2 sets added when AMPL linearizes nonconvex piecewise-
   linear terms.
*/

#include "asl.h"
#define SKIP_NL2_DEFINES
#include "nlp.h"
#include "nlp2.h"
#include "asl_pfg.h"
#include "asl_pfgh.h"
#undef cde

#ifdef __cplusplus
extern "C" {
static int  compar(const void*, const void*, void*);
static int rcompar(const void*, const void*, void*);
}
#endif

#undef nzc

 typedef struct
SOSinfo {
	SufDesc mysd[2];
	Char **zap[5];
	ASL *asl;
	SufDesc *grefd, *refd;
	int *col1, *g0, **gp0, **gp1, **gpe, *v0, *ve, *z;
	int nnc, nnv, nnz, nsos, nsos1, nsos2, nsosnz, nsosnz1, nzc;
	}
	SOSinfo;

 typedef struct
Coninfo {
	cgrad **cgp, *cg;
	int *z;
	real *clu, *cu, *lu, *u;
	int m, nbin;
	}
	Coninfo;

 static real pl_bigM, mpl_bigM;

 static SufDesc*
refd_copy(SufDesc *nu, SufDesc *ol, int n)
{
	if (!ol)
		return ol;
	memcpy(nu, ol, sizeof(SufDesc));
	memcpy(nu->u.r = (real*)Malloc(n*sizeof(real)), ol->u.r, n*sizeof(real));
	return nu;
	}

 static int
compar(const void *a, const void *b, void *v)
{
	int k;
	Not_Used(v);
	if ((k = *(int*)a - *(int*)b))
		return k;
	return ((int*)a)[1] - ((int*)b)[1];
	}

 static int
rcompar(const void *a, const void *b, void *v)
{
	real *r = (real *)v, t;
	t = r[*(int*)a] - r[*(int*)b];
	if (t == 0.)
		return 0;
	return t < 0. ? -1 : 1;
	}

 static int
refcomp(const void *a, const void *b, void *v)
{
	real *r = (real *)v, t;
	t = r[((int*)a)[1]] - r[((int*)b)[1]];
	if (t == 0.)
		return 0;
	return t < 0. ? -1 : 1;
	}

 static void
reorder(int *ind, real *ref, int j0, int k, int *p)
{
	int i, j, n, ti;
	real tr;

	ref += j0;
	ind += j0;
	n = k - j0;
	for(i = 0; i < n; i++)
		p[i] = i;
	qsortv(p, n, sizeof(int), rcompar, ref);
	for(i = 0; i < n; i++) {
		if ((j = p[i]) > i) {
			ti = ind[i];
			tr = ref[i];
			for(k = i;;) {
				ind[k] = ind[j];
				ref[k] = ref[j];
				j = p[k = j];
				p[k] = -1 - j;
				if (j == i) {
					ind[k] = ti;
					ref[k] = tr;
					break;
					}
				}
			}
		}
	}

 static void
turnon(int *z, int n)
{
	while(n-- > 0)
		*z++ = 1;
	}

 static int
nonbinary(int k, Coninfo *CI, real *LU)
{
	int *z;
	real *lu, *u;

	if ((u = CI->u)) {
		lu = CI->lu + k;
		u += k;
		}
	else {
		lu = CI->lu + 2*k;
		u = lu + 1;
		}
	if ((LU[1] = *u) > pl_bigM)
		LU[1] = pl_bigM;
	if ((LU[0] = *lu) < mpl_bigM) {
		LU[0] = mpl_bigM;
		return 1;
		}
	if (LU[0] != 0.)
		return 1;
	if (!(z = CI->z) || !z[k] || *u != 1.)
		return 1;
	CI->nbin++;
	return 0;
	}

/* Call sos_add() prior to the .nl reader:  sos_add scans incoming
   suffixes for sosno, ref, sos and sosref and if found, increments
   n_var and n_con in preparation for sos_finish(), which must be
   called after the .nl reader to add convexity constraints, etc.
   if sos_add has provided a nonzero return.

   We use a private suftab.  Users can omit the items in our suftab
   from theirs, but should provide "priority" if desired.
 */

 static SufDecl
suftab[] = {
	{ (char*)"ref", 0, ASL_Sufkind_var | ASL_Sufkind_real },
	{ (char*)"sos", 0, ASL_Sufkind_var },
	{ (char*)"sosno", 0, ASL_Sufkind_var | ASL_Sufkind_real },
	{ (char*)"sosref", 0, ASL_Sufkind_var | ASL_Sufkind_real }
	};

 Char*
sos_add_ASL(ASL *asl, FILE *f, int flags)
{
	Char **M1state1, **M1state2;
	EdRead ER, *R;
	SOSinfo *SI;
	SufDesc *gd, *grefd, *refd, *sds[4], *vd;
	int i, j, k, n, ng, nnc, nnv, nnz;
	int nsos, nsos1, nsos2, nsosnz, nsosnz1, tg;
	int *col1, *g0, *g, *g1, *ge, **gp, **gp0, **gp1, **gpe;
	int nss[5], *v, *v0, *ve, *z, *zg;
	long ft0;
	real *gn;
	size_t L;

	SI = 0;
	M1state1 = asl->i.Mbnext;
	M1state2 = asl->i.Mblast;

	/* Save and restore asl->i.suff stuff in case caller */
	/* has already called suf_declare(). */

	nss[0] = asl->i.nsuffixes;
	memcpy(nss+1, asl->i.nsuff, 4*sizeof(int));
	memcpy(sds, asl->i.suffixes, 4*sizeof(SufDesc*));
	suf_declare(suftab, sizeof(suftab)/sizeof(SufDecl));

	R = EdReadInit_ASL(&ER, asl, f, 0);
	if (binary_nl)
		xscanf = bscanf;
	else
		xscanf = ascanf;
	ft0 = ftell(f);
	while((i = edag_peek(R)) == 'S')
		Suf_read_ASL(R, 0);
	fseek(f, ft0, SEEK_SET);

	nsos1 = nsos2 = 0;
	gd = refd = grefd = 0;
	z = 0;
	if (!(flags & ASL_suf_sos_ignore_sosno)
	 && (gd = suf_get("sosno", ASL_Sufkind_var | ASL_Sufkind_input))
	 && (grefd = suf_get("ref", ASL_Sufkind_var | ASL_Sufkind_input)))
		nsos1 = 1;
	if (!(flags & ASL_suf_sos_ignore_amplsos)
	 && (vd = suf_get("sos", ASL_Sufkind_var | ASL_Sufkind_input))
	 && (refd = suf_get("sosref", ASL_Sufkind_var|ASL_Sufkind_input)))
		nsos2 = 1;
	else if (!nsos1)
		goto ret;
	nnc = nnv = nnz = nsos = nsosnz = nsosnz1 = 0;
	n = n_var;
	col1 = v0 = ve = 0;
	if (nsos2) {
		v = v0 = vd->u.i;
		ve = v + n;
		for(j = 0; j < n && v[j] <= 0; j++);
		col1 = v + j;
		for(; j < n; j++) {
			if ((i = v[j]) & 2 && !(i&1)) {
				nsosnz++;
				if (nsos < (i >>= 4))
					nsos = i;
				}
			}
		}
	gp0 = gp1 = gpe = 0;
	g0 = 0;
	if (nsos1) {
		nsos1 = 0;
		gn = gd->u.r;
		for(i = ng = 0; i < n; i++) {
			if ((j = gn[i]))
				ng++;
			}
		if (!ng) {
			/* should not happen */
			if (nsos2)
				goto havenz;
			goto ret;
			}
		j = (niv | nbv | nlvbi | nlvci | nlvoi) ? n : 0;
		L = (ng+1)*sizeof(int*) + (2*ng + j)*sizeof(int);
		gp = gp0 = (int**)Malloc(L);
		zg = g0 = (int*)(gpe = gp + ng + 1);
		for(i = 0; i < n; i++)
			if ((k = gn[i])) {
				*zg++ = k;
				*zg++ = i;
				}
		ge = zg;
		qsortv(g0, ng, 2*sizeof(int), compar, 0);
		z = 0;
		if (j) {
			memset(z = zg, 0, n*sizeof(int));
			if (nlvbi)
				turnon(z + (nlvb - nlvbi), nlvbi);
			if (nlvci)
				turnon(z + (nlvc - nlvci), nlvci);
			if (nlvoi)
				turnon(z + (nlvo - nlvoi), nlvoi);
			if ((j = niv + nbv))
				turnon(z + (n-j), j);
			}
		*gp++ = g = g1 = zg = g0;
		j = -((tg = *g) < 0);
		for(;;) {
			if ((g += 2) >= ge || *g != tg) {
				/* Ignore SOS1 sets of 1 element   */
				/* and SOS2 sets of <= 2 elements. */
				if ((j += ((g-zg)>>1)) >= 2) {
					nsosnz1 += j;
					nsos1++;
					if (!gp1 && tg > 0)
						gp1 = gp - 1;
					if (g1 == zg)
						g1 = g;
					else while(zg < g)
						*g1++ = *zg++;
					*gp++ = g1;
					}
				if (g >= ge)
					break;
				j = -((tg = *(zg = g)) < 0);
				}
			}
		nsosnz += nsosnz1;
		ge = g1;
		gpe = gp - 1;
		nsos += nsos1;
		}
 havenz:
	if (!nsos) {
		if (gp0)
			free(gp0);
		goto ret;
		}
	if (nsos1) {
		if (!gp1)
			gp1 = gpe;
		/* scan SOS2s */
		/* Here we must assume the worst case for nonbinary
		 * variables (since their bounds are not yet available):
		 * both lower and upper bounds nonzero.  In sos_finish(),
		 * we may discover that we need fewer new constraints
		 * if some bounds turn out to be zero.
		 */
		for(gp = gp0; gp < gp1; gp++) {
			g = gp[0];
			ge = gp[1] - 2;
			j = (ge - g) >> 1;
			nnv += j - 2;
			nnc += j;
			nnz += 4*j - 3;
			for(i = 0; i < 2; i++) {
				k = i ? ge[1] : g[1];
				if (!z || !z[k]) {
					nnv++;
					nnc += 2;
					nnz += 4;
					}
				}
			nnc += k = (ge - g) - 2;
			nnz += 3*k;
			}
		/* scan SOS1s */
		for(; gp < gpe; gp++) {
			g = gp[0];
			ge = gp[1];
			j = (ge - g) >> 1;
			nnc++;
			nnz += j;
			nnc += k = ge - g;
			nnv += k >> 1;
			nnz += k << 1;
			}
		}
	SI = (SOSinfo*)Malloc(sizeof(SOSinfo));
	memset(SI->zap, 0, sizeof(SI->zap));
	SI->asl = asl;
	if (v0) {
		L = n*sizeof(int);
		memcpy(SI->v0 = (int*)Malloc(L), v0, L);
		col1 = SI->v0 + (col1 - v0);
		ve = (v0 = SI->v0) + n;
		}
	SI->col1 = col1;
	SI->g0 = g0;
	SI->gp0 = gp0;
	SI->gp1 = gp1;
	SI->gpe = gpe;
	SI->grefd = grefd = refd_copy(&SI->mysd[0], grefd, n);
	asl->i.nsufext[ASL_Sufkind_con]  += SI->nnc = nnc;
	asl->i.nsufext[ASL_Sufkind_var]  += SI->nnv = nnv;
	asl->i.nsufext[ASL_Sufkind_prob] += SI->nnz = nnz;
	SI->nsos = nsos;
	SI->nsos1 = nsos1;
	SI->nsos2 = nsos2;
	SI->nsosnz = nsosnz;
	SI->nsosnz1 = nsosnz1;
	SI->refd = refd = refd_copy(&SI->mysd[1], refd, n);
	SI->v0 = v0;
	SI->ve = ve;
	SI->z = z;
	SI->nzc = asl->i.nzc_;
	asl->i.nzc_ += nnz;
 ret:
	M1free_ASL(&asl->i, M1state1, M1state2);
	if (SI) {
		SI->zap[0] = M1record(SI);
		if (gp0)
			SI->zap[1] = M1record(gp0);
		if (v0)
			SI->zap[2] = M1record(v0);
		if (grefd)
			SI->zap[3] = M1record(grefd->u.r);
		if (refd)
			SI->zap[4] = M1record(refd->u.r);
		}
	asl->i.nsuffixes = nss[0];
	memcpy(asl->i.nsuff, nss+1, 4*sizeof(int));
	memcpy(asl->i.suffixes, sds, 4*sizeof(SufDesc*));

	return (Char*)SI;
	}

 static real LUge[2];

 static cgrad**
newcon(Coninfo *CI, int ge)
{
	int m = CI->m++;
	real *lu;
	static real LU1[2] = { 0., 1. };

	lu = ge ? LUge : LU1;
	if (CI->cu) {
		CI->clu[m] = lu[0];
		CI-> cu[m] = lu[1];
		}
	else {
		m <<= 1;
		CI->clu[m]   = lu[0];
		CI->clu[m+1] = lu[1];
		}
	return CI->cgp++;
	}

 static void
newcoef(Coninfo *CI, cgrad ***p, int k, real t)
{
	cgrad *cg = CI->cg++;
	**p = cg;
	*p = &cg->next;
	cg->varno = k;
	cg->coef = t;
	}

 static void
Bound(Coninfo *CI, int j, int k, real *LU)
{
	cgrad **cgb;
	if (LU[1]) {
		cgb = newcon(CI, 1);
		if (j < k) {
			newcoef(CI, &cgb, j, -1.);
			newcoef(CI, &cgb, k, LU[1]);
			}
		else {
			newcoef(CI, &cgb, k, LU[1]);
			newcoef(CI, &cgb, j, -1.);
			}
		*cgb = 0;
		}
	if (LU[0]) {
		cgb = newcon(CI, 1);
		if (j < k) {
			newcoef(CI, &cgb, j, 1.);
			newcoef(CI, &cgb, k, -LU[0]);
			}
		else {
			newcoef(CI, &cgb, k, -LU[0]);
			newcoef(CI, &cgb, j, 1.);
			}
		*cgb = 0;
		}
	}

 static void
Bound2(Coninfo *CI, int j, int k0, int k, real *LU)
{
	cgrad **cgb;
	if (LU[1]) {
		cgb = newcon(CI, 1);
		if (j < k0) {
			newcoef(CI, &cgb, j, -1.);
			newcoef(CI, &cgb, k0, LU[1]);
			}
		else {
			newcoef(CI, &cgb, k0, LU[1]);
			newcoef(CI, &cgb, j, -1.);
			}
		newcoef(CI, &cgb, k,  LU[1]);
		*cgb = 0;
		}
	if (LU[0]) {
		cgb = newcon(CI, 1);
		if (j < k0) {
			newcoef(CI, &cgb, j, 1.);
			newcoef(CI, &cgb, k0, -LU[0]);
			}
		else {
			newcoef(CI, &cgb, k0, -LU[0]);
			newcoef(CI, &cgb, j, 1.);
			}
		newcoef(CI, &cgb, k,  -LU[0]);
		*cgb = 0;
		}
	}

 static void
debugchk(const char *what, int expected, int got, int exact)
{
	if (got != expected && (exact || got > expected)) {
		fprintf(Stderr, "sos_finish: expected %s = %d, got %d\n",
			what, expected, got);
		mainexit_ASL(2);
		}
	}

 int
sos_finish_ASL(ASL *asl, void **VP, int flags, int *nsosnz_p, int **sospri_p,
	int *copri, int **sosbeg_p, int **sosind_p, real **sosref_p)
{
	Char **vp;
	Coninfo CI;
	SOSinfo *SI;
	SufDesc *grefd, *pd, *refd;
	cde *Cde;
	cde2 *Cde2;
	cgrad **Cgrd, *cg, *cg0, *cg1, **cgb, **cgp0, **cgx;
	const char *s;
	expr_n *en;
	int f, i, j, j0, k, k0, m, m0, n, niv0, nnc, nnv, nnv0, nnz, ns;
	int nsos, nsos1, nsos2, nsosnz, nsosnz1, p1;
	int *col1, *cs, *g, *ge, **gp, **gp0, **gp1, **gpe, *myp[2], *p;
	int *pri, *sospri, *sosbeg, *sosbeg1, *sosind, *v, *v0, *ve, *z, **zg;
	ograd *og;
	real LU[2], LU1[2], *a, *lu, *mysr, *sosref, *sufref, *u, t, t1;
	size_t L;

	if (!VP || !(SI = (SOSinfo*)*VP))
		return 0;
	/* consistency check... */
	if (asl != SI->asl) {
		fprintf(Stderr, "Botched VP argument to sos_finish\n");
		mainexit_ASL(1);
		}
	*VP = 0;
	col1	= SI->col1;
	gp0	= SI->gp0;
	gp1	= SI->gp1;
	gpe	= SI->gpe;
	grefd	= SI->grefd;
	nnc	= SI->nnc;
	nnv	= SI->nnv;
	nnz	= SI->nnz;
	nsos	= SI->nsos;
	nsos1	= SI->nsos1;
	nsos2	= SI->nsos2;
	nsosnz	= SI->nsosnz;
	nsosnz1	= SI->nsosnz1;
	refd	= SI->refd;
	v0	= SI->v0;
	ve	= SI->ve;
	z	= SI->z;

	if (pl_bigM <= 0.) {
		if ((s = getenv("pl_bigM")))
			pl_bigM = strtod(s,0);
		if (pl_bigM <= 0.)
			pl_bigM = 1e6;
		mpl_bigM = -pl_bigM;
		}

	if (nsosnz_p)
		*nsosnz_p = nsosnz;
	L = nsos + 2*(nsos + 1)*sizeof(int)
		+ nsosnz*(sizeof(int) + sizeof(double));
	if (!copri) {
		if (sospri_p)
			*sospri_p = 0;
		sospri_p = 0;
		}
	else if (sospri_p)
		L += nsos*sizeof(int);
	mysr = 0;
	Cgrd = Cgrad;
	if (!sosref_p || !sosind_p || !sosbeg_p) {
		if (sosref_p)
			*sosref_p = 0;
		if (sosind_p)
			*sosind_p = 0;
		if (sosbeg_p)
			*sosbeg_p = 0;
		if (sospri_p) {
			*sospri_p = 0;
			sospri_p = 0;
			}
		flags |= ASL_suf_sos_explict_free;
		sosref_p = &mysr;
		sosind_p = &myp[0];
		sosbeg_p = &myp[1];
		}
	sosref = *sosref_p = flags & ASL_suf_sos_explict_free
				? (real*)Malloc(L)
				: (real*)M1alloc(L);
	sosind = *sosind_p = (int *)(sosref + nsosnz);
	sosbeg = *sosbeg_p = sosind + nsosnz;
	sosbeg1 = sosbeg + nsos + 1;	/* scratch */
	sospri = 0;
	if (sospri_p)
		sospri = *sospri_p = (int*)(sosbeg1 + nsos + 1);
	memset(sosbeg, 0, (nsos+1)*sizeof(int));
	f = Fortran;
	cg0 = 0;	/* silence buggy "not-initialized" warning */
	cgp0 = 0;	/* ditto */
	n = n_var;	/* ditto */
	k0 = nnv0 = 0;	/* ditto */
	m0 = n_con;
	if (nsos1) {
		LUge[0] = 0.;
		LUge[1] = Infinity;
		CI.z = z;
		CI.lu = LUv;
		CI.u = Uvx;
		CI.clu = LUrhs;
		CI.cu = Urhsx;
		CI.m = m0;
		CI.nbin = 0;
		if (Cgrd) {
			CI.cgp = Cgrd + m0;
			CI.cg = (cgrad*)M1alloc(nnz*sizeof(cgrad));
			}
		else {
			CI.cgp = (cgrad**)Malloc(nnc*sizeof(cgrad*));
			CI.cg = (cgrad*)Malloc(nnz*sizeof(cgrad));
			}
		cg0 = CI.cg;
		cgp0 = CI.cgp;
		pri = 0;
		if (sospri) {
			memset(sospri, 0, nsos1*sizeof(int));
			if ((pd = suf_get("priority",
					ASL_Sufkind_var | ASL_Sufkind_input)))
				pri = pd->u.i;
			}
		i = 0;
		for(gp = gp0; gp < gp1; gp++)
			sosbeg[i++] = ((gp[1] - gp[0]) >> 1) - 1;
		for(; gp < gpe; gp++)
			sosbeg[i++] = (gp[1] - gp[0]) >> 1;
		for(i = j = k = 0; i < nsos1; i++) { /* "k = 0" to omit an erroneous warning */
			k = sosbeg[i] + j;
			sosbeg[i] = j;
			j = k;
			}
		sosbeg[nsos1] = k;
		sufref = grefd->u.r;
		if (A_vals) {
			cs = A_colstarts;
			k = cs[i = n];
			j = n_var;
			while(i < j)
				cs[++i] = k;
			}
		p = get_vcmap_ASL(asl, ASL_Sufkind_var);
 		for(i = n, k = n + nnv; i < k; ++i)
			p[i] = -1;
		if (niv) {
			/* tell write_sol about moved integer variables */

			k = n + nnv;
			p1 = n;
			n -= niv;
			while(--k > p1)
				p[k] = p[k-nnv];
			for(i = nnv; i > 0; --i)
				p[k--] = -1;

			/* copy integer variables up */

			nnv0 = n + nnv;
			for(i = n_obj; --i >= 0; )
				for(og = Ograd[i]; og; og = og->next) {
					if (og->varno >= n)
						og->varno += nnv;
					}
			if (A_vals) {
				cs = A_colstarts;
				for(i = n + niv; --i >= n; )
					cs[i+nnv] = cs[i];
				j = cs[n];
				for(i = n + nnv; i > n;)
					cs[--i] = j;
				}
			else for(i = n_con - nnc; --i >= 0; )
				for(cg = Cgrd[i]; cg; cg = cg->next) {
					if (cg->varno >= n)
						cg->varno += nnv;
					}
			}
		niv0 = n;
		/* add SOS2 stuff */
		ns = 0;
		for(gp = gp0; gp < gp1; gp++) {
			cgx = newcon(&CI, 0);
			g = gp[0];
			ge = gp[1] - 4;
			/* Order on ref row values, then */
			/* average those values. */
			qsortv(g, ((int)(ge-g)>>1) + 2, 2*sizeof(int),
				refcomp, sufref);
			for(; g <= ge; g += 2) {
				j = g[1];
				sufref[j] += 0.5*(sufref[g[3]] + sufref[j]);
				}
			for(g = gp[0]; g <= ge; g += 2) {
				if ((i = g[1]) >= niv0)
					i += nnv;
				j = nonbinary(i, &CI, LU);
				if (g == gp[0]) {
					if (!j)
						k = i;
					else
						Bound(&CI, i, k = n++, LU);
					}
				else if (g == ge) {
					if (g[3] >= niv0)
						g[3] += nnv;
					if (nonbinary(g[3], &CI, LU1)) {
						Bound2(&CI, i,k0,k=n++,LU);
						Bound(&CI, g[3], k, LU1);
						}
					else
						Bound(&CI, i, g[3], LU);
					}
				else
					Bound2(&CI, i, k0, k = n++, LU);
				newcoef(&CI, &cgx, k, 1.);
				sosref[ns] = sufref[g[1]];
				sosind[ns++] = k + f;
				k0 = k;
				}
			*cgx = 0;
			}
		/* add SOS1 stuff */
		for(; gp < gpe; gp++) {
			cgx = newcon(&CI, 0);
			for(g = gp[0], ge = gp[1]; g < ge; g += 2) {
				if (g[1] >= niv0)
					g[1] += nnv;
				if (nonbinary(k = g[1], &CI, LU))
					Bound(&CI, g[1], k = n++, LU);
				newcoef(&CI, &cgx, k, 1.);
				sosref[ns] = sufref[g[1]];
				sosind[ns++] = k + f;
				}
			*cgx = 0;
			}
		if (pri)
			for(i = 0, gp = gp0; gp < gpe; gp++) {
				j = 0;
				for(g = gp[0], ge = gp[1]; g < ge; g += 2) {
					if (j < pri[g[1]])
						j = pri[g[1]];
					}
				sospri[i++] = j;
				}

		/* reorder each SOS1 set by sosref */
		p = SI->g0;
		for(i = k = 0; i < nsos1; i++) {
			j = j0 = k;
			k = sosbeg[i];
			t = sosref[j++];
			while(j < k) {
				t1 = sosref[j++];
				if (t1 < t) {
					reorder(sosind, sosref, j0, k, p);
					break;
					}
				t = t1;
				}
			}

		if (!nsos2)
			goto sosbeg_adjust;
		nsos -= nsos1;
		sosbeg += nsos1;
		*sosbeg = 0;
		sosind += nsosnz1;
		sosref += nsosnz1;
		if (sospri)
			sospri += nsos1;
		}
	for(v = col1; v < ve; )
		if (((i = *v++) & 3) == 2 && !sosbeg[j = i>>4]++) {
			if (sospri)
				sospri[j-1] = copri[(i & 4) >> 2];
			}
	for(j = 0, i = 1; i <= nsos; i++) {
		k = sosbeg[i] + j;
		sosbeg1[i] = sosbeg[i] = j;
		j = k;
		}
	sufref = refd->u.r;
	for(v = col1; v < ve; v++)
		if (((i = *v) & 3) == 2)
			sosind[sosbeg[i>>4]++] = (v - v0) + f;
	for(v = col1; v < ve; v++)
		if ((i = *v) && !(i & 2)) {
			j = i >> 4;
			if ((k = sosbeg1[j]++) < sosbeg[j]) {
				k0 = v - v0;
				sosref[k] = sufref[k0] +
					0.5*(sufref[k0+1] - sufref[k0]);
				}
			}
	if (nsos1)
		for(i = 0; i <= nsos; i++)
			sosbeg[i] += nsosnz1;
	nsos += nsos1;

 sosbeg_adjust:
	if (nsos1) {
		k = (int)(CI.cg - cg0);
		debugchk("n", n_var - niv - CI.nbin, n, 1);
		debugchk("m", n_con, CI.m, 0);
		debugchk("nnz", nnz, k, 0);

		if (niv && CI.nbin) {

			/* Correct overestimation of number */
			/* of new binary variables. */

			p = asl->i.vmap;
			for(i = n_var - niv - CI.nbin, k = i + niv; i < k; ++i)
				p[i] = p[i + CI.nbin];
			for(i = n_obj; --i >= 0; )
				for(og = Ograd[i]; og; og = og->next) {
					if (og->varno >= nnv0)
						og->varno -= CI.nbin;
					}
			if (A_vals) {
				cs = A_colstarts;
				for(i = n, j = n + niv; i <= j; i++)
					cs[i] = cs[i+CI.nbin];
				}
			else for(i = n_con - nnc; --i >= 0; )
				for(cg = Cgrd[i]; cg; cg = cg->next) {
					if (cg->varno >= nnv0)
						cg->varno -= CI.nbin;
					}
			}
		if (CI.nbin) {
			n_var -= CI.nbin;
			c_vars -= CI.nbin;
			o_vars -= CI.nbin;
			nbv -= CI.nbin;
			nnv -= CI.nbin;
			if ((zg = zerograds)) {
				j = n_var;
				for(i = n_obj; i > 0; --i) {
					for(p = *zg++; *p >= 0 && *p < j; p++);
					*p = -1;
					}
				}
			}

		nnc = CI.m - m0;
		n_conjac[1] = n_con = CI.m;
		asl->i.nzc_ -= nnz - k;
		nnz = k;

		en = (expr_n *)mem_ASL(asl, sizeof(expr_n));
		en->v = 0.;
		en->op = (efunc_n *)f_OPNUM;
		i = m0;
		switch(asl->i.ASLtype) {
		  case ASL_read_f:
		  case ASL_read_fg:
			Cde = ((ASL_fg*)asl)->I.con_de_;
			goto more_Cde;
		  case ASL_read_pfg:
			Cde = ((ASL_pfg*)asl)->I.con_de_;
 more_Cde:
			while(i < CI.m)
				Cde[i++].e = (expr*)en;
			break;
		  case ASL_read_fgh:
			Cde2 = ((ASL_fgh*)asl)->I.con2_de_;
			goto more_Cde2;
		  case ASL_read_pfgh:
			Cde2 = ((ASL_pfgh*)asl)->I.con2_de_;
 more_Cde2:
			while(i < CI.m)
				Cde2[i++].e = (expr2*)en;
		  }
		j = asl->i.n_con1;
		asl->i.n_con1 = j + nnc;
		if ((p = asl->i.cmap))
			for(i = m0; i < CI.m; ++i)
				p[i] = j++;
		}
	if (f)
		for(sosbeg = *sosbeg_p, i = 0; i <= nsos; i++)
			sosbeg[i] += f;

	if (!nnc) /* no new constraints */
		goto freeup;

	if (nnv) {
		j = asl->i.n_var1;
		asl->i.n_var1 = j + nnv;
		if ((p = asl->i.vminv)) {
			i = n_var;
			k = i + nnv;
			do p[i++] = j++; while(j < k);
			}
		lu = LUv;
		k = n;	/* n == n_var - niv */
		i = n + niv;
		if ((u = Uvx)) {
			while(i > k) {
				--i;
				lu[i] = lu[i-nnv];
				u[i] = u[i-nnv];
				}
			for(k -= nnv; --i >= k; ) {
				lu[i] = 0.;
				u[i] = 1.;
				}
			}
		else {
			for(i <<= 1, j = nnv << 1,  k <<= 1; --i >= k;)
				lu[i] = lu[i-j];
			k -= 2*nnv;
			do {
				lu[i--] = 1.;
				lu[i--] = 0.;
				}
				while(i >= k);
			}
		}
	n = n_var;
	if (Cgrd) {

		/* (re)compute goff fields */

		m = m0 + nnc;
		n += nnv;
		p = (int*)Malloc(L = n*sizeof(int));
		memset(p, 0, L);
		for(i = 0; i < m; i++)
			for(cg = Cgrd[i]; cg; cg = cg->next)
				p[cg->varno]++;
		for(i = k = 0; i < n; i++) {
			j = p[i];
			p[i] = k;
			k += j;
			}
		for(i = 0; i < m; i++)
			for(cg = Cgrd[i]; cg; cg = cg->next)
				cg->goff = p[cg->varno]++;
		free(p);
		}
	else {
		/* insert new nonzeros into A_val, A_rownos, A_colstarts */

		cgb = (cgrad**)Malloc(L = n*sizeof(cgrad*));
		memset(cgb, 0, L);
		for(i = 0; i < nnc; i++)
			for(cg = cgp0[i]; cg; cg = cg1) {
				cg1 = cg->next;
				j = cg->varno;
				cg->varno = i + m0;
				cg->next = cgb[j];
				cgb[j] = cg;
				}
		cs = A_colstarts;
		p = A_rownos;
		a = A_vals;
		k = asl->i.nzc_;
		for(i = n; i > 0; ) {
			j = cs[i];
			cs[i--] = k + f;
			for(cg = cgb[i]; cg; cg = cg->next) {
				a[--k] = cg->coef;
				p[k] = cg->varno + f;
				--nnz;
				}
			if (nnz <= 0)
				break;
			j -= cs[i];
			while(j-- > 0) {
				j0 = --k - nnz;
				a[k] = a[j0];
				p[k] = p[j0];
				}
			}
		free(cgb);
		free(cg0);
		free(cgp0);
		}
 freeup:
	if (mysr)
		free(mysr);
	for(i = 5; i--; )
		if ((vp = SI->zap[i])) {
			free(*vp);
			*vp = 0;
			}
	n_var += nnv;
	c_vars += nnv;
	o_vars += nnv;
	nbv += nnv;
	n_con += nnc;
	return nsos;
	}
