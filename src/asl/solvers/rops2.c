/****************************************************************
Copyright (C) 1997-1999, 2000-2001 Lucent Technologies
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

#include "nlp2.h"
#include "errchk.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef MULTIPLE_THREADS
#define MTNot_Used(x) Not_Used(x)
#else
#define MTNot_Used(x) /*nothing*/
#define asl cur_ASL
#endif

#define introuble(who, a, jv) introuble_ASL(asl, who, a, jv)
#define introuble2(who, a, b, jv) introuble2_ASL(asl, who, a, b, jv)
#define zero_div(L,s) zero_div_ASL(asl, L, s)

 real
f_OPPLUS(expr *e A_ASL)
{
	expr *L;
	/* e->dL = e->dR = 1.; */
	L = e->L.e;
	e = e->R.e;
	return (*L->op)(L K_ASL) + (*e->op)(e K_ASL);
	}

 static real
f_OPMINUS(expr *e A_ASL)
{
	expr *L;
	/* e->dL = 1.;  */
	/* e->dR = -1.; */
	L = e->L.e;
	e = e->R.e;
	return (*L->op)(L K_ASL) - (*e->op)(e K_ASL);
	}

 static real
f_OPMULT(expr *e A_ASL)
{
	expr *e1, *e2;
	e1 = e->L.e;
	e2 = e->R.e;
	return (e->dR = (*e1->op)(e1 K_ASL)) * (e->dL = (*e2->op)(e2 K_ASL));
	}

 static real
f_OPDIV(expr *e A_ASL)
{
	expr *e1;
	real L, R, dL, rv;
	e1 = e->L.e;
	L = (*e1->op)(e1 K_ASL);
	e1 = e->R.e;
	if (!(R = (*e1->op)(e1 K_ASL))
#ifdef WANT_INFNAN
	 && !L
#endif
		)
		zero_div(L, "/");
	rv = L / R;
	if (want_deriv) {
		e->dR = -rv * (dL = e->dL = 1. / R);
		/*e->dL2 = 0;*/	/* !! could special-case division to avoid */
				/* !! this and two multiplies in hv_back */
		e->dLR = -dL*dL;
		dL *= e->dR;
		e->dR2 = -(dL + dL);
		}
	return rv;
	}

 static real
f_OPREM(expr *e A_ASL)
{
	expr *e1;
	real L, R, rv;
	/* e->dL = 1.; */
	e1 = e->L.e;
	L = (*e1->op)(e1 K_ASL);
	e1 = e->R.e;
	R = (*e1->op)(e1 K_ASL);
	rv = fmod(L,R);
	if (errchk(rv))
		introuble2("fmod",L,R,1);
	else if (want_deriv) {
		e->dR = (rv - L) / R;
		e->dR2 = 0;
		}
	return rv;
	}

 real
f_OPPOW(expr *e A_ASL)
{
	expr *e1;
	real L, R, rv, xlog, xym1;
	e1 = e->L.e;
	L = (*e1->op)(e1 K_ASL);
	e1 = e->R.e;
	R = (*e1->op)(e1 K_ASL);
	rv = mypow(L,R);
	if (errchk(rv)) {
#ifdef WANT_INFNAN
		if (!L && R < 0.) {
			errno_set(0);
			if (want_deriv) {
				e->dL = e->dR = negInfinity;
				e->dL2 = e->dLR = e->dR2 = Infinity;
				}
			return Infinity;
			}
#endif
		introuble2("pow",L,R,1);
		}
	if (want_deriv) {
		if (L < 0.) {
 badpow:
			introuble2("pow'",L,R,2);
			}
		else if (L > 0.) {
			xym1 = rv / L;
			xlog = log(L);
			e->dL = R*xym1;
			e->dR = xlog * rv;
			e->dL2 = (R-1.)*(e->dL/L);
			e->dLR = xym1 * (1. + R*xlog);
			e->dR2 = xlog * e->dR;
			}
		else if (R > 1.) {
			e->dL = 0.;
			if (R >= 2.)
				goto dRzero;
			e->dR = 0.;
#ifdef WANT_INFNAN
			e->dL2 = Infinity;
			e->dLR = e->dR2 = 0.;
#else
			introuble2("pow\"",L,R,3);
#endif
			}
		else if (R == 1.) {
			e->dL = 1.;
 dRzero:
			e->dR = e->dL2 = e->dLR = e->dR2 = 0.;
			}
		else if (R == 0.) {
			e->dL = 0.;
			goto dRzero;
			}
		else {
			/* 0 < R < 1 */
#ifdef WANT_INFNAN
			e->dL = Infinity;
			e->dR = e->dR2 = 0.;
			e->dL2 = e->dLR = negInfinity;
#else
			goto badpow;
#endif
			}
		}
	return rv;
	}
 real
f_OP1POW(expr *e A_ASL)
	/* f_OPPOW for R = numeric constant */
{
	expr *e1;
	real L, R, rv;
	e1 = e->L.e;
	L = (*e1->op)(e1 K_ASL);
	R = ((expr_n *)e->R.e)->v;
	rv = mypow(L,R);
	if (errchk(rv)) {
#ifdef WANT_INFNAN
		if (!L && R < 0.) {
			errno_set(0);
			e->dL = negInfinity;
			return e->dL2 = Infinity;
			}
#endif
		introuble2("pow",L,R,1);
		}
	if (want_deriv) {
		if (L) {
			e->dL = R * (rv/L);
			e->dL2 = (R-1.)*(e->dL/L);
			}
		else if (R > 1.) {
			e->dL = 0.;
			if (R >= 2.)
				e->dL2 = 0.;
			else
#ifdef WANT_INFNAN
				e->dL2 = negInfinity;
#else
				introuble2("pow\"",L,R,3);
#endif
			}
		/* R == 1. does not arise here */
		else {
#ifdef WANT_INFNAN
			e->dL = negInfinity;
			e->dL2 = Infinity;
#else
			introuble2("pow'",L,R,2);
#endif
			}
		}
	return rv;
	}
 real
f_OP2POW(expr *e A_ASL)
	/* f_OPPOW for R = 2 */
{
	expr *e1;
	real L;
	e1 = e->L.e;
	L = (*e1->op)(e1 K_ASL);
	e->dL = L + L;
	return L*L;
	}

 real
f_OPCPOW(expr *e A_ASL)
		/* f_OPPOW for L = numeric constant */
{
	expr *e1;
	real L, R, rv;
	e1 = e->R.e;
	rv = mypow(L = e->L.en->v, R = (*e1->op)(e1 K_ASL));
	if (errchk(rv)) {
#ifdef WANT_INFNAN
		if (!L && R < 0.) {
			errno_set(0);
			e->dR = e->dR2 = 0;
			return Infinity;
			}
#endif
		introuble2("pow",L,R,1);
		}
	if (want_deriv) {
		if (L > 0.) {
			if (e->dL == 1)
				e->dL = log(L);	/* cache value */
			e->dR = e->dL * rv;
			e->dR2 = e->dR * e->dL;
			}
		else if (L == 0.)
			e->dR = e->dR2 = 0.;
		else {
#if defined(WANT_INFNAN) && defined(QNaN0)
			union { unsigned int x[2]; real r; } u;
			u.x[0] = QNaN0;
			u.x[1] = QNaN1;
			e->dR = e->dR2 = u.r;
#else
			introuble2("pow'",L,R,2);
#endif
			}
		}
	return rv;
	}

 static real
f_OPLESS(expr *e A_ASL)
{
	expr *e1;
	real L;
	e1 = e->L.e;
	L = (*e1->op)(e1 K_ASL);
	e1 = e->R.e;
	L -= (*e1->op)(e1 K_ASL);
	if (L < 0.)
		return e->dL = e->dR = 0;
	e->dL = 1.;
	e->dR = -1.;
	return L;
	}

 static real
f_MINLIST(expr *e0 A_ASL)
{
	de *d, *d1;
	real t, rv;
	expr *e1, *e2;
	derp *D;
	expr_va *e = (expr_va *)e0;

	d = e->L.d;
	e1 = e2 = d->e;
	rv = (*e1->op)(e1 K_ASL);
	for(d1 = d++; (e1 = d->e); d++) {
		t = (*e1->op)(e1 K_ASL);
		if (rv > t) {
			rv = t;
			d1 = d;
			e2 = e1;
			}
		}
	if ((D = e->R.D)) {
		D->a.rp = d1->dv.rp;
		D->next = d1->d;
		}
	e->val  = e2;
	e->vale = d1->ee;
	e->valf = d1->ef;
	return rv;
	}

 static real
f_MAXLIST(expr *e0 A_ASL)
{
	de *d, *d1;
	real t, rv;
	expr *e1, *e2;
	derp *D;
	expr_va *e = (expr_va *)e0;

	d = e->L.d;
	e1 = e2 = d->e;
	rv = (*e1->op)(e1 K_ASL);
	for(d1 = d++; (e1 = d->e); d++) {
		t = (*e1->op)(e1 K_ASL);
		if (rv < t) {
			rv = t;
			d1 = d;
			}
		}
	if ((D = e->R.D)) {
		D->a.rp = d1->dv.rp;
		D->next = d1->d;
		}
	e->val  = e2;
	e->vale = d1->ee;
	e->valf = d1->ef;
	return rv;
	}

 static real
f_FLOOR(expr *e A_ASL)
{
	/* e->dL = 0.; */
	e = e->L.e;
	return floor((*e->op)(e K_ASL));
	}

 static real
f_CEIL(expr *e A_ASL)
{
	/* e->dL = 0.; */
	e = e->L.e;
	return ceil((*e->op)(e K_ASL));
	}

 static real
f_ABS(expr *e A_ASL)
{
	real rv;
	expr *e1;
	e1 = e->L.e;
	rv = (*e1->op)(e1 K_ASL);
	if (rv < 0.) {
		e->dL = -1.;
		return -rv;
		}
	e->dL = 1.;
	return rv;
	}

 static real
f_OPUMINUS(expr *e A_ASL)
{
	/* e->dL = -1.; */
	e = e->L.e;
	return -(*e->op)(e K_ASL);
	}

 static real
f_OP_tanh(expr *e A_ASL)
{
	real rv, t, t1;
	expr *e1;
	e1 = e->L.e;
	rv = tanh(t = (*e1->op)(e1 K_ASL));
	if (errchk(rv))
		introuble("tanh",t,1);
	if (want_deriv) {
		t1 = cosh(t);
		if (errchk(t1))
			introuble("tanh'", t, 2);
		else {
			t1 = 1. / t1;
			e->dL2 = -(rv+rv) * (e->dL = t1*t1);
			}
		}
	return rv;
	}

 static real
f_OP_tan(expr *e A_ASL)
{
	real rv, t, t1;
	expr *e1;
	e1 = e->L.e;
	rv = tan(t = (*e1->op)(e1 K_ASL));
	if (errchk(rv))
		introuble("tan",t,1);
	if (want_deriv) {
		t1 = cos(t);
		if (errchk(t1) || !t1)
			introuble("tan'",t,2);
		else {
			t1 = 1. / t1;
			e->dL2 = (rv+rv) * (e->dL = t1*t1);
			}
		}
	return rv;
	}

 static real
f_OP_sqrt(expr *e A_ASL)
{
	real t, rv;
	expr *e1;
	e1 = e->L.e;
	t = (*e1->op)(e1 K_ASL);
	if (t < 0. || (rv = sqrt(t), errchk(rv))) {
		introuble("sqrt",t,1);
		rv = 0.; /* not reached */
		}
	if (want_deriv) {
		if (rv <= 0.)
			introuble("sqrt'",t,2);
		else {
			e->dL = 0.5 / rv;
			e->dL2 = -0.5 * e->dL / t;
			}
		}
	return rv;
	}

 static real
f_OP_sinh(expr *e A_ASL)
{
	real t, rv;
	expr *e1;
	e1 = e->L.e;
	rv = sinh(t = (*e1->op)(e1 K_ASL));
	if (errchk(rv))
		introuble("sinh",t,1);
	if (want_deriv) {
		e->dL = cosh(t);
		if (errchk(e->dL))
			introuble("sinh'",t,2);
		else
			e->dL2 = rv;
		}
	return rv;
	}

 static real
f_OP_sin(expr *e A_ASL)
{
	real t, rv;
	expr *e1;
	e1 = e->L.e;
	rv = sin(t = (*e1->op)(e1 K_ASL));
	if (errchk(rv))
		introuble("sin",t,1);
	if (want_deriv) {
		e->dL = cos(t);
		if (errchk(e->dL))
			introuble("sin'",t,2);
		else
			e->dL2 = -rv;
		}
	return rv;
	}

 static real
f_OP_log10(expr *e A_ASL)
{
	real t, rv;
	expr *e1;
	static real Le10;

	e1 = e->L.e;
	rv = log10(t = (*e1->op)(e1 K_ASL));
	if (errchk(rv))
		introuble("log10",t,1);
	if (want_deriv) {
		if (!Le10)
			Le10 = 1. / log(10.);
		e->dL = Le10 / t;
		e->dL2 = -e->dL / t;
		}
	return rv;
	}

 static real
f_OP_log(expr *e A_ASL)
{
	real t, rv;
	expr *e1;
	e1 = e->L.e;
	rv = log(t = (*e1->op)(e1 K_ASL));
	if (errchk(rv))
		introuble("log",t,1);
	if (want_deriv) {
		t = e->dL = 1. / t;
		e->dL2 = -t * t;
		}
	return rv;
	}

 static real
f_OP_exp(expr *e A_ASL)
{
	real t, rv;
	expr *e1;
	e1 = e->L.e;
	rv = e->dL = e->dL2 = exp(t = (*e1->op)(e1 K_ASL));
	if (errchk(rv)) {
		if (t >= 0.)
			introuble("exp",t,1);
		else {
			errno_set(0);
			rv = 0.;
			}
		}
	return rv;
	}

 static real
f_OP_cosh(expr *e A_ASL)
{
	real t, rv;
	expr *e1;
	e1 = e->L.e;
	rv = cosh(t = (*e1->op)(e1 K_ASL));
	if (errchk(rv))
		introuble("cosh",t,1);
	if (want_deriv) {
		e->dL = sinh(t);
		if (errchk(e->dL))
			introuble("cosh'",t,2);
		else
			e->dL2 = rv;
		}
	return rv;
	}

 static real
f_OP_cos(expr *e A_ASL)
{
	real t, rv;
	expr *e1;
	e1 = e->L.e;
	rv = cos(t = (*e1->op)(e1 K_ASL));
	if (errchk(rv))
		introuble("cos",t,1);
	if (want_deriv) {
		e->dL = -sin(t);
		if (errchk(e->dL))
			introuble("cos'",t,2);
		else
			e->dL2 = -rv;
		}
	return rv;
	}

 static real
f_OP_atanh(expr *e A_ASL)
{
	real rv, t, t1;
	expr *e1;

	e1 = e->L.e;
	t = (*e1->op)(e1 K_ASL);
	if (t <= -1. || t >= 1.) {
		errno_set(EDOM);
 bad_atanh:
		rv = 0.;
		introuble("atanh",t,1);
		}
	else {
		rv = 0.5*log((1. + t) / (1. - t));
		if (errchk(rv))
			goto bad_atanh;
		}
	if (want_deriv) {
		e->dL = t1 = 1. / (1. - t*t);
		e->dL2 = (t+t)*t1*t1;
		}
	return rv;
	}

 static real
f_OP_atan2(expr *e A_ASL)
{
	real L, R, rv, t, t1;
	expr *e1;

	e1 = e->L.e;
	L = (*e1->op)(e1 K_ASL);
	e1 = e->R.e;
	R = (*e1->op)(e1 K_ASL);
	rv = atan2(L,R);
	if (errchk(rv))
		introuble2("atan2",L,R,1);
	if (want_deriv) {
/*
		t = L / R;
		t1 = 1. / (1. + t*t);
		e->dL = t1 /= R;
		e->dR = -t*t1;
 */

		t = 1. / (L*L + R*R);
		e->dL =  t * R;
		e->dR = -t * L;
		t *= t;
		t1 = L*R;
		e->dL2 = -(e->dR2 = t * (t1+t1));
		e->dLR = t * (L*L - R*R);
		}
	return rv;
	}

 static real
f_OP_atan(expr *e A_ASL)
{
	real rv, t, t1;
	expr *e1;

	e1 = e->L.e;
	rv = atan(t = (*e1->op)(e1 K_ASL));
	if (errchk(rv))
		introuble("atan",t,1);
	if (want_deriv) {
		e->dL = t1 = 1. / (1. + t*t);
		e->dL2 = -(t+t)*t1*t1;
		}
	return rv;
	}

 static real
f_OP_asinh(expr *e A_ASL)
{
	expr *e1;
	real rv, s, t, t1, t2;

	e1 = e->L.e;
	t = (*e1->op)(e1 K_ASL);
	s = 1.;
	if (t < 0.)
		s = -1.;
	rv = log(s*t + (t1 = sqrt(t2 = t*t + 1.)));
	if (errchk(rv))
		introuble("asinh",t,1);
	if (want_deriv)
		e->dL2 = -(t/t2) * (e->dL = 1. / t1);
	return s*rv;
	}

 static real
f_OP_asin(expr *e A_ASL)
{
	expr *e1;
	real rv, t, t1;

	e1 = e->L.e;
	rv = asin(t = (*e1->op)(e1 K_ASL));
	if (errchk(rv))
		introuble("asin",t,1);
	if (want_deriv) {
		if ((t1 = 1. - t*t) <= 0.)
			introuble("asin'",t,2);
		else
			e->dL2 = t * (e->dL = 1. / sqrt(t1)) / t1;
		}
	return rv;
	}

 static real
f_OP_acosh(expr *e A_ASL)
{
	expr *e1;
	real rv, t, t1, t2;

	e1 = e->L.e;
	t = (*e1->op)(e1 K_ASL);
	if (t < 1.) {
		errno_set(EDOM);
 bad_acosh:
		rv = t1 = t2 = 0.;
		introuble("acosh",t,1);
		}
	else {
		rv = log(t + (t1 = sqrt(t2 = t*t - 1.)));
		if (errchk(rv))
			goto bad_acosh;
		}
	if (want_deriv) {
		if (t2 <= 0.)
			introuble("acosh'",t,2);
		else
			e->dL2 = -t * (e->dL = 1. / t1) / t2;
		}
	return rv;
	}

 static real
f_OP_acos(expr *e A_ASL)
{
	expr *e1;
	real rv, t, t1;

	e1 = e->L.e;
	rv = acos(t = (*e1->op)(e1 K_ASL));
	if (errchk(rv))
		introuble("acos",t,1);
	if (want_deriv) {
		if ((t1 = 1. - t*t) <= 0.)
			introuble("acos'",t,2);
		else
			e->dL2 = t * (e->dL = -1. / sqrt(1. - t*t)) / t1;
		}
	return rv;
	}

 static real
f_OPIFnl(expr *e A_ASL)
{
	expr_if *eif = (expr_if *)e;
	derp *D;

	e = eif->e;
	if ((*e->op)(e K_ASL)) {
		eif->val = e = eif->T;
		eif->vale = eif->Te;
		eif->valf  = eif->Tf;
		if ((D = eif->D)) {
			D->a.rp = eif->Tv.rp;
			D->next = eif->dT;
			}
		}
	else {
		eif->val = e = eif->F;
		eif->vale = eif->Fe;
		eif->valf  = eif->Ff;
		if ((D = eif->D)) {
			D->a.rp = eif->Fv.rp;
			D->next = eif->dF;
			}
		}
	return (*e->op)(e K_ASL);
	}

 static real
f_OPOR(expr *e A_ASL)
{
	expr *e2;
	e2 = e->R.e;
	e = e->L.e;
	return (*e->op)(e K_ASL) || (*e2->op)(e2 K_ASL) ? 1. : 0.;
	}

 static real
f_OPAND(expr *e A_ASL)
{
	expr *e2;
	e2 = e->R.e;
	e = e->L.e;
	return (*e->op)(e K_ASL) && (*e2->op)(e2 K_ASL) ? 1. : 0.;
	}

 static real
f_LT(expr *e A_ASL)
{
	expr *e2;
	e2 = e->R.e;
	e = e->L.e;
	return (*e->op)(e K_ASL) < (*e2->op)(e2 K_ASL) ? 1. : 0;
	}

 static real
f_LE(expr *e A_ASL)
{
	expr *e2;
	e2 = e->R.e;
	e = e->L.e;
	return (*e->op)(e K_ASL) <= (*e2->op)(e2 K_ASL) ? 1. : 0;
	}

 static real
f_EQ(expr *e A_ASL)
{
	expr *e2;
	e2 = e->R.e;
	e = e->L.e;
	return (*e->op)(e K_ASL) == (*e2->op)(e2 K_ASL) ? 1. : 0;
	}

 static real
f_GE(expr *e A_ASL)
{
	expr *e2;
	e2 = e->R.e;
	e = e->L.e;
	return (*e->op)(e K_ASL) >= (*e2->op)(e2 K_ASL) ? 1. : 0;
	}

 static real
f_GT(expr *e A_ASL)
{
	expr *e2;
	e2 = e->R.e;
	e = e->L.e;
	return (*e->op)(e K_ASL) > (*e2->op)(e2 K_ASL) ? 1. : 0;
	}

 static real
f_NE(expr *e A_ASL)
{
	expr *e2;
	e2 = e->R.e;
	e = e->L.e;
	return (*e->op)(e K_ASL) != (*e2->op)(e2 K_ASL) ? 1. : 0;
	}

 static real
f_OPNOT(expr *e A_ASL)
{
	e = e->L.e;
	return (*e->op)(e K_ASL) ? 0. : 1.;
	}

 static real
f_ANDLIST(expr *e A_ASL)
{
	expr **ep, **epe;

	ep = e->L.ep;
	epe = e->R.ep;
	do {
		e = *ep++;
		if ((*e->op)(e K_ASL) == 0.)
			return 0.;
		}
		while(ep < epe);
	return 1.;
	}

 static real
f_ORLIST(expr *e A_ASL)
{
	expr **ep, **epe;

	ep = e->L.ep;
	epe = e->R.ep;
	do {
		e = *ep++;
		if ((*e->op)(e K_ASL) != 0.)
			return 1.;
		}
		while(ep < epe);
	return 0.;
	}

#define f_OPIMPELSE f_OPIMPELSE_rops2_ASL
 real
f_OPIMPELSE(expr *e A_ASL)
{
	expr_if *eif = (expr_if *)e;

	e = eif->e;
	e = ((*e->op)(e K_ASL) != 0.) ? eif->T : eif->F;
	return (*e->op)(e K_ASL);
	}

 static real
f_OP_IFF(expr *e A_ASL)
{
	expr *e1;
	int a, b;

	e1 = e->L.e;
	a = (*e1->op)(e1 K_ASL) != 0.;
	e1 = e->R.e;
	b = (*e1->op)(e1 K_ASL) != 0.;
	return a == b ? 1. : 0.;
	}

 real
f_OPSUMLIST(expr *e A_ASL)
{
	expr **ep, **epe;
	real x;
	ep = e->L.ep;
	epe = e->R.ep;
	e = *ep;
	x = (*e->op)(e K_ASL);
	while(++ep < epe) {
		e = *ep;
		x += (*e->op)(e K_ASL);
		}
	return x;
	}

 real
f_OPPLTERM(expr *e A_ASL)
{
	int n, z;
	plterm *p;
	real *bs;
	real r, t;
	MTNot_Used(asl);

	p = e->L.p;
	n = p->n;
	z = p->z;
	bs = p->bs + z;
	z >>= 1;
	r = ((expr_v *)e->R.e)->v;
	if (r >= 0) {
		n -= z;
		if (n <= 1 || r <= bs[1])
			return r*(e->dL = *bs);
		for(t = bs[0]*bs[1]; --n > 1 && r > bs[3]; bs += 2)
			t += (bs[3]-bs[1])*bs[2];
		return t + (r-bs[1])*(e->dL = bs[2]);
		}
	if (z <= 0)
		return r*(e->dL = bs[0]);
	for(t = bs[0]*bs[-1]; --z > 0 && r < bs[-3]; bs -= 2) {
		t += bs[-2]*(bs[-3] - bs[-1]);
		}
	return t + (r - bs[-1])*(e->dL = bs[-2]);
	}

 char *
f_OPHOL(expr *e A_ASL)
{
	MTNot_Used(asl);
	return ((expr_h *)e)->sym;
	}

 real
f_OPVARVAL(expr *e A_ASL)
{
	MTNot_Used(asl);
	return ((expr_v *)e)->v;
	}

 char *
f_OPIFSYM(expr *e A_ASL)
{
	expr_if *eif = (expr_if *)e;
	e = eif->e;
	e = (*e->op)(e K_ASL) ? eif->T : eif->F;
	return (*(sfunc*)e->op)(e K_ASL);
	}

 real
f_OPFUNCALL(expr *e A_ASL)
{
	TMInfo T, *T1, *T1prev;
	arglist *al;
	argpair *ap, *ape;
	const char *s;
	expr_f *f = (expr_f *)e;
	func_info *fi = f->fi;
	real rv;

	for(ap = f->ap, ape = f->ape; ap < ape; ap++) {
		e = ap->e;
		*ap->u.v = (*e->op)(e K_ASL);
		}
	for(ap = f->sap, ape = f->sape; ap < ape; ap++) {
		e = ap->e;
		*ap->u.s = (*(sfunc*)e->op)(e K_ASL);
		}
	T.u.prev = 0;
	al = f->al;
	al->TMI = &T;
	al->Errmsg = 0;
	rv = (*fi->funcp)(al);
	errno_set(0);
	if ((s = al->Errmsg))
		fintrouble_ASL(asl, fi, s, &T);
	for(T1 = T.u.prev; T1; T1 = T1prev) {
		T1prev = T1->u.prev;
		free(T1);
		}
	return rv;
	}

 static real
f_OPintDIV(expr *e A_ASL)
{
	expr *e1;
	real L, R;

	e1 = e->L.e;
	L = (*e1->op)(e1 K_ASL);
	e1 = e->R.e;
	if (!(R = (*e1->op)(e1 K_ASL)))
		zero_div(L, " div ");
	return (L /= R) >= 0 ? floor(L) : ceil(L);
	}

 static real
f_OPprecision(expr *e A_ASL)
{
	expr *e1;
	real L, R;
	char buf[32];

	e1 = e->L.e;
	L = (*e1->op)(e1 K_ASL);
	e1 = e->R.e;
	R = (*e1->op)(e1 K_ASL);
	g_fmtp(buf, L, (int)R);
	return strtod(buf, (char **)0);
	}

 static real
Round(real x, int prec)
{
	char *b, *s, *s0, *se;
	int decpt, L, sign;
	char buf[96];

	if (!x)
		return x;
	s = dtoa(x, 3, prec, &decpt, &sign, &se);
	if (decpt == 9999) {
 zreturn:
		freedtoa(s);
		return x;
		}
	L = se - s;
	if (L <= 0) {
		x = 0.;
		goto zreturn;
		}
	if (L > 80)
		se = s + 80;
	s0 = s;
	b = buf;
	if (sign)
	*b++ = '-';
	*b++ = '.';
	while(s < se)
		*b++ = *s++;
	*b = 0;
	freedtoa(s0);
	if (decpt)
		snprintf(b, buf + sizeof(buf) - b, "e%d", decpt);
	return strtod(buf, (char **)0);
	}

 static real
f_OPround(expr *e A_ASL)
{
	expr *e1;
	real L, R;

	e1 = e->L.e;
	L = (*e1->op)(e1 K_ASL);
	e1 = e->R.e;
	R = (*e1->op)(e1 K_ASL);
	return Round(L, (int)R);
	}

 static real
f_OPtrunc(expr *e A_ASL)
{
	expr *e1;
	real L, R;
	int prec;

	e1 = e->L.e;
	L = (*e1->op)(e1 K_ASL);
	e1 = e->R.e;
	if (!(R = (*e1->op)(e1 K_ASL)))
		return L >= 0. ? floor(L) : ceil(L);
	R = Round(L, prec = (int)R);
	if (R != L) {
		R = 0.5*mypow(10., (real)-prec);
		R = Round(L > 0 ? L - R : L + R, prec);
		}
	return R;
	}


 static real
f_OPCOUNT(expr *e A_ASL)
{
	expr **ep, **epe;
	real x;
	ep = e->L.ep;
	epe = e->R.ep;
	e = *ep++;
	if ((x = (*e->op)(e K_ASL)))
		x = 1.;
	do {
		e = *ep++;
		if ((*e->op)(e K_ASL))
			x++;
		}
		while(ep < epe);
	return x;
	}

 typedef struct
jb_st {
	jmp_buf jb;
	} jb_st;

 static int
rcompj(const void *a, const void *b, void *v)
{
	jb_st *J;
	real t = *(real *)a - *(real *)b;

	if (!t) {
		J = (jb_st*)v;
		longjmp(J->jb, 1);
		}
	return t < 0 ? -1 : 1;
	}

 static real
f_OPALLDIFF(expr *e A_ASL)
{
	int n;
	expr **ep, **epe;
	jb_st J;
	real *r, *r1, rbuf[128], t;

	ep = e->L.ep;
	epe = e->R.ep;
	n = epe - ep;
	r = rbuf;
	if (n > sizeof(rbuf)/sizeof(real))
		r = (real*)Malloc(n*sizeof(real));
	r1 = r;
	while(ep < epe) {
		e = *ep++;
		*r1++ = (*e->op)(e K_ASL);
		}
	t = 1.;
	if (setjmp(J.jb)) {
		t = 0.;
		goto done;
		}
	qsortv(r, n, sizeof(real), rcompj, &J);
 done:
	if (r != rbuf)
		free(r);
	return t;
	}

 static real
f_OPSOMESAME(expr *e A_ASL)
{
	real t = f_OPALLDIFF(e K_ASL);
	return t == 0. ? 1. : 0.;
	}

 static real
f_OPNUMBEROF(expr *e A_ASL)
{
	expr **ep, **epe;
	real n, t;

	ep = e->L.ep;
	epe = e->R.ep;
	n = 0.;
	e = *ep++;
	t = (*e->op)(e K_ASL);
	while(ep < epe) {
		e = *ep++;
		if (t == (*e->op)(e K_ASL))
			n++;
		}
	return n;
	}

#define f_OPATLEAST f_LE
#define f_OPATMOST  f_GE
#define f_OPEXACTLY f_EQ
#define f_OPNOTATLEAST f_GT
#define f_OPNOTATMOST  f_LT
#define f_OPNOTEXACTLY f_NE
#define efunc efunc2
#define f_OPNUMBEROFs f_OPNUMBEROF /*TEMPORARY*/

 efunc *
r_ops[] = {
#include "r_op.hd"
	};

#ifdef __cplusplus
}
#endif
