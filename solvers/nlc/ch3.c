#include "math.h"
#include "errno.h"
#ifndef fint
#ifndef Long
#include "arith.h"	/* for Long */
#ifndef Long
#define Long long
#endif
#endif
#define fint Long
#endif
#ifndef real
#define real double
#endif
#ifdef __cplusplus
extern "C" {
#endif
 real acosh_(real *);
 real asinh_(real *);
 real acoshd_(real *, real *);
 real asinhd_(real *, real *);
 void in_trouble(char *, real);
 void in_trouble2(char *, real, real);
 void domain_(char *, real *, fint);
 void zerdiv_(real *);
 fint auxcom_[1] = { 2 /* nlc */ };
 fint funcom_[19] = {
	3 /* nvar */,
	1 /* nobj */,
	3 /* ncon */,
	9 /* nzc */,
	0 /* densejac */,

	/* objtype (0 = minimize, 1 = maximize) */

	0,

	/* colstarts */
	1,
	4,
	7,
	10,

	/* rownos */
	1,
	2,
	3,
	1,
	2,
	3,
	1,
	2,
	3 };

 real boundc_[1+6+6] /* Infinity, variable bounds, constraint bounds */ = {
		1.7e308,
		-1.7e308,
		1.7e308,
		-1.7e308,
		1.7e308,
		-1.7e308,
		1.7e308,
		-0.3333333333333333,
		-0.3333333333333333,
		-1.,
		-1.,
		1.,
		1.};

 real x0comn_[3] = {
		0.25,
		0.5,
		0.75 };

 static real pd[30];

 static void
funnelb(real *x)
{
	real dv[1];

	/*** funnel ***/

	dv[0] = pd[0]*2.;
	dv[0] *= 2.;
	dv[0] += pd[9]*2.;
	pd[15] = dv[0];

	/*** funnel ***/

	dv[0] = pd[1]*2.;
	dv[0] *= 2.;
	dv[0] += pd[10]*2.;
	pd[16] = dv[0];

	/*** funnel ***/

	dv[0] = pd[2]*2.;
	dv[0] *= 2.;
	dv[0] += pd[11]*2.;
	pd[17] = dv[0];

	/*** funnel ***/

	dv[0] = pd[3]*2.;
	dv[0] *= 2.;
	dv[0] += pd[12]*pd[15];
	pd[18] = dv[0];

	/*** funnel ***/

	dv[0] = pd[4]*2.;
	dv[0] *= 2.;
	dv[0] += pd[13]*pd[16];
	pd[19] = dv[0];

	/*** funnel ***/

	dv[0] = pd[5]*2.;
	dv[0] *= 2.;
	dv[0] += pd[14]*pd[17];
	pd[20] = dv[0];
	}
static real old_x[3];
static int xkind = -1;

 static int
xcheck(real *x)
{
	real *x0 = x, *x1 = old_x, *xe = x + 3;
	real v[2];
	errno = 0;
	if (xkind >= 0) {
		while(*x0++ == *x1++)
			if (x0 == xe)
				return 0;
		--x0, --x1;
		}
	do *x1++ = *x0++;
		while(x0 < xe);
	xkind = 0;

	/*** defined variable 1 ***/

	pd[0] = -1. + 2.*x[0];

	/*** defined variable 2 ***/

	pd[1] = -1. + 2.*x[1];

	/*** defined variable 3 ***/

	pd[2] = -1. + 2.*x[2];

	/*** defined variable 4 ***/

	v[0] = 2. * x[0];
	v[1] = v[0] - 1.;
	pd[9] = 2. * v[1];
	v[1] = pd[9] * pd[0];
	v[0] = v[1] + -1.;
	pd[3] = v[0];

	/*** defined variable 5 ***/

	v[0] = 2. * x[1];
	v[1] = v[0] - 1.;
	pd[10] = 2. * v[1];
	v[1] = pd[10] * pd[1];
	v[0] = v[1] + -1.;
	pd[4] = v[0];

	/*** defined variable 6 ***/

	v[0] = 2. * x[2];
	v[1] = v[0] - 1.;
	pd[11] = 2. * v[1];
	v[1] = pd[11] * pd[2];
	v[0] = v[1] + -1.;
	pd[5] = v[0];

	/*** defined variable 7 ***/

	v[0] = 2. * x[0];
	v[1] = v[0] - 1.;
	pd[12] = 2. * v[1];
	v[1] = pd[12] * pd[3];
	pd[6] = v[1];

	/*** defined variable 8 ***/

	v[0] = 2. * x[1];
	v[1] = v[0] - 1.;
	pd[13] = 2. * v[1];
	v[1] = pd[13] * pd[4];
	pd[7] = v[1];

	/*** defined variable 9 ***/

	v[0] = 2. * x[2];
	v[1] = v[0] - 1.;
	pd[14] = 2. * v[1];
	v[1] = pd[14] * pd[5];
	pd[8] = v[1];
	return 1;
	}
 real
feval_(fint *nobj, fint *needfg, real *x, real *g)
{
	real v[9], dv[9];
	fint wantfg = *needfg;
	if (xcheck(x) && wantfg == 2)
		wantfg = 3;
	if (wantfg & 2) {
		funnelb(x);
		}

	if (wantfg & 1) {

	/*** defined variable 10 ***/

	v[6] = pd[6] + 1.;
	v[6] = v[6] - 2.*x[0];

	/*** defined variable 11 ***/

	v[7] = pd[7] + 1.;
	v[7] = v[7] - 2.*x[1];

	/*** defined variable 12 ***/

	v[8] = pd[8] + 1.;
	v[8] = v[8] - 2.*x[2];

  /***  objective ***/

	pd[21] = pd[0] + pd[1];
	pd[21] += pd[2];
	v[3] = 0.3333333333333333 * pd[21];
	pd[22] = v[3] * v[3];
	pd[23] = v[3] + v[3];
	v[3] = 0.5 * pd[22];
	pd[24] = pd[3] + pd[4];
	pd[24] += pd[5];
	v[4] = 0.3333333333333333 * pd[24];
	v[5] = v[4] - -0.3333333333333333;
	pd[25] = v[5] * v[5];
	pd[26] = v[5] + v[5];
	v[5] = 0.5 * pd[25];
	v[3] += v[5];
	pd[27] = v[6] + v[7];
	pd[27] += v[8];
	v[5] = 0.3333333333333333 * pd[27];
	pd[28] = v[5] * v[5];
	pd[29] = v[5] + v[5];
	v[5] = 0.5 * pd[28];
	v[3] += v[5];
	;}

	if (wantfg & 2) {
	dv[0] = 0.5*pd[29];
	dv[0] *= 0.3333333333333333;
	dv[1] = dv[0];
	dv[2] = dv[0];
	dv[3] = 0.5*pd[26];
	dv[3] *= 0.3333333333333333;
	dv[4] = dv[3];
	dv[5] = dv[3];
	dv[6] = 0.5*pd[23];
	dv[6] *= 0.3333333333333333;
	dv[7] = dv[6];
	dv[8] = dv[6];
	g[2] = -dv[1]*2.;
	g[1] = -dv[2]*2.;
	g[0] = -dv[0]*2.;
	g[2] += dv[1]*pd[20];
	g[1] += dv[2]*pd[19];
	g[0] += dv[0]*pd[18];
	g[2] += dv[4]*pd[17];
	g[1] += dv[5]*pd[16];
	g[0] += dv[3]*pd[15];
	g[2] += dv[7]*2.;
	g[1] += dv[8]*2.;
	g[0] += dv[6]*2.;
	}

	return v[3];
}

 void
ceval_(fint *needfg, real *x, real *c, real *J)
{
	real v[2], dv[1];
	real t1;
	fint wantfg = *needfg;
	if (xcheck(x) && wantfg == 2)
		wantfg = 3;

	if (wantfg & 1) {

  /***  constraint 1  ***/

	v[0] = 0.3333333333333333 * pd[5];
	v[1] = 0.3333333333333333 * pd[4];
	v[0] += v[1];
	v[1] = 0.3333333333333333 * pd[3];
	v[0] += v[1];
	c[0] = v[0];

  /***  constraint 2  ***/

	v[0] = 0.3333333333333333 * pd[8];
	v[1] = 0.3333333333333333 * pd[7];
	v[0] += v[1];
	v[1] = 0.3333333333333333 * pd[6];
	v[0] += v[1];
	t1 = v[0] + -0.6666666666666666*x[0];
	t1 += -0.6666666666666666*x[1];
	t1 += -0.6666666666666666*x[2];
	c[1] = t1;

  /***  constraint 3  ***/

	t1 = 0.6666666666666666*x[0];
	t1 += 0.6666666666666666*x[1];
	t1 += 0.6666666666666666*x[2];
	c[2] = t1;
	;}
	if (wantfg & 2) {
	funnelb(x);

   /*** derivatives for constraint 1 ***/

	J[6] = 0.3333333333333333*pd[17];
	J[3] = 0.3333333333333333*pd[16];
	J[0] = 0.3333333333333333*pd[15];

   /*** derivatives for constraint 2 ***/

	J[7] = 0.3333333333333333*pd[20] + -0.6666666666666666;
	J[4] = 0.3333333333333333*pd[19] + -0.6666666666666666;
	J[1] = 0.3333333333333333*pd[18] + -0.6666666666666666;

   /*** derivatives for constraint 3 ***/

	J[2] = 0.6666666666666666;
	J[5] = 0.6666666666666666;
	J[8] = 0.6666666666666666;
	}
}
#ifdef __cplusplus
	}
#endif
