/****************************************************************
Copyright (C) 1997 Lucent Technologies
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

/* Sample mex function (in MATLAB 4.x format) for getting functions,
   gradients, and sparse Hessians from an AMPL .nl file.  Start with

	[x,bl,bu,v,cl,cu] = spamfunc('stub')

   to read in a problem (discarding the previous problem, if any).
   The return values are:

	x = primal initial guess
	v = dual initial guess
	bl, bu = lower and upper bounds on x
	cl, cu = lower and upper bounds on c (the constraint bodies).

   Then

	[f,c] = spamfunc(x,0)

   gives the function (f) and constraint bodies (c) at x;

	[g,Jac] = spamfunc(x,1)

   gives the gradient g of f, the Jacobian matrix J of c, and

	W = spamfunc(v)

   gives the Hessian W of the Lagrangian function L = f + v*c
   (at the last x at which spamfunc(x,1) was called).

   After finding optimal values for x and v,

	spamfunc('solution message',x,v)

   to write a stub.sol file.
*/

#include "mex.h"
#undef printf
#include "asl_pfgh.h"

static char msgbuf[256];

 static real*
sizechk(Matrix *mp, char *who, fint m)
{
	int m1, n1;
	m1 = mxGetM(mp);
	n1 = mxGetN(mp);
	if (m1 != m || n1 != 1 && m1) {
		sprintf(msgbuf,
			"Expected %s to be %d x 1 rather than %d x %d\n",
			who, m, m1, n1);
		mexErrMsgTxt(msgbuf);
		}
	return mxGetPr(mp);
	}

 static void
at_end(void)
{
	if (cur_ASL)
		ASL_free(&cur_ASL);
	}

 static void
usage(void)
{
	mexErrMsgTxt("spamfunc usage:\n\n\
	[x,bl,bu,v,cl,cu] = spamfunc('stub')\nor\n\
	[f,c] = spamfunc(x,0)\nor\n\
	[g,Jac] = spamfunc(x,1)\nor\n\
	W = spamfunc(v)\nor\n\
	spamfunc('solution message',x,v)\nwith\n\
	x = primal, v = dual variables");
	}

 void
mexFunction(int nlhs, Matrix **plhs, int nrhs, Matrix **prhs)
{
	FILE *nl;
	char *buf1, buf[512], *what;
	static fint n, nc, nz;
	fint nerror;
	real *J1, *W, *c, *f, *g, *v, *t, *x;
	static real *J;
	cgrad *cg, **cgp;
	static size_t Jsize;
	Jmp_buf err_jmp0;
	ASL_pfgh *asl = (ASL_pfgh*)cur_ASL;
	static fint nhnz;
	static real *Hsp;
	real *H, *He;
	int *Ir, *Jc;
	fint *hcs, *hr, i;

	if (nrhs == 1 && mxIsString(prhs[0])) {
		if (nlhs != 6)
			usage();
		if (mxGetString(prhs[0], buf1 = buf, sizeof(buf)))
			mexErrMsgTxt("Expected 'stub' as argument\n");
		at_end();
		mexAtExit(at_end);
		asl = (ASL_pfgh*)ASL_alloc(ASL_read_pfgh);
		return_nofile = 1;
		if (!(nl = jac0dim(buf1,strlen(buf)))) {
			sprintf(msgbuf, "Can't open %.*s\n",
				sizeof(msgbuf)-20, buf);
			mexErrMsgTxt(msgbuf);
			}
		if (n_obj <= 0)
			printf("Warning: objectve == 0\n");
		n = n_var;
		nc = n_con;
		nz = nzc;
		J = (real *)M1alloc(nz*sizeof(real));
		X0 = mxGetPr(plhs[0] = mxCreateFull(n, 1, REAL));
		LUv = mxGetPr(plhs[1] = mxCreateFull(n, 1, REAL));
		Uvx = mxGetPr(plhs[2] = mxCreateFull(n, 1, REAL));
		pi0 = mxGetPr(plhs[3] = mxCreateFull(nc, 1, REAL));
		LUrhs = mxGetPr(plhs[4] = mxCreateFull(nc, 1, REAL));
		Urhsx = mxGetPr(plhs[5] = mxCreateFull(nc, 1, REAL));
		pfgh_read(nl, ASL_findgroups);
		Jsize = nc*n*sizeof(real);

		/* Arrange to compute the whole sparese Hessian */
		/* of the Lagrangian function (both triangles). */

		nhnz = sphsetup(0, 0, nc > 0, 0);
		Hsp = (real *)M1alloc(nhnz*sizeof(real));
		return;
		}

	if (!filename)
		mexErrMsgTxt("spamfunc(\"stub\") has not been called\n");
	nerror = -1;
	err_jmp1 = &err_jmp0;
	if (nlhs == 2) {
		if (nrhs != 2)
			usage();
		x = sizechk(prhs[0],"x",n);
		t = sizechk(prhs[1],"0 or 1", 1);
		if (t[0] == 0.) {
			f = mxGetPr(plhs[0] = mxCreateFull(1, 1, REAL));
			c = mxGetPr(plhs[1] = mxCreateFull(nc, 1, REAL));
			if (setjmp(err_jmp0.jb)) {
				sprintf(msgbuf, "Trouble evaluating %s\n",
					what);
				mexErrMsgTxt(msgbuf);
				}
			what = "f";
			*f = objval(0, x, &nerror);
			what = "c";
			conval(x, c, &nerror);
			return;
			}
		g = mxGetPr(plhs[0] = mxCreateFull(n, 1, REAL));
		J1 = mxGetPr(plhs[1] = mxCreateSparse(nc, n, nz, REAL));
		what = "g";
		objgrd(0, x, g, &nerror);
		if (nc) {
			what = "J";
			jacval(x, J1, &nerror);
			Ir = mxGetIr(plhs[1]);
			memcpy(mxGetJc(plhs[1]), A_colstarts, (n+1)*sizeof(int));
			cgp = Cgrad;
			for(i = 0; i < nc; i++)
				for(cg = *cgp++; cg; cg = cg->next)
					Ir[cg->goff] = i;
			}
		return;
		}
	if (nlhs == 0 && nrhs == 3) {
		/* eval2('solution message', x, v): x = primal, v = dual */
		if (!mxIsString(prhs[0]))
			usage();
		x = sizechk(prhs[1],"x",n);
		v = sizechk(prhs[2],"v",nc);
		if (mxGetString(prhs[0], buf, sizeof(buf)))
			mexErrMsgTxt(
			 "Expected 'solution message' as first argument\n");
		write_sol(buf, x, v, 0);
		return;
		}
	if (nlhs != 1 || nrhs != 1)
		usage();
	v = sizechk(prhs[0],"v",nc);
	W = mxGetPr(plhs[0] = mxCreateSparse(n, n, nhnz, REAL));

	what = "W";
	sphes(H = Hsp, 0, 0, v);

	/* Expand the Hessian lower triangle into the full Hessian... */

	Ir = mxGetIr(plhs[0]);
	Jc = mxGetJc(plhs[0]);
	hcs = sputinfo->hcolstarts;
	hr = sputinfo->hrownos;
	for(i = 0; i <= n; i++)
		Jc[i] = hcs[i];
	He = H + hcs[n];
	while(H < He) {
		*W++ = *H++;
		*Ir++ = *hr++;
		}
	}
