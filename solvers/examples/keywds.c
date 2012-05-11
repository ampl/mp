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

#include "rvmsg.h"
#ifdef PSEDREAD
#include "asl_pfg.h"
#else
#ifdef PSHVREAD
#include "asl_pfgh.h"
real cweight = 1.;
#endif
#endif

#ifdef Std_dev
char *stddev_file;
int want_stddev;
#endif

int showstats;
int nprob = 1;

#ifdef PSHVREAD
int wantgr = 3;
int hesprint, sparsetype, wanthvc;
#endif

#ifdef Timing
int nTimes, nhTimes;
#endif

 static keyword
keywds[] = {	/* must be in alphabetical order */
 KW("afctol", D_U, 30, "absolute function convergence tolerance"),
 KW("bias", D_U, 42, "bias from dogleg to double dogleg step"),
 KW("cosmin", D_U, 46, "min. allowed abs. cosine in updating s"),
 KW("covprt", L_U, 13, "1 means print covariance matrix"),
 KW("covreq", L_U, 14, "nonzero means try to compute covariance at soln"),
#ifdef PSHVREAD
 KW("cweight", D_val, &cweight, "weight for squared constraints (aug. Lagr.)"),
#endif
 KW("d0init", D_U, 39, "if positive, value for all d0 array components"),
 KW("decfac", D_U, 21, "factor for decreasing radius if iv(toobig) .ne. 0"),
 KW("delta0", D_U, 43, "used in picking f.d. step for cov. comp. (func. and grad"),
 KW("dfac", D_U, 40, "d updates set d(i).ge.v(dfac)*d(i) for iv(dtype)=1"),
 KW("dinit", D_U, 37, "nonneg. means initialize d to v(dinit)"),
 KW("dltfdc", D_U, 41, "used in picking f.d. step for computing cov. (func. only)"),
 KW("dltfdj", D_U, 42, "used in picking f.d. step for estimating jacobian (nl2sn)"),
 KW("dradpr", L_U, 100, "1 means print drops and adds"),
 KW("dtinit", D_U, 38, "if positive, initialize dtol array to it"),
 KW("dtype", L_U, 15, "tells how scale vector d is to be chosen"),
 KW("epslon", D_U, 18, "frac. by which quad. model may exceed min. value"),
 KW("eta0", D_U, 41, "relative error in f (used by sgrad2)"),
 KW("fuzz", D_U, 44, "factor by which other model must improve on current one"),
#ifdef PSHVREAD
 KW("groups", I_val, &wantgr, "find group sep. struc. (1 = obj, 2 = constrs, dflt = 3)"),
 KW("hesprint", I_val, &hesprint, "print Hessian after -t output"),
 KW("hffactor", DA_val, voffset_of(ASL,p.hffactor), "factor (default 1) in Hessian funnel decision"),
#endif
 KW("huberc", D_U, 47, "constant in huber*s robust criterion function"),
#ifdef PSHVREAD
 KW("hvcomp", I_val, &wanthvc, "Hessian times vector computation of full Hessian"),
#endif
 KW("incfac", D_U, 22, "factor by which to increase radius"),
 KW("inith", L_U, 24, "0 means initialize scaled hessian approx. to i"),
 KW("inits", L_U, 24, "0 means initialize s to zero and start with gn"),
 KW("iterlim", L_U, 17, "max. no. of iterations (same as mxiter)"),
 KW("lmax0", D_U, 34, "max. 2-norm allowed for very first step"),
 KW("lmaxs", D_U, 35, "radius for singular convergence test"),
 KW("maxfwd", IA_val, voffset_of(ASL,p.maxfwd_), "# of partials to forward recur; default = 5"),
 KW("maxit", L_U, 17, "max. no. of iterations (same as mxiter)"),
#ifdef PSEDREAD
 KW("merge", IA_val, voffset_of(ASL_pfg,P.merge), "do not merge elements with the same internal variables"),
#endif
 KW("mxfcal", L_U, 16, "max. no. of function (i.e., residual) evaluations"),
 KW("mxiter", L_U, 17, "max. no. of iterations"),
#ifdef Timing
 KW("nhtimes", I_val, &nhTimes, "repetitions of Hessian timing loop"),
 KW("ntimes", I_val, &nTimes, "repetitions of non-Hessian timing loops"),
#endif
 KW("objno", I_val, &nprob, "objective choice: 1 (default) = 1st"),
 KW("outlev", L_U, 18, "0 means do not summarize iterations"),
 KW("parprt", L_U, 19, "1 means print nondefault and changed v values"),
 KW("phmnfc", D_U, 19, "step computed by lmstep and gqtst must have"),
 KW("phmxfc", D_U, 20, "see phmnfc"),
 KW("prunit", L_U, 20, "i/o unit used for printing by itsum, parck"),
 KW("rdfcmn", D_U, 23, "min. allowed value for v(radfac)"),
 KW("rdfcmx", D_U, 24, "max. allowed value for v(radfac)"),
 KW("rdreq", L_U, 56, "kind of regression diagnostics required"),
 KW("rfctol", D_U, 31, "relative function convergence tolerance"),
 KW("rlimit", D_U, 45, "largest allowed norm of residual (ow treat as overflow)"),
 KW("rsptol", D_U, 48, "tolerance for splitting residual (into qr and normal eq)"),
 KW("sctol", D_U, 36, "tolerance for singular convergence test"),
 KW("sigmin", D_U, 49, "minimum sigma allowed in huber robust regression"),
 KW("solprt", L_U, 21, "1 means print x, g, d at convergence"),
#ifdef PSHVREAD
 KW("sparse", I_val, &sparsetype, "0 (default) = duthes; 1 = sputhes; 2 = fullhes"),
#endif
 KW("statpr", L_U, 22, "1 means print summary statistics at convergence"),
 KW("stats", I_val, &showstats, "show timing and memory statistics"),
#ifdef Std_dev
 KW("stddev", I_val, &want_stddev, "show solution and stddev"),
 KW("stddev_file", C_val, &stddev_file, "file for stddev (as AMPL param)"),
#endif
 KW("timing", I_val, &showstats, "show timing -- same as showstats"),
 KW("tuner1", D_U, 25, "used by assst to decide if func. dec. was tiny"),
 KW("tuner2", D_U, 26, "used by assst to decide if func. dec. was large"),
 KW("tuner3", D_U, 27, "used to decide if radius should be increased"),
 KW("tuner4", D_U, 28, "used to decide if predicted gradient change is"),
 KW("tuner5", D_U, 29, "also used to decide whether to increase radius"),
 KW("wantsol", WS_val, 0, WSu_desc_ASL+5),
 KW("x0prt", L_U, 23, "1 means print initial x"),
 KW("xctol", D_U, 32, "x-convergence tolerance (on v(reldx))"),
 KW("xftol", D_U, 33, "false convergence tolerance (on v(reldx))")
 };

Option_Info Oinfo = { 0, 0, 0, keywds, nkeywds, 1 };

/* To suppress by default echoing (by write_sol) of the solution
   message with invocations without -AMPL, change the above line to

	Option_Info Oinfo = { 0, 0, 0, keywds, nkeywds, 1,0,0,0,0,0,0,2 };

*/
