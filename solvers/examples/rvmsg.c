/****************************************************************
Copyright (C) 1997-1998 Lucent Technologies
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
#define asl cur_ASL

 void
vivvals(char *sname, char *opname, char **argv, fint *IV, real *V)
{
	U_info uinfo;

	Oinfo.uinfo = (char*)&uinfo;
	uinfo.L = IV;
	uinfo.D = V;
	nprob = 1;
	Oinfo.sname = sname;
	Oinfo.bsname = sname;
	Oinfo.opname = opname;
	if (getopts(argv, &Oinfo))
		exit(1);
	--nprob;
	}

#define F 9
#define F0 12
#define NFCALL 5
#define NFCOV 51
#define NGCALL 29
#define NGCOV 52
#define NREDUC 5
#define PREDUC 6
#define RELDX 16

 char *
rvmsg(char *who, fint *IV, real *V, int wantcov)
{
	typedef struct { char *msg; int code; } Sol_info;
	static Sol_info solinfo[] = {
		{ "X-Convergence", 2 }, /* 3 */
		{ "Relative Function Convergence", 1 },
		{ "X- and Relative Function Convergence", 3 },
		{ "Absolute Function Convergence", 4 },
		{ "Singular Convergence", 5 },
		{ "False Convergence", 500 },
		{ "Function Evaluation Limit", 400 },
		{ "Iteration Limit", 401 },
		{ "STOPX", 401 }, /* 11 */
		{ "Initial f(x) cannot be computed", 501 }, /* 63 */
		{ "Bad parameters to ASSESS", 511 },
		{ "Gradient could not be computed", 502 }
		};
	static char buf[256];
	char *b;
	real f0, f1, nreldf, preldf;
	fint nfcov, ngcov, rc;

	if ((rc = IV[0]) >= 63)
		rc += 12 - 63;
	if (rc < 3 || rc > 14) {
		solve_result_num = 510;
		b = buf + Sprintf(buf, "%s:\t??? return code %ld", who, rc);
		}
	else {
		rc -= 3;
		solve_result_num = solinfo[rc].code;
		b = buf + Sprintf(buf, "%s:\t%s", who, solinfo[rc].msg);
		}
	if (IV[0] == 63)
		return buf;
	b += Sprintf(b, "; function = ");
	f1 = V[F];
	b += g_fmtop(b, f1);
	b += Sprintf(b, "\n\tRELDX = ");
	b += g_fmtp(b, V[RELDX], 3);
	nreldf = preldf = 0;
	f0 = V[F0];
	if (f0 < 0)
		f0 = -f0;
	if (f1 < 0)
		f1 = -f1;
	if (f0 < f1)
		f0 = f1;
	if (f0 > 0) {
		nreldf = V[NREDUC]/f0;
		preldf = V[PREDUC]/f0;
		}
	b += Sprintf(b, "; PRELDF = ");
	b += g_fmtp(b, preldf, 3);
	b += Sprintf(b, "; NPRELDF = ");
	b += g_fmtp(b, nreldf, 3);
	nfcov = ngcov = 0;
	if (wantcov) {
		nfcov = IV[NFCOV];
		ngcov = IV[NGCOV];
		}
	b += Sprintf(b, "\n\t%ld func. evals", IV[NFCALL] - nfcov);
	if (nfcov)
		b += Sprintf(b, " + %ld for covariance", nfcov);
	b += Sprintf(b, "%s%ld grad. evals", nfcov ? "\n\t" : "; ",
			IV[NGCALL] - ngcov);
	if (ngcov)
		b += Sprintf(b, " + %ld for covariance", ngcov);
	return buf;
	}

 char *
D_U(Option_Info *oi, keyword *kw, char *value)
{
	U_info *u = (U_info *)oi->uinfo;
	return Dval_ASL(oi, kw, value, u->D + Intcast kw->info);
	}

 char *
L_U(Option_Info *oi, keyword *kw, char *value)
{
	U_info *u = (U_info *)oi->uinfo;
	return Lval_ASL(oi, kw, value, u->L + Intcast kw->info);
	}
