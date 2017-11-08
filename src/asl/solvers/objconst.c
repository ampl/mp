/*******************************************************************
Copyright (C) 2017 AMPL Optimization, Inc.; written by David M. Gay.

Permission to use, copy, modify, and distribute this software and its
documentation for any purpose and without fee is hereby granted,
provided that the above copyright notice appear in all copies and that
both that the copyright notice and this permission notice and warranty
disclaimer appear in supporting documentation.

The author and AMPL Optimization, Inc. disclaim all warranties with
regard to this software, including all implied warranties of
merchantability and fitness.  In no event shall the author be liable
for any special, indirect or consequential damages or any damages
whatsoever resulting from loss of use, data or profits, whether in an
action of contract, negligence or other tortious action, arising out
of or in connection with the use or performance of this software.
*******************************************************************/

/* For LPs (and IPs and MIPs), objconst(n) is the constant term
 * for objective n (first objective has n = 0).
 */
#define SKIP_NL2_DEFINES
#include "nlp.h"
#include "nlp2.h"
#include "asl_pfg.h"
#include "asl_pfgh.h"
#include "obj_adj.h"
#undef f_OPNUM
#include "r_opn0.hd"	/* for f_OPNUM */

 real
#ifdef KR_headers
objconst_ASL(asl, n) ASL *asl; int n;
#else
objconst_ASL(ASL *asl, int n)
#endif
{
	Objrep *r, **rp;
	expr_n *e;
	efunc_n *opnum = f_OPNUM_ASL;
	real c;
	static char who[] = "objconst";

	if (!asl)
		badasl_ASL(asl,0,who);
	else if (asl->i.ASLtype < ASL_read_f || asl->i.ASLtype >ASL_read_pfgh)
		badasl_ASL(asl,ASL_read_f,who);

	c = 0.;
	if (n >= 0 && n < n_obj) {
		if ((rp = asl->i.Or) && (r = rp[n]))
			c = r->c0a;
		switch(asl->i.ASLtype) {
		  case ASL_read_fgh:
			e = (expr_n*)(((ASL_fgh*)asl)->I.obj2_de_ + n)->e;
			break;
		  case ASL_read_pfg:
			e = (expr_n*)(((ASL_pfg*)asl)->I.obj_de_ + n)->e;
			opnum = (efunc_n*)f_OPNUM;
			break;
		  case ASL_read_pfgh:
			e = (expr_n*)(((ASL_pfgh*)asl)->I.obj2_de_ + n)->e;
			opnum = (efunc_n*)f_OPNUM;
			break;
		  default:
			e = (expr_n*)(((ASL_fg*)asl)->I.obj_de_ + n)->e;
		  }
		if (e->op == opnum || e->op == (efunc_n *)f_OPNUM)
			return e->v;
		}
	return c;
	}
