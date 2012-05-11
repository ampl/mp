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

/* nl_obj(n) = 1 if objective n is nonlinear, 0 otherwise. */
#define SKIP_NL2_DEFINES
#include "nlp.h"
#include "nlp2.h"
#include "asl_pfg.h"
#include "asl_pfgh.h"

 int
#ifdef KR_headers
nl_obj_ASL(asl, n) ASL *asl; int n;
#else
nl_obj_ASL(ASL *asl, int n)
#endif
{
	expr_n *e;
	static char who[] = "nl_obj";

	if (!asl)
		badasl_ASL(asl,0,who);
	else if (asl->i.ASLtype < ASL_read_f || asl->i.ASLtype >ASL_read_pfgh)
		badasl_ASL(asl,ASL_read_f,who);

	if (n >= 0 && n < n_obj) {
		switch(asl->i.ASLtype) {
		  case ASL_read_fgh:
			e = (expr_n*)(((ASL_fgh*)asl)->I.obj2_de_ + n)->e;
			break;
		  case ASL_read_pfg:
			e = (expr_n*)(((ASL_pfg*)asl)->I.obj_de_ + n)->e;
			break;
		  case ASL_read_pfgh:
			e = (expr_n*)(((ASL_pfgh*)asl)->I.obj2_de_ + n)->e;
			break;
		  default:
			e = (expr_n*)(((ASL_fg*)asl)->I.obj_de_ + n)->e;
		  }
		if (e->op != f_OPNUM_ASL)
			return 1;
		}
	return 0;
	}
