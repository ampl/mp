/****************************************************************
Copyright (C) 1997, 2000-2001 Lucent Technologies
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

#ifndef JAC2DIM_H_included
#define JAC2DIM_H_included
#ifndef NLP_H2_included
#include "nlp2.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif
 extern void hv2comp_ASL ANSI((ASL*,real*hv,real*p,int nobj,real*ow,real*y));
 extern real con2ival_ASL ANSI((ASL*, int i, real *X, fint *ne));
 extern void con2grd_ASL  ANSI((ASL*, int i, real *X, real *G, fint *nerror));
 extern void con2val_ASL  ANSI((ASL*, real *X, real *F, fint *nerror));
 extern void jac2val_ASL  ANSI((ASL*, real *X, real *J, fint *nerror));
 extern int  lcon2val_ASL ANSI((ASL*, int i, real *X, fint *ne));
 extern real obj2val_ASL  ANSI((ASL*, int nobj, real *X, fint *nerror));
 extern void obj2grd_ASL  ANSI((ASL*, int nobj, real *X, real *G, fint *ne));
 extern int  x2_check_ASL ANSI((ASL_fgh*, real*));
 extern void x2known_ASL  ANSI((ASL*, real*, fint*));
#ifdef __cplusplus
	}
#endif

#endif /* JAC2DIM_H_included */
