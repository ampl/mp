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

#ifndef JACPDIM_H_included
#define JACPDIM_H_included
#include "asl_pfgh.h"

#define conpival conpival_ASL
#define conpgrd  conpgrd_ASL
#define conpval  conpval_ASL
#define jacpval  jacpval_ASL
#define lconpval  lconpval_ASL
#define objpgrd  objpgrd_ASL
#define objpval  objpval_ASL
#define hvpcomp  hvpcomp_ASL
#define xp2known xp2known_ASL

#ifdef __cplusplus
extern "C" {
#endif
 extern real conpival ANSI((ASL*, int nc, real *X, fint *ne));
 extern void conpgrd ANSI((ASL*, int nc, real *X, real *G, fint *nerror));
 extern void conpval ANSI((ASL*, real *X, real *F, fint *nerror));
 extern void jacpval ANSI((ASL*, real *X, real *JAC, fint *nerror));
 extern int  lconpval ANSI((ASL*, int nc, real *X, fint *ne));
 extern void objpgrd ANSI((ASL*, int nobj, real *X, real *G, fint *nerror));
 extern real objpval ANSI((ASL*, int nobj, real *X, fint *nerror));
 extern void hvpcomp ANSI((ASL*,real *hv,real *p,int nobj,real *ow,real *y));
 extern int xp_check_ASL ANSI((ASL_pfgh*,real*));
 extern void xp2known ANSI((ASL*, real*, fint*));
#ifdef __cplusplus
	}
#endif
#endif /* JACPDIM_H_included */
