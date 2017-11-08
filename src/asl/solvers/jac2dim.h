/*******************************************************************
Copyright (C) 2016 AMPL Optimization, Inc.; written by David M. Gay.

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

#ifndef JAC2DIM_H_included
#define JAC2DIM_H_included
#ifndef NLP_H2_included
#include "nlp2.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif
 extern void hv2comp_ASL(ASL*, real *hv, real *p, int nobj, real *ow, real *y);
 extern void hv2compd_ASL(ASL*, real *hv, real *p, int co);
 extern varno_t hv2comps_ASL(ASL*, real *hv, real *p, int co, varno_t nz, varno_t *z);
 extern real con2ival_ASL(ASL*, int i, real *X, fint *ne);
 extern void con2grd_ASL (ASL*, int i, real *X, real *G, fint *nerror);
 extern void con2val_ASL (ASL*, real *X, real *F, fint *nerror);
 extern void jac2val_ASL (ASL*, real *X, real *J, fint *nerror);
 extern int  lcon2val_ASL(ASL*, int i, real *X, fint *ne);
 extern real obj2val_ASL (ASL*, int nobj, real *X, fint *nerror);
 extern void obj2grd_ASL (ASL*, int nobj, real *X, real *G, fint *ne);
 extern int  x2_check_ASL(ASL_fgh*, real*);
 extern int x2known_ASL (ASL*, real*, fint*);
#ifdef __cplusplus
	}
#endif

#endif /* JAC2DIM_H_included */
