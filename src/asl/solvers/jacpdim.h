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

#ifndef JACPDIM_H_included
#define JACPDIM_H_included
#include "asl_pfgh.h"

#define conpival conpival_ASL
#define conpgrd  conpgrd_ASL
#define conpval  conpval_ASL
#define jacpval  jacpval_ASL
#define lconpval lconpval_ASL
#define objpgrd  objpgrd_ASL
#define objpval  objpval_ASL
#define hvpcomp  hvpcomp_ASL
#define hvpcompd hvpcompd_ASL
#define hvpcomps hvpcomps_ASL
#define xp2known xp2known_ASL

#ifdef __cplusplus
extern "C" {
#endif
 extern real conpival(ASL*, int nc, real *X, fint *ne);
 extern real conpival_nomap_ASL(ASL*, int nc, real *X, fint *ne);
 extern void conpgrd(ASL*, int nc, real *X, real *G, fint *nerror);
 extern void conpgrd_nomap_ASL(ASL*, int nc, real *X, real *G, fint *nerror);
 extern void conpval(ASL*, real *X, real *F, fint *nerror);
 extern void jacpval(ASL*, real *X, real *JAC, fint *nerror);
 extern int  lconpval(ASL*, int nc, real *X, fint *ne);
 extern void objpgrd(ASL*, int nobj, real *X, real *G, fint *nerror);
 extern real objpval(ASL*, int nobj, real *X, fint *nerror);
 extern void hvpcomp(ASL*, real *hv, real *p, int nobj, real *ow, real *y);
 extern void hvpcompd(ASL*,real *hv, real *p, int co);
 extern varno_t hvpcomps(ASL*, real *hv, real *p, int co, varno_t nz, varno_t *z);
 extern int xp_check_ASL(ASL_pfgh*, real*);
 extern int xp2known(ASL*, real*, fint*);
#ifdef __cplusplus
	}
#endif
#endif /* JACPDIM_H_included */
