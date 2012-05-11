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

#ifndef RVMSG_H_included
#define RVMSG_H_included

#ifndef GETSTUB_H_included
#include "getstub.h"
#endif

 typedef struct U_info { real *D; fint *L; } U_info;
 typedef int (*U_fp)(VOID);

#ifdef __cplusplus
extern "C" {
#endif
extern Kwfunc D_U, L_U;
extern void divset_ ANSI((fint *L1, fint *IV, fint *LIV, fint *LV, real *V));
extern char *rvmsg ANSI((char *who, fint *IV, real *V, int wantcov));
extern void vivvals ANSI((char *sname, char *ename, char **av, fint *iv, real *v));

extern int hesprint, nprob, showstats, sparsetype, wanthvc, wantgr;
extern int want_stddev;
extern Option_Info Oinfo;
extern int xargc;
extern char *stddev_file, **xargv;
#ifdef __cplusplus
	}
#endif

#endif /* RVMSG_H_included */
