/****************************************************************
Copyright (C) 1997, 2000 Lucent Technologies
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

#include "asl.h"

#ifdef __cplusplus
extern "C" {
extern void
jacinc1_(fint *M, fint *N, fint *NO, fint *NZ, fint *JP,
	fint *JI, real *X, real *L, real *U,
	real *Lrhs, real *Urhs, real *Inf, fint *Objdir);
}
#endif

 void
jacinc1_(fint *M, fint *N, fint *NO, fint *NZ, fint *JP,
	fint *JI, real *X, real *L, real *U,
	real *Lrhs, real *Urhs, real *Inf, fint *Objdir)
{
	int n;
	cgrad *gr, **grp;
	ASL *asl = cur_ASL;

	mnnzchk_ASL(asl, M, N, *NZ, "jacinc");
	*Inf = Infinity;
	if ((n = n_con)) {
		LUcopy_ASL(n, Lrhs, Urhs, LUrhs);
		grp = Cgrad + n;
		for(; n > 0; --n)
			for(gr = *--grp; gr; gr = gr->next) {
				JI[gr->goff] = n;
				JP[gr->varno] = gr->goff + 1;
				}
		JP[c_vars] = nzjac + 1;
		}
	LUcopy_ASL(c_vars, L, U, LUv);
	memcpy(X, X0, asl->i.n_var0 * sizeof(real));
	for(n = n_obj; --n >= 0; )
		Objdir[n] = objtype[n];
	}
