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

/* qsortv -- qsort with an extra argument v (of type void*) that is
   passed to the comparison function, to facilitate complicated
   comparisons in environments with multiple threads.
   Derived (by dmg) from code by Jon Bentley and Doug McIlroy, with
   modifications to avoid comparing an element to itself.
   See "Engineering a Sort Function" by Jon L. Bentley and M. Douglas
   McIlroy, Software--Practice and Experience 23 (11), Nov. 1993,
   pp. 1249-1265.
 */

#ifdef KR_headers
typedef int cmpfunc();
typedef unsigned int size_t;
#else
typedef int cmpfunc(const void *, const void*, void*);
#include "stddef.h"
#ifdef __cplusplus
extern "C" void qsortv(void*, size_t, size_t, cmpfunc*, void*);
#endif
#endif

#define SWAPINIT(a, es) swaptype =         \
    ((char*)a-(char*)0 | es) % sizeof(long) ? 2 : es > sizeof(long) ? 1 : 0;
#define swapcode(TYPE, parmi, parmj, n) {  \
    register TYPE *pi = (TYPE *) (parmi);  \
    register TYPE *pj = (TYPE *) (parmj);  \
    do {                                   \
        register TYPE t = *pi;             \
        *pi++ = *pj;                       \
        *pj++ = t;                         \
    } while ((n -= sizeof(TYPE)) > 0);     \
}

 static void
#ifdef KR_headers
swapfunc(a, b, n, swaptype) char *a, *b; size_t n; int swaptype;
#else
swapfunc(char *a, char *b, size_t n, int swaptype)
#endif
{   if (swaptype <= 1) swapcode(long, a, b, n)
    else swapcode(char, a, b, n)
}
#define swap(a, b)                         \
    if (swaptype == 0) {                   \
        t = *(long*)(a);                   \
        *(long*)(a) = *(long*)(b);         \
        *(long*)(b) = t;                   \
    } else                                 \
        swapfunc(a, b, es, swaptype)

#define vecswap(a, b, n) if (n > 0) swapfunc(a, b, n, swaptype)

#define min(x, y) ((x)<=(y) ? (x) : (y))

 static char *
#ifdef KR_headers
med3(a, b, c, cmp, v) char *a, *b, *c; cmpfunc *cmp; char *v;
#else
med3(char *a, char *b, char *c, cmpfunc *cmp, void *v)
#endif
{	return (*cmp)(a,b,v) < 0 ?
		  ((*cmp)(b,c,v) < 0 ? b : (*cmp)(a,c,v) < 0 ? c : a)
		: ((*cmp)(b,c,v) > 0 ? b : (*cmp)(a,c,v) > 0 ? c : a);
}

 void
#ifdef KR_headers
qsortv(a, n, es, cmp, v) char *a; size_t n; size_t es; cmpfunc *cmp; char *v;
#else
#define a ((char*)a0)
qsortv(void *a0, size_t n, size_t es, cmpfunc *cmp, void *v)
#endif
{
	char *ae, *pa, *pb, *pc, *pd, *pl, *pm, *pn;
	int r, swaptype;
	size_t s;
	long t;

	SWAPINIT(a, es);
	if (n < 7) {	 /* Insertion sort on smallest arrays */
		for (pm = a + es, ae = a + n*es; pm < ae; pm += es)
			for (pl = pm; pl > a && (*cmp)(pl-es, pl, v) > 0; pl -= es)
				swap(pl, pl-es);
		return;
	}
	pm = a + (n/2)*es;    /* Small arrays, middle element */
	if (n > 7) {
		pl = a;
		pn = a + (n-1)*es;
		if (n > 40) {    /* Big arrays, pseudomedian of 9 */
			s = (n/8)*es;
			pl = med3(pl, pl+s, pl+2*s, cmp, v);
			pm = med3(pm-s, pm, pm+s, cmp, v);
			pn = med3(pn-2*s, pn-s, pn, cmp, v);
		}
		pm = med3(pl, pm, pn, cmp, v); /* Mid-size, med of 3 */
	}
	pa = pb = a;
	pc = pd = a + (n-1)*es;
	for (;;) {
		for ( ; pb <= pc; pb += es) {
			if (pb != pm)
				if ((r = (*cmp)(pb, pm, v)) > 0) break;
				else if (r < 0) continue;
			swap(pa, pb); pm = pa; pa += es;
		}
		for ( ; pb <= pc; pc -= es) {
			if (pc != pm)
				if ((r = (*cmp)(pc, pm, v)) < 0) break;
				else if (r > 0) continue;
			swap(pc, pd); pm = pd; pd -= es;
			}
		if (pb > pc) break;
		swap(pb, pc);
		pb += es;
		pc -= es;
	}
	pn = a + n*es;
	s = min(pa-a,  pb-pa   ); vecswap(a,  pb-s, s);
	s = min(pd-pc, pn-pd-es); vecswap(pb, pn-s, s);
	if ((s = pb-pa) > es) qsortv(a,    s/es, es, cmp, v);
	if ((s = pd-pc) > es) qsortv(pn-s, s/es, es, cmp, v);
}
