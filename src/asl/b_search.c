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

/* Binary search for keywords. */
/* Assuming
 *
 *	typedef struct { char *name; ... } option_word;
 *	option_word keywords[] = { ... };
 *
 * b_search has calling sequence
 *
 *	b_search((char *)keywords, (int)sizeof(option_word),
 *		(int)(number of keywords = sizeof(keywords)/sizeof(option_word),
 *		char **sp, char **peq)
 *
 * where *sp is a string that may contain a keyword.
 * No keyword ==> b_search returns 0 with **sp = 0.
 * Bad keyword ==> b_search returns 0 with *sp = the bad keyword.
 * Keyword matched ==> b_search returns the matched option_word*,
 * sets *sp to the next token (skipping '=' if present), and sets
 * *peq to "=" if '=' is present and to " " otherwise.
 */

#ifdef __cplusplus
extern "C" {
#endif


#ifdef Use_tolower
#include "ctype.h"
#define Tolower(x) tolower(x)
#else
static unsigned char lc[256];

 static void
lc_init(void)
{
	int i;
	const char *s;
	for(i = 0; i < 256; i++)
		lc[i] = i;
	for(s = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"; *s; s++)
		lc[(int)*s] = *s + 'a' - 'A';
	}
#define Tolower(x) lc[x]
#endif

 void *
b_search_ASL(void *ow, int owsize, int n, char **sp, char **peq)
{
	int c, c1, c2, n1;
	char *s, *s1, *s2;
	void *ow1;
	static char Blank[] = " ", Eq[] = "=";
#ifndef Use_tolower
	static int first = 1;
	if (first) {
		lc_init();
		first = 0;
		}
#endif

	for(s = *sp; (c = *(unsigned char *)s) <= ' '; s++)
		if (!c)
			goto no_ow1;

	/* binary search */

	while(n > 0) {
		ow1 = (char*)ow + (n1 = n >> 1)*owsize;
		s2 = *(char **)ow1;
		for(s1 = s;; s1++) {
			c1 = Tolower(*(unsigned char *)s1);
			if (!(c2 = *s2++)) {
				if (c1 <= ' ' || c1 == '=')
					goto found;
				break;
				}
			if (c1 != c2)
				break;
			}
		if (c1 == '=' || c1 < c2)
			n = n1;
		else {
			n -= n1 + 1;
			ow = (char*)ow1 + owsize;
			}
		continue;
 found:
		*peq = Blank;
		while(*s1 && *s1 <= ' ')
			s1++;
		if (*s1 == '=') {
			*peq = Eq;
			while(*++s1 && *s1 <= ' ');
			}
		*sp = s1;
		return ow1;
		}
 no_ow1:
	ow1 = 0;
	*sp = s;
	return 0;
	}

#ifdef __cplusplus
	}
#endif
