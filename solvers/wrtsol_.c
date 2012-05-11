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

#include "getstub.h"

 static int
#ifdef KR_headers
slen(s, len) register char *s; ftnlen len;
#else
slen(register char *s, ftnlen len)
#endif
{
	register char *s1 = s + len;
	while(s1 > s)
		if (*--s1 > ' ') {
			s1++;
			break;
			}
	return s1 - s;
	}

 void
#ifdef KR_headers
wrsolw_(msg, nmsg, x, y, ws, msg_len) char *msg; fint *nmsg, *ws;
				  real *x, *y; ftnlen msg_len;
#else
wrsolw_(char *msg, fint *nmsg, real *x, real *y, fint *ws, ftnlen msg_len)
#endif
{
	char *b, *buf, *me;
	int i, len, nm;
	ASL *asl;
	Option_Info oi;

	if (!(asl = cur_ASL))
		badasl_ASL(asl,0,"wrtsol");
	nm = (int)*nmsg;
	len = nm + 1;
	me = msg + nm*msg_len;
	for(b = msg; b < me; b += msg_len)
		len += slen(b,msg_len);
	b = buf = (char *)Malloc(len);
	if (nm)
		for(;;) {
			if (i = slen(msg, msg_len)) {
				strncpy(b, msg, i);
				b += i;
				}
			msg += msg_len;
			if (msg >= me)
				break;
			*b++ = '\n';
			}
	*b = 0;
	if (!(oi.wantsol = *ws & 7))
		oi.wantsol = 1;
	write_sol_ASL(asl, buf, x, y, &oi);
	free(buf);
	}

 void
#ifdef KR_headers
wrtsol_(msg, nmsg, x, y, msg_len) char *msg; fint *nmsg;
				  real *x, *y; ftnlen msg_len;
#else
wrtsol_(char *msg, fint *nmsg, real *x, real *y, ftnlen msg_len)
#endif
{
	static fint ws = 7;
	wrsolw_(msg, nmsg, x, y, &ws, msg_len);
	}
