/****************************************************************
Copyright (C) 1990 Lucent Technologies
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


.text
	entry rnd_prod
rnd_prod:
	using	rnd_prod,15
	ld	0,0+64(13)
	mxd	0,8+64(13)
	lrdr	0,0
	b	2(,14)
	drop
	entry	rnd_quot
rnd_quot:
	using	rnd_quot,15
	ld	0,0+64(13)
	ldr	2,0
	ld	4,8+64(13)
	ddr	2,4
	std	2,32(13)
	mxdr	4,2
	sdr	2,2
	sxr	0,4
	dd	0,8+64(13)
	sdr	2,2
	ld	4,32(13)
	sdr	6,6
	axr	0,4
	lrdr	0,0
	b	2(,14)
