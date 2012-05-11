/****************************************************************
Copyright (C) 1994 Lucent Technologies
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

#ifndef sgi
#ifdef sun
#define sprintf Sprintf
#else
#define ftncls_  ftncls
#define ftnopn_  ftnopn
#define ekkbasi_ ekkbasi
#define ekkbaso_ ekkbaso
#define ekkbcdo_ ekkbcdo
#define ekkbslv_ ekkbslv
#define ekkcrsh_ ekkcrsh
#define ekkdsca_ ekkdsca
#define ekkdscm_ ekkdscm
#define ekkiget_ ekkiget
#define ekkimdl_ ekkimdl
#define ekkiset_ ekkiset
#define ekklmdl_ ekklmdl
#define ekkmpre_ ekkmpre
#define ekkmset_ ekkmset
#define ekkmsgu_ ekkmsgu
#define ekkmslv_ ekkmslv
#define ekkname_ ekkname
#define ekknslv_ ekknslv
#define ekknget_ ekknget
#define ekkprsl_ ekkprsl
#define ekkpssl_ ekkpssl
#define ekkqmdl_ ekkqmdl
#define ekkqslv_ ekkqslv
#define ekkrget_ ekkrget
#define ekkrset_ ekkrset
#define ekkscal_ ekkscal
#define ekksslv_ ekksslv
#ifdef WATCOM
#pragma aux ftncls "FTNCLS"
#pragma aux ftnopn "FTNOPN"
#pragma aux ekkbasi "EKKBASI"
#pragma aux ekkbaso "EKKBASO"
#pragma aux ekkbcdo "EKKBCDO"
#pragma aux ekkbslv "EKKBSLV"
#pragma aux ekkcrsh "EKKCRSH"
#pragma aux ekkdsca "EKKDSCA"
#pragma aux ekkdscm "EKKDSCM"
#pragma aux ekkiget "EKKIGET"
#pragma aux ekkimdl "EKKIMDL"
#pragma aux ekkiset "EKKISET"
#pragma aux ekklmdl "EKKLMDL"
#pragma aux ekkmpre "EKKMPRE"
#pragma aux ekkmset "EKKMSET"
#pragma aux ekkmsgu "EKKMSGU"
#pragma aux ekkmslv "EKKMSLV"
#pragma aux ekkname "EKKNAME"
#pragma aux ekknslv "EKKNSLV"
#pragma aux ekknget "EKKNGET"
#pragma aux ekkprsl "EKKPRSL"
#pragma aux ekkpssl "EKKPSSL"
#pragma aux ekkqmdl "EKKQMDL"
#pragma aux ekkqslv "EKKQSLV"
#pragma aux ekkrget "EKKRGET"
#pragma aux ekkrset "EKKRSET"
#pragma aux ekkscal "EKKSCAL"
#pragma aux ekksslv "EKKSSLV"
#endif /*WATCOM*/
#endif
#endif
