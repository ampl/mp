/****************************************************************
Copyright (C) 1994, 1999 Lucent Technologies
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

#ifdef OSLMSDLL
#define OSLLINKAGE __stdcall
#define OSLLIBAPI __declspec(dllimport)
#define WIN32
#else
#ifdef OSLMSSTAT
#define OSLLINKAGE __stdcall
#define WIN32
#endif
#endif

#ifdef _AIX
#define No_s
#endif

#ifdef WIN32
#define USE_UC
#endif

#ifdef sun
#define sprintf Sprintf
#endif

#ifdef USE_UC
#define ekkbasi_ EKKBASI
#define ekkbaso_ EKKBASO
#define ekkbcdo_ EKKBCDO
#define ekkbmpr_ EKKBMPR
#define ekkbmps_ EKKBMPS
#define ekkbslv_ EKKBSLV
#define ekkclcb  EKKCLCB
#define ekkcrsh_ EKKCRSH
#define ekkdsca_ EKKDSCA
#define ekkdscm_ EKKDSCM
#define ekkfcls_ EKKFCLS
#define ekkfopn_ EKKFOPN
#define ekkiget_ EKKIGET
#define ekkimdl_ EKKIMDL
#define ekkiset_ EKKISET
#define ekklmdl_ EKKLMDL
#define ekkinit_ EKKINIT
#define ekkmpre_ EKKMPRE
#define ekkmset_ EKKMSET
#define ekkmslv_ EKKMSLV
#define ekkname_ EKKNAME
#define ekknfes_ EKKNFES
#define ekknslv_ EKKNSLV
#define ekknget_ EKKNGET
#define ekknwmt_ EKKNWMT
#define ekkprsl_ EKKPRSL
#define ekkpssl_ EKKPSSL
#define ekkqmdl_ EKKQMDL
#define ekkqpar_ EKKQPAR
#define ekkqslv_ EKKQSLV
#define ekkrget_ EKKRGET
#define ekkrgcb  EKKRGCB
#define ekkrset_ EKKRSET
#define ekksbnd_ EKKSBND
#define ekkscal_ EKKSCAL
#define ekksobj_ EKKSOBJ
#define ekksos_  EKKSOS
#define ekkspar_ EKKSPAR
#define ekksslv_ EKKSSLV
#endif

#ifdef No_s
#define ekkbasi_ ekkbasi
#define ekkbaso_ ekkbaso
#define ekkbcdo_ ekkbcdo
#define ekkbmpr_ ekkbmpr
#define ekkbmps_ ekkbmps
#define ekkbslv_ ekkbslv
#define ekkcrsh_ ekkcrsh
#define ekkdsca_ ekkdsca
#define ekkdscm_ ekkdscm
#define ekkfcls_ ekkfcls
#define ekkfopn_ ekkfopn
#define ekkiget_ ekkiget
#define ekkimdl_ ekkimdl
#define ekkiset_ ekkiset
#define ekklmdl_ ekklmdl
#define ekkinit_ ekkinit
#define ekkmpre_ ekkmpre
#define ekkmset_ ekkmset
#define ekkmslv_ ekkmslv
#define ekkname_ ekkname
#define ekknfes_ ekknfes
#define ekknslv_ ekknslv
#define ekknget_ ekknget
#define ekknwmt_ ekknwmt
#define ekkprsl_ ekkprsl
#define ekkpssl_ ekkpssl
#define ekkqmdl_ ekkqmdl
#define ekkqpar_ ekkqpar
#define ekkqslv_ ekkqslv
#define ekkrget_ ekkrget
#define ekkrset_ ekkrset
#define ekksbnd_ ekksbnd
#define ekkscal_ ekkscal
#define ekksobj_ ekksobj
#define ekksos_  ekksos
#define ekkspar_ ekkspar
#define ekksslv_ ekksslv
#endif /* No_s */
