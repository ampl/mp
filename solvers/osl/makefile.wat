# /****************************************************************
# Copyright (C) 1997, 1999 Lucent Technologies
# All Rights Reserved
#
# Permission to use, copy, modify, and distribute this software and
# its documentation for any purpose and without fee is hereby
# granted, provided that the above copyright notice appear in all
# copies and that both that the copyright notice and this
# permission notice and warranty disclaimer appear in supporting
# documentation, and that the name of Lucent or any of its entities
# not be used in advertising or publicity pertaining to
# distribution of the software without specific, written prior
# permission.
#
# LUCENT DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
# INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.
# IN NO EVENT SHALL LUCENT OR ANY OF ITS ENTITIES BE LIABLE FOR ANY
# SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER
# IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,
# ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
# THIS SOFTWARE.
# ****************************************************************/


# For making osl.exe with WATCOM C/C++ 11.x under NT on Intel
# machines.

# Invoke with "wmake -u -f makefile.wat" .

CC = wcc386
S = ..
# "S = .." assumes this directory is solvers/osl .
CFLAGS = -fpd -bt=nt -DOSLMSDLL -DOSLCALLBACK=__cdecl
M = osllib\program\oslmd202.lib
# Adjust "M = ..." appropriately.
L = $S\amplsolv.lib $M

.c.obj:
	$(CC) $(CFLAGS) -I$S $*.c

o = osl.obj $S\funcadd0.obj

osl.exe: $o
	wcl386 -l=nt -bt=nt -fe=osl.exe $o $L

osl.obj: osl.c osl.h $S/asl.h
