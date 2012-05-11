# /****************************************************************
# Copyright (C) 1997 Lucent Technologies
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

.SUFFIXES: .c .f .o
CC = cc
CFLAGS = -O
F77 = f77
FFLAGS =
S = ..
# "S = .." assumes this directory is solvers/osl .

# Add -DKR_headers to CFLAGS if your C compiler does not
# understand ANSI C function headers, e.g.
#	CFLAGS = -O -DKR_headers
# If things don't run right, you may need to change -O to -g
# so you can poke around with a debugger.
# See $S/makefile for other comments about CFLAGS.

# For RS6000 systems, do not use -O unless you have cc version
# at least 01.02.0000.0009 -- earlier versions are too buggy.

.c.o:
	$(CC) -c $(CFLAGS) -I$S $*.c

#### We assume the f2c calling conventions:

.f.o:
	$(F77) -c $(FFLAGS) $*.f

o = osl.o ftncls.o ftnopn.o version.o $S/funcadd0.o $S/amplsolver.a

osl: $o
	$(F77) -o osl $o -losl

# The above works on an RS6000 that has a Fortran compiler.
# For an RS6000 with just the Fortran library libxlf.a,
#       $(CC) -bloadmap:zork -o osl $o -losl -bnso \
#	-bI:/lib/syscalls.exp libxlf.a -lm
# For Sun 4,
#	$(F77) -dalign -cg89 $o -losl1 -losl2 -L/usr/lpp/osl/lib -o osl

osl.o: $S/asl.h

$S/amplsolver.a:
	cd $S; make amplsolver.a

# make xsum.out to check for transmission errors.
# This assumes you have the xsum program, whose source
# you can get by asking research!netlib to
#	send xsum.c from f2c/src

xs0 = README.1st README.osl ftncls.f ftnopn.f makefile\
 osl.c osl.h version.c

xsum.out: $(xs0)
	xsum $(xs0) >xsum1.out
	cmp xsum0.out xsum1.out && mv xsum1.out xsum.out || diff xsum[01].out
