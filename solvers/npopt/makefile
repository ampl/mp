# /****************************************************************
# Copyright (C) 1997-1998 Lucent Technologies
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
S = ..
# "S = .." assumes this directory is solvers/npopt .

# Add -DKR_headers to CFLAGS if your C compiler does not
# understand ANSI C function headers, e.g.
#	CFLAGS = -O -DKR_headers
# If things don't run right, you may need to change -O to -g
# so you can poke around with a debugger.
# You may also need to add -D_SIZE_T to CFLAGS, or to
# add -D_SIZE_T to the CFLAGS assignment above (if you get
# complaints about size_t already being defined).  This is
# necessary on at least some Sun systems.

# For Sun systems, try
# CFLAGS = -O -DKR_headers -D_SIZE_T

.c.o:
	$(CC) -c $(CFLAGS) -I$S $*.c

#### We assume the f2c calling conventions:

.f.o:
	f2c -c -g $*.f
	$(CC) -c $(CFLAGS) $*.c
	rm $*.c

n = npopt.o version.o $S/amplsolver.a npopt.a

npopt: $n
	$(CC) -o npopt $n -lf2c -lm

npopt.o: $S/asl.h $S/getstub.h

$S/amplsolver.a:
	cd ..;make amplsolver.a

# make xsum.out to check for transmission errors.
# This assumes you have the xsum program, whose source
# you can get by asking research!netlib to
#	send xsum.c from f2c/src

xs0 = README makefile npopt.c version.c

xsum.out: $(xs0)
	xsum $(xs0) >xsum1.out
	cmp xsum0.out xsum1.out && mv xsum1.out xsum.out || diff xsum[01].out
