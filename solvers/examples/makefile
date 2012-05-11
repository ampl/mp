# /****************************************************************
# Copyright (C) 1997-1999, 2001 Lucent Technologies
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

# For compiling solvers/examples on Unix and Linux systems.

# Invoke with "make -f makefile.u",
# or "cp makefile.u makefile", edit makefile as necessary,
# and just invoke with "make".

.SUFFIXES: .c .f .o
CC = cc
S = ..
# "S = .." assumes this directory is solvers/examples .

CFLAGS = -g
# Add -DKR_headers to CFLAGS if your C compiler does not
# understand ANSI C function headers, e.g.
#	CFLAGS = -O -DKR_headers
# If things don't run right, you may need to change -O to -g
# so you can poke around with a debugger.

# See $S/makefile for other comments about CFLAGS.

.c.o:
	$(CC) -c $(CFLAGS) -I$S $*.c

## On HP systems, add -ldld to the end of the L = line below.
## On linux and Solaris systems, add -ldl to the end of the L = line below.

L = -lm
Lf2c = -lf2c
cport = cport.a
# $(cport) = PORT subroutine library.
# Ask netlib@netlib.org to
#	send dmngb dmnhb dn2g dn2gb dq7rgs from port
# to get source for the relevant (small) part of PORT.
# -lf2c is for f2c's version of Fortran I/O: we assume
# $(cport) was compiled by f2c and cc.
Lp = $(cport) $(Lf2c)

FC = f77
FFLAGS =
# Fortran compiler (fc script provided with f2c)

.f.o:
	$(FC) -c $(FFLAGS) $*.f

A = $S/amplsolver.a
Af = $S/funcadd0.o $A

all1 = evalchk lin1 lin2 lin3 linrc dualconv nlcopy qtest
all = mng1 fmng1 nl21 fnl21 mng mnh nl2 tn v8 ve08

all: $(all1) $(all)

# "make all" requires various routines from netlib's opt and port
# directories; "make all1" requires just files distributed with
# this directory.

all1: $(all1)

evalchk = evalchk.o $A
evalchk: $(evalchk)
	$(CC) - o evalchk $(evalchk)

lin1 = lin1.o $(Af)
lin1: $(lin1)
	$(CC) -o lin1 $(lin1)

lin2 = lin2.o $(Af)
lin2: $(lin2)
	$(CC) -o lin2 $(lin2)

lin3 = lin3.o $(Af)
lin3: $(lin3)
	$(CC) -o lin3 $(lin3)

linrc = linrc.o $(Af)
linrc: $(linrc)
	$(CC) -o linrc $(linrc)

mng1 = mng1.o $A

mng1: $(mng1)
	$(CC) -o mng1 $(mng1) $(Lp) $L

nl21 = nl21.o $A

nl21: $(nl21)
	$(CC) -o nl21 $(nl21) $(Lp) $L

mng = mng.o rvmsg.o keywds.o $A

mng: $(mng)
	$(CC) -o mng $(mng) $(Lp) $L

mnh = mnh.o rvmsg.o hkeywds.o $A

mnh: $(mnh)
	$(CC) -o mnh $(mnh) $(Lp) $L

nl2 = nl2.o rvmsg.o keywds.o $A

nl2: $(nl2)
	$(CC) -o nl2 $(nl2) $(Lp) $L

fmng1 = fmng1.o $A

fmng1: $(fmng1)
	$(CC) -o fmng1 $(fmng1) $(Lp) $L

fnl21 = fnl21.o $A

fnl21: $(fnl21)
	$(CC) -o fnl21 $(fnl21) $(Lp) $L

# tn.f and blas.f are available from netlib: ask for "tn from opt".

tn = tnmain.o tn.o blas.o $A

tn: $(tn)
	$(CC) -o tn $(tn) $(Lf2c) $L

# ve08ad.f is available from netlib: ask for "ve08 from opt".
# See comments in README about a possible adjustment to ve08ad.f .

v8.o: ve08.c

v8 = v8.o ve08ad.o $A

v8: $(v8)
	$(CC) -o v8 $(v8) $(Lp) $L

ve08 = ve08.o ve08ad.o $A

ve08: $(ve08)
	$(CC) -o ve08 $(ve08) $(Lf2c) $L

dualconv = dualconv.o $(Af)
dualconv: $(dualconv)
	$(CC) -o dualconv $(dualconv)

nlcopy = nlcopy.o $A
nlcopy: $(nlcopy)
	$(CC) -o nlcopy $(nlcopy) $L

qtest: qtest.o $(Af)
	$(CC) -o qtest qtest.o $(Af) $L

$S/amplsolver.a:
	cd $S; make amplsolver.a

clean:
	rm -rf $(all1) $(all) *.o dietd.duw *.mex* t bad bad1

# make xsum.out to check for transmission errors.
# This assumes you have the xsum program, whose source
# you can get by asking research!netlib to
#	send xsum.c from f2c/src

xs0 = README README.gjh amplfun4.c amplfunc.c dminos dualconv.c\
 enewt.m evalchk.c evalf.m evalg.m evalw.m fmng1.f fnl21.f gjh.c\
 hist.c hkeywds.c init.m keywds.c lin1.c lin2.c\
 lin3.c linrc.c makefile.sy makefile.u makefile.vc makefile.wat\
 mng.c mng1.c mnh.c nl2.c nl21.c nlcopy.c qtest.c rvmsg.c rvmsg.h\
 savesol.m spamfunc.c tnmain.c v8.c ve08.c\
 diet.amp diet.nl diet.row diet.col dietl1.out dietl2.out dietl3.out\
 dietd.nl diet.sol dietd.sol dietu.sol ch3.amp ch3.nl\
 ch3mng.out ch3mng.sol ch3mnh.out ch3mnh.sol ch3nl21.out ch3nl2.out\
 ch3nl2.sol ch3qtest.out ch3tn.out ch3tn.sol hs100.amp hs100.nl\
 hs100.so0 lin1diet.out lin2diet.out lin3diet.out linrcdie.out\
 v8ch3.out ve08ch3.out


xsum.out: $(xs0)
	xsum $(xs0) >xsum1.out
	cmp xsum0.out xsum1.out && mv xsum1.out xsum.out || diff xsum[01].out

### "make test" to compare outputs; equivalent to "make test1 test2".

### "make test1" to compare outputs that do not require things
### from other netlib directories.

### Anything left in directory "bad" or "bad1" may indicate something
### amiss.  Trivial differences in function and gradient values in the
### *.out files (particularly ve08ch3.out) may be due to differences in
### special-function implementations among math libraries.  Other
### differences are cause for concern.

test: test1 test2

test1: all1
	mkdir bad1
	lin1 diet >foo
	cmp foo lin1diet.out || mv foo bad1/lin1diet.out
	lin2 diet >foo
	cmp foo lin2diet.out || mv foo bad1/lin2diet.out
	lin3 diet >foo
	cmp foo lin3diet.out || mv foo bad1/lin3diet.out
	linrc diet >foo
	cmp foo linrcdie.out || mv foo bad1/linrcdie.out
	dualconv diet foo
	cmp foo.nl dietd.nl || mv foo.nl bad1/dietd.nl
	mv foo.duw dietd.duw
	dualconv -u dietd foo
	cmp foo.sol dietu.sol || mv foo.sol bad1/dietu.sol
	rm -f foo.nl foo.sol
	qtest ch3 >foo
	cmp foo ch3qtest.out || mv foo bad1/ch3qtest.out
	mkdir t
	nlcopy -g hs100 t/hs100
	cmp hs100.nl t/hs100.nl
	rm -rf t
	rmdir bad1

test2: all
	mkdir bad
	mng ch3 >foo 2>&1
	cmp foo ch3mng.out || mv foo bad/ch3mng.out
	mnh ch3 >foo 2>&1
	cmp foo ch3mnh.out || mv foo bad/ch3mnh.out
	nl2 ch3 >foo 2>&1
	cmp foo ch3nl2.out || mv foo bad/ch3nl2.out
	nl21 ch3 >foo 2>&1
	cmp foo ch3nl21.out || mv foo bad/ch3nl21.out
	tn ch3 >foo 2>&1
	cmp foo ch3tn.out || mv foo bad/ch3tn.out
	mng ch3 -AMPL
	cmp ch3.sol ch3mng.sol || mv ch3.sol bad/ch3mng.sol
	mnh ch3 -AMPL
	cmp ch3.sol ch3mnh.sol || mv ch3.sol bad/ch3mnh.sol
	nl2 ch3 -AMPL
	cmp ch3.sol ch3nl2.sol || mv ch3.sol bad/ch3nl2.sol
	tn ch3 -AMPL
	cmp ch3.sol ch3tn.sol || mv ch3.sol bad/ch3tn.sol
	rm -f ch3.sol
	v8 ch3 >foo 2>&1
	cmp foo v8ch3.out || mv foo bad/v8ch3.out
	ve08 ch3 >foo 2>&1
	cmp foo ve08ch3.out || mv foo bad/ve08ch3.out
	rmdir bad
	echo tests done.

## For use with MATLAB 5.x or 6.x, "make amplfunc.mex" or
## "make spamfunc.mex".
## Make sure your PATH includes MATLAB's "bin" directory.

af = amplfunc.c $S/funcaddr.o $A

amplfunc.mex: $(af)
	mex -I$S $(af)
	touch amplfunc.mex

saf = spamfunc.c $S/funcaddr.o $A

spamfunc.mex: $(saf)
	mex -I$S $(saf)
	touch spamfunc.mex

$S/funcaddr.o: $S/funcadd1.c
	cd $S; make funcaddr.o

## Variants for MATLAB 4.x: "make amplfun4.mex" or "make spamfun4.mex".

af4 = amplfun4.c $A

amplfun4.mex: $(af4)
	cmex -I$S $(af4)
	touch amplfun4.mex

saf4 = spamfun4.c $A

spamfun4.mex: $(saf4)
	cmex -I$S $(saf4)
	touch spamfun4.mex

gjh: gjh.o
	$(CC) -o gjh gjh.o $A $L

# To make "hist", you must have the GNU history and readline packages installed.
# When using "hist" with "ampl", invoke "hist ampl -vi2 ..." to see prompts.
hist: hist.o
	$(CC) -o hist hist.o -lhistory -lreadline
