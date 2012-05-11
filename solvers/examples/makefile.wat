# /****************************************************************
# Copyright (C) 1997-1999 Lucent Technologies
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

# For compiling solvers/examples with WATCOM C/C++.
# In the CFLAGS line, -bt=nt is for NT console binaries; for MSDOS,
# adjust -bt suitably.

# Invoke with "wmake -u -f makefile.wat" ,
# or "copy makefile.wat makefile", possibly edit makefile,
# and just invoke with "wmake -u".

.SUFFIXES: .f

CC = wcc386
S = ..
# "S = .." assumes this directory is solvers/examples .

CFLAGS = -fpd -I$S -bt=nt
LFLAGS = -bt=nt -l=nt

.c.obj:
	$(CC) $(CFLAGS) $*.c

.f.obj:
	f2cx -A $*.f
	$(CC) $(CFLAGS) $*.c
	del $*.c

Lf2c = $S\watf2c.lib
Lp = cport.lib $(Lf2c)
# cport.lib = PORT subroutine library.
# Ask netlib@netlib.org to
#	send dmngb dmnhb dn2g dn2gb dq7rgs from port
# to get source for the relevant (small) part of PORT.
# f2c.lib is for f2c's version of Fortran I/O: we assume
# cport.lib was compiled by f2c and $CC.

A = $S\amplsolv.lib
Af = $S\funcadd0.obj $A
all1 = lin1.exe lin2.exe lin3.exe linrc.exe dualconv.exe qtest.exe
all = mng1.exe fmng1.exe nl21.exe fnl21.exe mng.exe mnh.exe nl2.exe tn.exe v8.exe ve08.exe

all: $(all1) $(all)
	echo all made >all

# "make all" requires various routines from netlib's opt and port
# directories; "make all1" requires just files distributed with
# this directory.

all1: $(all1)

lin1 = lin1.obj $(Af)
lin1.exe: $(lin1)
	wcl386 $(LFLAGS) -fe=lin1 $(lin1)

lin2 = lin2.obj $A
lin2.exe: $(lin2)
	wcl386 $(LFLAGS) -fe=lin2 $(lin2)

lin3 = lin3.obj $A
lin3.exe: $(lin3)
	wcl386 $(LFLAGS) -fe=lin3 $(lin3)

linrc = linrc.obj $A
linrc.exe: $(linrc)
	wcl386 $(LFLAGS) -fe=linrc $(linrc)

mng1 = mng1.obj $(Af)

mng1.exe: $(mng1)
	wcl386 $(LFLAGS) -fe=mng1 $(mng1) $(Lp)

nl21 = nl21.obj $(Af)

nl21.exe: $(nl21)
	wcl386 $(LFLAGS) -fe=nl21 $(nl21) $(Lp)

mng = mng.objbj rvmsg.obj keywds.obj $A

mng.exe: $(mng)
	wcl386 $(LFLAGS) -fe=mng $(mng) $(Lp)

mnh = mnh.objbj rvmsg.obj hkeywds.obj $A

mnh.exe: $(mnh)
	wcl386 $(LFLAGS) -fe=mnh $(mnh) $(Lp)

nl2 = nl2.objbj rvmsg.obj keywds.obj $A

nl2.exe: $(nl2)
	wcl386 $(LFLAGS) -fe=nl2 $(nl2) $(Lp)

fmng1 = fmng1.obj $(Af)

fmng1.exe: $(fmng1)
	wcl386 $(LFLAGS) -fe=fmng1 $(fmng1) $(Lp)

fnl21 = fnl21.obj $(Af)

fnl21.exe: $(fnl21)
	wcl386 $(LFLAGS) -fe=fnl21 $(fnl21) $(Lp)

# tn.f and blas.f are available from netlib.exe: ask for "tn from opt".

tn = tnmain.obj tn.obj blas.obj $(Af)

tn.exe: $(tn)
	wcl386 $(LFLAGS) -fe=tn $(tn) $(Lf2c)

# ve08ad.f is available from netlib.exe: ask for "ve08 from opt".
# See comments in README about a possible adjustment to ve08ad.f .

v8.o.exe: ve08.c

v8 = v8.obj ve08ad.obj $(Af)

v8.exe: $(v8)
	wcl386 $(LFLAGS) -fe=v8 $(v8) $(Lp)

ve08 = ve08.obj ve08ad.obj $(Af)

ve08.exe: $(ve08)
	wcl386 $(LFLAGS) -fe=ve08 $(ve08) $(Lp)

dualconv = dualconv.obj $A
dualconv.exe: $(dualconv)
	wcl386 $(LFLAGS) -fe=dualconv $(dualconv)

nlcopy = nlcopy.obj $A
nlcopy.exe: $(nlcopy)
	wcl386 $(LFLAGS) -fe=nlcopy $(nlcopy)

qtest.exe: qtest.obj
	wcl386 $(LFLAGS) -fe=qtest qtest.obj $(Af)

$S\amplsolv.lib.exe:
	cd $S; make amplsolv.lib

clean:
	deltree /Y *.exe *.obj all test test1 dietd.duw

### "make test" to compare outputs.

### "make test1" to compare outputs that do not require things
### from other netlib directories.

test1:
	lin1 diet | rmcr >foo
	diff foo lin1diet.out
	lin2 diet | rmcr >foo
	diff foo lin2diet.out
	lin3 diet | rmcr >foo
	diff foo lin3diet.out
	linrc diet | rmcr >foo
	diff foo linrcdie.out
	dualconv diet foo
	rmcr foo.nl >goo
	diff goo dietd.nl
	copy foo.duw dietd.duw
	dualconv -u dietd foo
	rmcr foo.sol >goo
	diff goo dietu.sol
	del foo.nl
	del foo.sol
	qtest ch3 | rmcr >foo
	diff foo ch3qtest.out
	echo test1 done >test1

test: all test1
	mng ch3 | rmcr >foo
	diff foo ch3mng.out
	mnh ch3 | rmcr >foo
	diff foo ch3mnh.out
	nl2 ch3 | rmcr >foo
	diff foo ch3nl2.out
	nl21 ch3 | rmcr >foo
	diff foo ch3nl21.out
	tn ch3 | rmcr >foo
	diff foo ch3tn.out
	mng ch3 -AMPL
	rmcr ch3.sol >foo
	diff foo ch3mng.sol
	mnh ch3 -AMPL
	rmcr ch3.sol >foo
	diff foo ch3mnh.sol
	nl2 ch3 -AMPL
	rmcr ch3.sol >foo
	diff foo ch3nl2.sol
	tn ch3 -AMPL
	rmcr ch3.sol >foo
	diff foo ch3tn.sol
	del ch3.sol
	v8 ch3 | rmcr >foo
	diff foo v8ch3.out
	ve08 ch3 | rmcr >foo
	diff foo ve08ch3.out
	del foo
	echo tests done. >test
