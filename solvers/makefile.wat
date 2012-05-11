# /****************************************************************
# Copyright (C) 1997-2001 Lucent Technologies
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

# For making amplsolv.lib with WATCOM C/C++.
# In the CFLAGS line, -bt=nt is for NT console binaries; for MSDOS,
# adjust -bt suitably, change "-DWIN32" to "-DNO_FUNCADD" in the
# funcadd1.obj rule below, and change fpinitmt.obj to fpinit.obj
# in the "a = ..." assignment.

# Invoke with "wmake -u -f makefile.wat" .

# NOTE: WATCOM C/C++ 11.0 mishandles some typedefs.  (This bug
# was not present in version 10.6 and is fixed in version 11.0a.)
# To overcome this bug, in nlp.h, change all occurrences of
# "efunc *op" to "real (*op)(expr*)", and in nlp2.h, change
# all occurrences of "efunc2 *op" to "real (*op)(expr2*)".

CC = wcc386
CFLAGS = -fpd -DWATCOM -bt=nt

.c.obj:
	$(CC) $(CFLAGS) $*.c

all: arith.h stdio1.h amplsolv.lib funcadd0.obj

a = \
	asldate.obj \
	atof.obj \
	auxinfo.obj \
	avldelete.obj \
	avltree.obj \
	b_search.obj \
	basename.obj \
	bscanf.obj \
	com2eval.obj \
	comeval.obj \
	con1ival.obj \
	con2ival.obj \
	con2val.obj \
	conadj.obj \
	conpval.obj \
	conscale.obj \
	conval.obj \
	derprop.obj \
	details.obj \
	dtoa1.obj \
	duthes.obj \
	dynlink.obj \
	f_read.obj \
	fg_read.obj \
	fg_write.obj \
	fgh_read.obj \
	fpecatch.obj \
	fpinitmt.obj \
	fullhes.obj \
	func_add.obj \
	funcadd1.obj \
	g_fmt.obj \
	genrowno.obj \
	getenv.obj \
	getstub.obj \
	htcl.obj \
	jac0dim.obj \
	jac2dim.obj \
	jacdim.obj \
	jacinc.obj \
	jacinc1.obj \
	mach.obj \
	mainexit.obj \
	mip_pri.obj \
	misc.obj \
	mypow.obj \
	names.obj \
	nl_obj.obj \
	nqpcheck.obj \
	obj2val.obj \
	obj_prec.obj \
	objconst.obj \
	objval.obj \
	objval_.obj \
	op_type.obj \
	pfg_read.obj \
	pfghread.obj \
	printf.obj \
	pshvprod.obj \
	punknown.obj \
	qp_read.obj \
	qpcheck.obj \
	qsortv.obj \
	readsol.obj \
	repwhere.obj \
	rops.obj \
	rops2.obj \
	sigcatch.obj \
	sos_add.obj \
	sphes.obj \
	sscanf.obj \
	stderr.obj \
	studchk0.obj \
	suf_sos.obj \
	value.obj \
	writesol.obj \
	wrtsol_.obj \
	ws_desc.obj \
	wsu_desc.obj \
	x2check.obj \
	xectim.obj \
	xp1known.obj \
	xp2known.obj

amplsolv.lib: $a
	wlib -c amplsolv.lib @amplsolv

Aslh = arith.h asl.h funcadd.h stdio1.h
auxinfo.obj: funcadd.h stdio1.h
mach.obj: arith.h
avldelete.obj avltree.obj bscanf.obj conscale.obj derprop.obj dynlink.obj func_add.obj\
 funcadd.obj funcadd1.obj funcaddk.obj funcaddr.obj funcadd0.obj g_fmt.obj\
 genrowno.obj jac0dim.obj jacdim.obj jac2dim.obj jacinc.obj jacinc1.obj names.obj\
 obj_prec.obj objval_.obj repwhere.obj sigcatch.obj sjac0dim.obj studchk0.obj:\
 $(Aslh)
avldelete.obj avltree.obj: avltree.h
xp1known.obj: asl_pfg.h psinfo.h nlp.h $(Aslh)
duthes.obj fullhes.obj htcl.obj sphes.obj: asl_pfgh.h psinfo.h nlp2.h $(Aslh)
getstub.obj value.obj writesol.obj wrtsol_.obj: getstub.h $(Aslh)
com2eval.obj con2ival.obj con2val.obj obj2val.obj\
 x2check.obj: jac2dim.h nlp2.h $(Aslh)
conpval.obj pshvprod.obj xp2known.obj:\
	jacpdim.h asl_pfgh.h psinfo.h nlp2.h $(Aslh)
comeval.obj con1ival.obj conval.obj mip_pri.obj objval.obj qpcheck.obj\
 readsol.obj: nlp.h $(Aslh)
misc.obj nl_obj.obj sos_add.obj suf_sos.obj:\
	nlp.h nlp2.h asl_pfg.h asl_pfgh.h psinfo.h $(Aslh)
op_type.obj: op_type.hd op_typeb.hd
fgh_read.obj: jac2dim.h opnos.hd op_type.hd dvalue.hd nlp2.h $(Aslh)
rops.obj: nlp.h errchk.h $(Aslh)
rops2.obj: nlp2.h errchk.h $(Aslh)
conadj.obj fg_write.obj qp_read.obj: nlp.h r_opn.hd $(Aslh)
f_read.obj fg_read.obj: nlp.h r_opn.hd dvalue.hd $(Aslh)
objconst.obj: r_opn0.hd nlp.h nlp2.h asl_pfg.h asl_pfgh.h psinfo.h $(Aslh)
pfg_read.obj: asl_pfg.h r_opn0.hd dvalue.hd nlp.h psinfo.h $(Aslh)
pfghread.obj: jacpdim.h asl_pfgh.h opnos.hd r_opn0.hd dvalue.hd\
	psinfo.h nlp2.h $(Aslh)
nqpcheck.obj: nlp.h r_qp.hd $(Aslh)
printf.obj punknown.obj sscanf.obj: stdio1.h
dtoa1.obj: dtoa.c arith.h stdio1.h

# Use CFLAGS in compiling arithchk.c in case something in CFLAGS affects
# the number of bits in integral data types.  (It's probably best not to
# add such options to CFLAGS.)

arith.h: arithchk.c
	comptry.bat wcl386 -DNO_FPINIT arithchk.c
	arithchk >arith.h
	del arithchk.exe
	del arithchk.obj

stdio1.h: stdio1.h0
	copy stdio1.h0 stdio1.h

funcadd1.obj: funcadd1.c
	$(CC) $(CFLAGS) -DWIN32 -DFUNCADD="funcadd_ASL_" funcadd1.c

stderr.obj: stderr.c
	$(CC) $(CFLAGS) -DSTDERR=stdout stderr.c

xectim.obj: xectim.c
	$(CC) $(CFLAGS) -DNO_RUSAGE xectim.c

details.c: details.c0
	echo create details.c by suitably editing details.c0
