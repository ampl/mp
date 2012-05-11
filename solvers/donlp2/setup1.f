* /****************************************************************
* Copyright (C) 1997-9 Lucent Technologies
* All Rights Reserved
*
* Permission to use, copy, modify, and distribute this software and
* its documentation for any purpose and without fee is hereby
* granted, provided that the above copyright notice appear in all
* copies and that both that the copyright notice and this
* permission notice and warranty disclaimer appear in supporting
* documentation, and that the name of Lucent or any of its entities
* not be used in advertising or publicity pertaining to
* distribution of the software without specific, written prior
* permission.
*
* LUCENT DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
* INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.
* IN NO EVENT SHALL LUCENT OR ANY OF ITS ENTITIES BE LIABLE FOR ANY
* SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
* WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER
* IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,
* ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
* THIS SOFTWARE.
* ****************************************************************/

	subroutine setup1(b, da, L, name0, x0a)
	double precision da(2), x0a(*)
	integer b(*), L(15)
	character*(*) name0
	include 'O8COMM.INC'
	include 'O8FINT.INC'
	integer i, j, k, lineq, linge, nsbeq, nsbge

	name = name0

	del0 = da(1)
	tau0 = da(2)
	analyt = .true.
	cold = .true.
	meu = L(1)
	n = L(2)
	ng = L(3)
	nh = L(4)
	nreset = L(5)
	prou = L(6)
	silent = .false.
	if (L(7) .gt. 0) silent = .true.
	intakt = .false.
	if (L(15) .gt. 0) intakt = .true.
	te0 = .false.
	if (mod(L(13),2) .eq. 1) te0 = .true.
	te1 = .false.
	if (mod(L(13),4) .ge. 2) te1 = .true.
	te2 = .false.
	if (mod(L(13),8) .ge. 4) te2 = .true.
	te3 = .false.
	if (mod(L(13),16) .ge. 8) te3 = .true.
	do 10 i = 1, n
 10		x(i) = x0a(i)
	lineq = L(8)
	linge = L(9)
	nsbeq = L(10)
	nsbge = L(11)
	gunit(1,0) = -1
	gunit(2,0) = 0
	gunit(3,0) = 0
	gconst(0) = .false.
	if (L(12) .ne. 0) gconst(0) = .true.
	do 20 i = 1, nh-lineq
		gunit(1,i) = -1
		gunit(2,i) = 0
		gunit(3,i) = 0
		gconst(i) = .false.
 20		continue
	do 30 i = nh-lineq+1, nh-nsbeq
		gunit(1,i) = -1
		gunit(2,i) = 0
		gunit(3,i) = 0
		gconst(i) = .true.
 30		continue
	do 40 i = nh-nsbeq+1, nh
		gunit(1,i) = 1
		gunit(2,i) = b(i) + 1
		gunit(3,i) = 1
		gconst(i) = .true.
 40		continue
	j = nh + ng - linge
	do 50 i = nh+1, j
		gunit(1,i) = -1
		gunit(2,i) = 0
		gunit(3,i) = 0
		gconst(i) = .false.
 50		continue
	k = ng + nh - nsbge
	do 60 i = j+1, k
		gunit(1,i) = -1
		gunit(2,i) = 0
		gunit(3,i) = 0
		gconst(i) = .true.
 60		continue
	do 70 i = k+1, ng+nh
		gunit(1,i) = 1
		if (b(i) .lt. 0) then
			gunit(2,i) = -b(i)
			gunit(3,i) = -1
		else
			gunit(2,i) = b(i) + 1
			gunit(3,i) = 1
			endif
		gconst(i) = .true.
 70		continue
	call setup3(icf, icgf, cres(nh+1), cgres(nh+1), cres, cgres,xsc)
	end

	subroutine endinf
	include 'O8COMM.INC'
	double precision dinfo(6), umin
	integer j, k(2)
	umin = 0
	do j = nh+1, nres
		umin = min(umin,u(j))
		enddo
	k(1) = int(optite) + 11
	k(2) = itstep
	dinfo(1) = fx
	dinfo(2) = scf
	dinfo(3) = gfn
	dinfo(4) = b2n
	dinfo(5) = upsi
	dinfo(6) = umin
	call outinf(k, dinfo, x, u)
	end

	subroutine setup2(itmax, e, numsm1, epsph1)
	integer itmax, numsm1
	double precision e, epsph1
	include 'O8COMM.INC'

	if (intakt .and. silent) then
		silent = .false.
		open(meu,file='/dev/null',status='UNKNOWN',err=10)
		open(prou,file='/dev/null',status='UNKNOWN')
		go to 20
 10		open(meu,file='NUL',status='UNKNOWN',err=20)
		open(prou,file='NUL',status='UNKNOWN')
		endif
 20	if (itmax .gt. 0) iterma = itmax
	if (e .gt. 0.d0) epsx = e
	if (numsm1 .gt. 0) numsm = numsm1
	if (epsph1 .gt. 0.d0) epsphi = epsph1
	end

	subroutine errset(i)
	integer i
	include 'O8COMM.INC'

	if (i .gt. 0) then
		cfuerr(i) = .true.
	else
		ffuerr = .true.
		end if
	end
	subroutine user_eval(x,mode)
	integer mode
	double precision x(*)
	end
	integer function nmax()
	include 'O8PARA.INC'
	nmax = nx
	end
	integer function mmax()
	include 'O8PARA.INC'
	mmax = nresm
	end
