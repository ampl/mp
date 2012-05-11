* /****************************************************************
* Copyright (C) 1997 Lucent Technologies
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

	subroutine calcr(n, p, x, nf, r, ui, ur, uf)
	integer n, p, nf, ui(*)
*	n observations, p parameters
	double precision x(p), r(n), ur(*)
	external uf

	integer i, nerror

	nerror = 0
	call conval(n, p, x, r, nerror)
	if (nerror .ne. 0) then
		nf = 0
	else
		do 10 i = 1, n
 10			r(i) = r(i) - ur(i)
		endif
	end

	subroutine calcj(n, p, x, nf, jac, ui, ur, uf)
	integer n, p, nf, ui(*)
	double precision x(p), jac(n,p), ur(*)
	external uf

	integer nerror

	nerror = 0
	call jacval(n, p, ui(1), x, jac, nerror)
	if (nerror .ne. 0) nf = 0
	end

	program fnl21
	integer liv, lv, nmax, pmax, nzmax
	parameter(nmax = 200, pmax = 20, nzmax = nmax*pmax)
	parameter(liv = 82 + 4*pmax)
	parameter(lv = 105 + pmax*(nmax + 2*pmax + 21) + 2*nmax)
	integer i, n, no, nz, p
	character*80 stub
	double precision b(2,pmax), d(pmax), v(lv), x(pmax)
	integer iv(liv)
	external calcj, calcr
	character*40 endmsg(2)
	character*5 arg2
	logical amplfl

	integer prunit
	parameter (prunit = 21)

*** The following are only needed to satisfy calling sequences...
	integer mxrow, mxcol
	integer jp(nmax+1)
	integer*2 ji(nzmax)
	double precision l(pmax), u(pmax), lrhs(nmax), urhs(nmax)
	double precision inf

	call getarg(1,stub)
	call jacdim(stub, n, p, no, nz, mxrow, mxcol)
	if (n .gt. nmax .or. p .gt. pmax) then
		write(*,*) 'Problem too large: ',n,' equations,',
     1				p, ' variables'
		stop
		endif
	call jacinc(n, p, nz, jp, ji, x, l, u, lrhs, urhs, inf)
	call densej
	do 10 i = 1, n
		d(i) = 1
		b(1,i) = l(i)
		b(2,i) = u(i)
 10		continue
	call divset(1, iv, liv, lv, v)
	amplfl = .false.
	call getarg(2,arg2)
	if (arg2 .eq. '-AMPL') then
		amplfl = .true.
		iv(prunit) = 0
		endif
	call dn2gb(n, p, x, b, calcr, calcj, iv, liv, lv, v, nz, lrhs,
     1		   calcr)
	write(endmsg, '(''dn2gb return code ='',i3/''final f ='',
     1 1p,g18.10)') iv(1), v(10)
	if (amplfl) then
		call wrtsol(endmsg,2,x,urhs)
	else
		write(*,*) endmsg
		endif
	end
