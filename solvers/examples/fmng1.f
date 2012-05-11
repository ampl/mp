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
	subroutine calcf(n, x, nf, f, ui, ur, uf)
	integer n, nf, ui(*)
	double precision x(n), f, ur(*)
	external uf
	integer nerror
	double precision objval

	nerror = 0
	f = objval(n, x, 0, nerror)
	if (nerror .ne. 0) nf = 0
	end

	subroutine calcg(n, x, nf, g, ui, ur, uf)
	integer n, nf, ui(*)
	double precision x(n), g(n), ur(*)
	external uf
	integer nerror

	nerror = 0
	call objgrd(n, x, 0, g, nerror)
	if (nerror .ne. 0) nf = 0
	end

	program fmng1
	integer liv, lv, mmax, nmax, nzmax
	parameter(mmax = 20, nmax = 20, nzmax = 100)
	parameter(liv = nmax + 59, lv = 71 + nmax*(nmax+21)/2)
	integer i, m, n, no, nz
	character*80 stub
	double precision b(2,nmax), d(nmax), v(lv), x(nmax)
	integer iv(liv)
	external calcf, calcg
	character*40 endmsg(2)
*** The following are only needed to satisfy calling sequences...
	integer mxrow, mxcol
	integer jp(nmax+1), ui(1)
	integer*2 ji(nzmax)
	integer ui(1)
	double precision l(nmax), u(nmax), lrhs(mmax), urhs(mmax)
	double precision inf, ur(1)

	call getarg(1,stub)
	call jacdim(stub, m, n, no, nz, mxrow, mxcol)
	if (m .gt. 0) then
		if (nz .le. nzmax .and. m .le. mmax) then
			write(*,*) 'Ignoring ',m,' constraints.'
		else
			write(*,*) 'Too many constraints to ignore'
			stop
			endif
		endif
	if (n .gt. nmax) then
		write(*,*) 'Too many variables (',n,') -- more than ',
     1			nmax
		stop
		endif

	call jacinc(m, n, nz, jp, ji, x, l, u, lrhs, urhs, inf)
	do 10 i = 1, n
		d(i) = 1
		b(1,i) = l(i)
		b(2,i) = u(i)
 10		continue
	iv(1) = 0
	call dmngb(n, d, x, b, calcf, calcg, iv, liv, lv, v,
     1		   ui, ur, calcf)
	write(endmsg,'(''dmngb return code ='',i3/''final f ='',
     1 1p,g18.10)') iv(1), v(10)
	call wrtsol(endmsg,2,x,urhs)
	end
