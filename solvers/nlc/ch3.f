	block data
	common /funcom/ nvar, nobj, ncon, nzc, densej, colrow
	integer nvar, nobj, ncon, nzc, densej, colrow(14)

	common /boundc/ bounds
	double precision bounds(13)
	common /x0comn/ x0
	double precision x0(3)
	common /auxcom/ nlc
	integer nlc
	common /xkindc/ xkind
	integer xkind
	save /xkindc/
	common /pdcomn/ pd
	double precision pd(30)
	save /pdcomn/
	data nvar/3/, nobj/1/, ncon/3/, nzc/9/, densej/0/

*	*** objtype (0 = minimize, 1 = maximize) ***
	data colrow(1)/0/

*	*** colstarts ***
	data colrow(2)/1/
	data colrow(3)/4/
	data colrow(4)/7/
	data colrow(5)/10/

*	*** rownos ***
	data colrow(6)/1/
	data colrow(7)/2/
	data colrow(8)/3/
	data colrow(9)/1/
	data colrow(10)/2/
	data colrow(11)/3/
	data colrow(12)/1/
	data colrow(13)/2/
	data colrow(14)/3/
	data x0(1)/2.5d-01/
	data x0(2)/5.d-01/
	data x0(3)/7.5d-01/
	data bounds(1)/1.7e+308/
	data bounds(2)/-1.7e308/
	data bounds(3)/1.7e308/
	data bounds(4)/-1.7e308/
	data bounds(5)/1.7e308/
	data bounds(6)/-1.7e308/
	data bounds(7)/1.7e308/
	data bounds(8)/-3.333333333333333d-01/
	data bounds(9)/-3.333333333333333d-01/
	data bounds(10)/-1.d+00/
	data bounds(11)/-1.d+00/
	data bounds(12)/1.d+00/
	data bounds(13)/1.d+00/
	data nlc/2/
	data xkind/-1/
	end

	subroutine funelb(x)
	double precision x(3)
	common /xkindc/ xkind
	integer xkind
	save /xkindc/
	common /pdcomn/ pd
	double precision pd(30)
	save /pdcomn/
	double precision dv(2)

*	/*** funnel ***/

	dv(1) = pd(1)*2.d+00
	dv(1) = dv(1)*2.d+00
	dv(1) = dv(1) + pd(10)*2.d+00
	pd(16) = dv(1)

*	/*** funnel ***/

	dv(1) = pd(2)*2.d+00
	dv(1) = dv(1)*2.d+00
	dv(1) = dv(1) + pd(11)*2.d+00
	pd(17) = dv(1)

*	/*** funnel ***/

	dv(1) = pd(3)*2.d+00
	dv(1) = dv(1)*2.d+00
	dv(1) = dv(1) + pd(12)*2.d+00
	pd(18) = dv(1)

*	/*** funnel ***/

	dv(1) = pd(4)*2.d+00
	dv(1) = dv(1)*2.d+00
	dv(1) = dv(1) + pd(13)*pd(16)
	pd(19) = dv(1)

*	/*** funnel ***/

	dv(1) = pd(5)*2.d+00
	dv(1) = dv(1)*2.d+00
	dv(1) = dv(1) + pd(14)*pd(17)
	pd(20) = dv(1)

*	/*** funnel ***/

	dv(1) = pd(6)*2.d+00
	dv(1) = dv(1)*2.d+00
	dv(1) = dv(1) + pd(15)*pd(18)
	pd(21) = dv(1)
	end
	logical function xcheck(x)
	double precision x(3)
	double precision oldx(3)
	double precision v(2)
	common /xkindc/ xkind
	integer xkind
	save /xkindc/
	common /pdcomn/ pd
	double precision pd(30)
	save /pdcomn/
	if (xkind .lt. 0) then
		i = 1
		goto 20
		endif
	do 10 i = 1, 3
		if (x(i) .ne. oldx(i)) goto 20
 10		continue
	xcheck = .false.
	return
 20 	do 30 i = i, 3
 30		oldx(i) = x(i)
	xkind = 0

*	/*** defined variable 1 ***/

	pd(1) = -1.d+00 + 2.d+00*x(1)

*	/*** defined variable 2 ***/

	pd(2) = -1.d+00 + 2.d+00*x(2)

*	/*** defined variable 3 ***/

	pd(3) = -1.d+00 + 2.d+00*x(3)

*	/*** defined variable 4 ***/

	v(1) = 2.d+00 * x(1)
	v(2) = v(1) - 1.d+00
	pd(10) = 2.d+00 * v(2)
	v(2) = pd(10) * pd(1)
	v(1) = v(2) + -1.d+00
	pd(4) = v(1)

*	/*** defined variable 5 ***/

	v(1) = 2.d+00 * x(2)
	v(2) = v(1) - 1.d+00
	pd(11) = 2.d+00 * v(2)
	v(2) = pd(11) * pd(2)
	v(1) = v(2) + -1.d+00
	pd(5) = v(1)

*	/*** defined variable 6 ***/

	v(1) = 2.d+00 * x(3)
	v(2) = v(1) - 1.d+00
	pd(12) = 2.d+00 * v(2)
	v(2) = pd(12) * pd(3)
	v(1) = v(2) + -1.d+00
	pd(6) = v(1)

*	/*** defined variable 7 ***/

	v(1) = 2.d+00 * x(1)
	v(2) = v(1) - 1.d+00
	pd(13) = 2.d+00 * v(2)
	v(2) = pd(13) * pd(4)
	pd(7) = v(2)

*	/*** defined variable 8 ***/

	v(1) = 2.d+00 * x(2)
	v(2) = v(1) - 1.d+00
	pd(14) = 2.d+00 * v(2)
	v(2) = pd(14) * pd(5)
	pd(8) = v(2)

*	/*** defined variable 9 ***/

	v(1) = 2.d+00 * x(3)
	v(2) = v(1) - 1.d+00
	pd(15) = 2.d+00 * v(2)
	v(2) = pd(15) * pd(6)
	pd(9) = v(2)
	xcheck = .true.
	end

	double precision function feval(nobj, needfg, x, g)
	integer nobj, needfg
	double precision x(3), g(3)
	integer wantfg
	logical xcheck
	external xcheck

	double precision v(9), dv(9)
	common /xkindc/ xkind
	integer xkind
	save /xkindc/
	common /pdcomn/ pd
	double precision pd(30)
	save /pdcomn/
	wantfg = needfg
	if (xcheck(x) .and. wantfg .eq. 2) wantfg = 3
	if (wantfg .ge. 2) then
		call funelb(x)
		endif
	if (mod(wantfg,2) .eq. 1) then

*	/*** defined variable 10 ***/

	v(7) = pd(7) + 1.d+00
	v(7) = v(7) - 2.d+00*x(1)

*	/*** defined variable 11 ***/

	v(8) = pd(8) + 1.d+00
	v(8) = v(8) - 2.d+00*x(2)

*	/*** defined variable 12 ***/

	v(9) = pd(9) + 1.d+00
	v(9) = v(9) - 2.d+00*x(3)

*  /***  objective ***/

	pd(22) = pd(1) + pd(2)
	pd(22) = pd(22) + pd(3)
	v(4) = 3.333333333333333d-01 * pd(22)
	pd(23) = v(4) * v(4)
	pd(24) = v(4) + v(4)
	v(4) = 5.d-01 * pd(23)
	pd(25) = pd(4) + pd(5)
	pd(25) = pd(25) + pd(6)
	v(5) = 3.333333333333333d-01 * pd(25)
	v(6) = v(5) - -3.333333333333333d-01
	pd(26) = v(6) * v(6)
	pd(27) = v(6) + v(6)
	v(6) = 5.d-01 * pd(26)
	v(4) = v(4) + v(6)
	pd(28) = v(7) + v(8)
	pd(28) = pd(28) + v(9)
	v(6) = 3.333333333333333d-01 * pd(28)
	pd(29) = v(6) * v(6)
	pd(30) = v(6) + v(6)
	v(6) = 5.d-01 * pd(29)
	v(4) = v(4) + v(6)
	feval = v(4)

	endif

	if (wantfg .gt. 1) then
	dv(1) = 5.d-01*pd(30)
	dv(1) = dv(1)*3.333333333333333d-01
	dv(2) = dv(1)
	dv(3) = dv(1)
	dv(4) = 5.d-01*pd(27)
	dv(4) = dv(4)*3.333333333333333d-01
	dv(5) = dv(4)
	dv(6) = dv(4)
	dv(7) = 5.d-01*pd(24)
	dv(7) = dv(7)*3.333333333333333d-01
	dv(8) = dv(7)
	dv(9) = dv(7)
	g(3) = -dv(2)*2.d+00
	g(2) = -dv(3)*2.d+00
	g(1) = -dv(1)*2.d+00
	g(3) = g(3) + dv(2)*pd(21)
	g(2) = g(2) + dv(3)*pd(20)
	g(1) = g(1) + dv(1)*pd(19)
	g(3) = g(3) + dv(5)*pd(18)
	g(2) = g(2) + dv(6)*pd(17)
	g(1) = g(1) + dv(4)*pd(16)
	g(3) = g(3) + dv(8)*2.d+00
	g(2) = g(2) + dv(9)*2.d+00
	g(1) = g(1) + dv(7)*2.d+00

	endif
	end

	subroutine ceval(needfg, x, c, J)
	integer needfg
	double precision x(3), c(3), J(1)
	integer wantfg
	logical xcheck
	external xcheck

	double precision v(2), dv(1)
	double precision t1
	common /xkindc/ xkind
	integer xkind
	save /xkindc/
	common /pdcomn/ pd
	double precision pd(30)
	save /pdcomn/
	wantfg = needfg
	if (xcheck(x) .and. wantfg .eq. 2) wantfg = 3

	if (mod(wantfg,2) .gt. 0) then

*  /***  constraint 1  ***/

	v(1) = 3.333333333333333d-01 * pd(6)
	v(2) = 3.333333333333333d-01 * pd(5)
	v(1) = v(1) + v(2)
	v(2) = 3.333333333333333d-01 * pd(4)
	v(1) = v(1) + v(2)
	c(1) = v(1)

*  /***  constraint 2  ***/

	v(1) = 3.333333333333333d-01 * pd(9)
	v(2) = 3.333333333333333d-01 * pd(8)
	v(1) = v(1) + v(2)
	v(2) = 3.333333333333333d-01 * pd(7)
	v(1) = v(1) + v(2)
	t1 = v(1) + -6.666666666666666d-01*x(1)
	t1 = t1 + -6.666666666666666d-01*x(2)
	t1 = t1 + -6.666666666666666d-01*x(3)
	c(2) = t1

*  /***  constraint 3  ***/

	t1 = 6.666666666666666d-01*x(1)
	t1 = t1 + 6.666666666666666d-01*x(2)
	t1 = t1 + 6.666666666666666d-01*x(3)
	c(3) = t1
	endif

	if (wantfg .gt. 1) then
		call funelb(x)

*   /*** derivatives for constraint 1 ***/

	J(7) = 3.333333333333333d-01*pd(18)
	J(4) = 3.333333333333333d-01*pd(17)
	J(1) = 3.333333333333333d-01*pd(16)

*   /*** derivatives for constraint 2 ***/

	J(8) = 3.333333333333333d-01*pd(21) + -6.666666666666666d-01
	J(5) = 3.333333333333333d-01*pd(20) + -6.666666666666666d-01
	J(2) = 3.333333333333333d-01*pd(19) + -6.666666666666666d-01

*   /*** derivatives for constraint 3 ***/

	J(3) = 6.666666666666666d-01
	J(6) = 6.666666666666666d-01
	J(9) = 6.666666666666666d-01
	endif
	end
