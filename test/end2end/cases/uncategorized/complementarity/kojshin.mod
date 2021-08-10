#==> kojshin.mod

# This simple four-variable problem was given by:
# M. Kojima and S. Shindo, "Extensions of Newton and Quasi-Newton
# Method to PC^1 equations", Journal of Operations Research Society of
# Japan (29), pp. 352-374.
#
# Two solutions: x1 = (1.2247, 0, 0, 0.5), x2 = (1, 0, 3, 0).

set Rn	:=  1 .. 4 ;

var x {j in Rn};
var sx {j in Rn: j <= 2} = x[j]**2;

subject to f1:
	0 <= x[1] complements
	0 <= 3*sx[1] + 2*x[1]*x[2] + 2*sx[2] + x[3] + 3*x[4] - 6;

subject to f2:
	0 <= x[2] complements
	0 <= 2*sx[1] + x[1] + sx[2] + 10*x[3] + 2*x[4] - 2;

subject to f3:
	0 <= x[3] complements
	0 <= 3*sx[1] + x[1]*x[2] + 2*sx[2] + 2*x[3] + 9*x[4] - 9;

subject to f4:
	0 <= x[4] complements
	0 <= sx[1] + 3*sx[2] + 2*x[3] + 3*x[4] - 3;
