# hs54_inspired_3_qp.mod
# Derived from hs54.mod to test summing up and out-multiplication
# of quadratics / polynomials

# hs54.mod	OLR2-MN-6-13
# Original AMPL coding by Elena Bobrovnikova (summer 1996 at Bell Labs).

# Ref.: W. Hock and K. Schittkowski, Test Examples for Nonlinear Programming
# Codes.  Lecture Notes in Economics and Mathematical Systems, v. 187,
# Springer-Verlag, New York, 1981, p. 77.

# Number of variables:  6
# Number of constraints:  13
# Objective nonseparable
# Objective nonconvex
# Linear constraints

# There is apparently a mistake in the formulation in the book which we
# corrected by substituting 6.4E+7 by 6.4E+13 in h(x). The correct result
# is different from the solution of the book.

var xx{1..6};
var x = xx[1];
var y = xx[2];
var z = xx[3];
var t = xx[4];
var h =  ((5*t-2)^2) / 4
                        + ((2*x-8)^2) - (y-z)*z + (x-3)*(x-2*z+5)
                        + (2*x -2*y - 5*z +12) * (1.2*x +3*t +4.7*z +0.3);

minimize Obj:
      # -exp(-h);
      h;

s.t. C1:
     xx[1] + 40 * xx[2] - 1.76E+2 = 0;
s.t. B1:
     0 <= xx[1] <= 200;
s.t. B2:
     -10 <= xx[2] <= 10;
s.t. B3:
     -5 <= xx[3] <= 12;
s.t. B4:
     -15 <= xx[4] <= 20;
s.t. B5:
     -100 <= xx[5] <= 1;
s.t. B6:
     -5 <= xx[6] <= 22;
