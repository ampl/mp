# hs54_inspired.mod
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

var x{1..6};
var h = ((x[1] - 1E+1)^2 / 6.4E+3 + (x[1] - 1E+1) * (x[2] - 1) / 2E+2 +
         (x[2] - 1)^2 / 3E+3) * (x[3] - 2E+2)^2 / (.96 * 4.9E+3) +
         (x[4] - 10)^2 / 2.5E+2 + (x[5] - 1E-1)^2 / 2.5E+2 +
         (x[6] - 1E+2)^2 / 2.5E+3;

minimize Obj:
       -exp(-h/2);

s.t. C1:
     x[1] + 40 * x[2] - 1.76E+2 = 0;
s.t. B1:
     0 <= x[1] <= 200;
s.t. B2:
     -10 <= x[2] <= 10;
s.t. B3:
     -5 <= x[3] <= 12;
s.t. B4:
     -15 <= x[4] <= 20;
s.t. B5:
     -100 <= x[5] <= 1;
s.t. B6:
     -5 <= x[6] <= 22;
