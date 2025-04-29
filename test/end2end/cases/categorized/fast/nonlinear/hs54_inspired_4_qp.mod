# hs54_inspired_4_qp.mod
# Derived from hs54.mod to test summing up and out-multiplication
# of quadratics / polynomials
# -- Just 3 variables in the QP terms to invoke QP2Pass by default
# -- Using proper vars to do invoke QP2Pass (it avoids defvars)

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

var xx{5..6};
var x;
var y;
var z;
var t;
var h = ((2*x-8)^2) - (y-z)*z + (x-3)*(x-2*z+5);

minimize Obj:
      # -exp(-h);
      h;

s.t. C1:
     x + 40 * y - 1.76E+2 = 0;
s.t. B1:
     0 <= x <= 200;
s.t. B2:
     -10 <= y <= 10;
s.t. B3:
     -5 <= z <= 12;
s.t. B4:
     -15 <= t <= 20;
s.t. B5:
     -100 <= xx[5] <= 1;
s.t. B6:
     -5 <= xx[6] <= 22;
