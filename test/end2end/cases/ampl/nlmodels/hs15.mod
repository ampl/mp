# hs15.mod	OQR2-MN-2-3
# Original AMPL coding by Elena Bobrovnikova (summer 1996 at Bell Labs).

# Ref.: W. Hock and K. Schittkowski, Test Examples for Nonlinear Programming
# Codes.  Lecture Notes in Economics and Mathematical Systems, v. 187,
# Springer-Verlag, New York, 1981, p. 38.

# Number of variables:  2
# Number of constraints:  3
# Objective nonseparable
# Objective nonconvex
# Nonlinear constraints


var x{1..2};

minimize Obj:
         100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2;

s.t. Ineq1:
     x[1] * x[2] - 1 >= 0;
s.t. Ineq2:
     x[1] + x[2] * x[2] >= 0;
s.t. Bound1:
     x[1] <= .5;

data;

var x :=
    1   .4
    2    1.5;
