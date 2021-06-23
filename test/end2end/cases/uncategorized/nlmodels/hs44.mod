# hs44.mod	QLR2-MN-4-10
# Original AMPL coding by Elena Bobrovnikova (summer 1996 at Bell Labs).

# Ref.: W. Hock and K. Schittkowski, Test Examples for Nonlinear Programming
# Codes.  Lecture Notes in Economics and Mathematical Systems, v. 187,
# Springer-Verlag, New York, 1981, p. 114.

# Number of variables:  4
# Number of constraints:  10
# Objective nonseparable quadratic
# Objective nonconvex
# Linear constraints


var x{1..4} >= 0;

minimize f:
    x[1] - x[2] - x[3] - x[1] * x[3] + x[1] * x[4] + x[2] * x[3] - x[2] * x[4];

s.t. C1:
     8 - x[1] - 2 * x[2] >= 0;
s.t. C2:
     12 - 4 * x[1] - x[2] >= 0;
s.t. C3:
     12 - 3 * x[1] - 4 * x[2] >= 0;
s.t. C4:
     8 - 2 * x[3] - x[4] >= 0;
s.t. C5:
     8 - x[3] - 2 * x[4] >= 0;
s.t. C6:
     5 - x[3] - x[4] >= 0;


data;

var x :=
    1   0
    2   0
    3   0
    4   0;
