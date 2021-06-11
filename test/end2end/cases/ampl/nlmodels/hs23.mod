# hs23.mod	QQR2-MN-2-9
# Original AMPL coding by Elena Bobrovnikova (summer 1996 at Bell Labs).

# Ref.: W. Hock and K. Schittkowski, Test Examples for Nonlinear Programming
# Codes.  Lecture Notes in Economics and Mathematical Systems, v. 187,
# Springer-Verlag, New York, 1981, p. 46.

# Number of variables:  2
# Number of constraints:  9
# Objective convex separable quadratic
# Nonlinear constraints


var x{1..2} <= 50, >= -50;

minimize f:
         x[1] * x[1] + x[2] * x[2];

s.t. C1:
     x[1] + x[2] - 1 >= 0;
s.t. C2:
     x[1] * x[1] + x[2] * x[2] - 1 >= 0;
s.t. C3:
     9 * x[1] * x[1] + x[2] * x[2] - 9 >= 0;
s.t. C4:
     x[1] * x[1] - x[2] >= 0;
s.t. C5:
     x[2] * x[2] - x[1] >= 0;

data;

var x :=
    1  3
    2  1;
