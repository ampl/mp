# hs111.mod	OOR2-MN-10-23
# Original AMPL coding by Elena Bobrovnikova (summer 1996 at Bell Labs).

# Ref.: W. Hock and K. Schittkowski, Test Examples for Nonlinear Programming
# Codes.  Lecture Notes in Economics and Mathematical Systems, v. 187,
# Springer-Verlag, New York, 1981, p. 120.

# Number of variables:  10
# Number of constraints:  23
# Objective nonseparable
# Objective nonconvex
# Nonlinear constraints

set I := 1 .. 10;

param c{I};
var x{I} <= 100, >= -100, := -2.3;

minimize f:
         sum {i in I} exp(x[i]) * (c[i] + x[i] - log(sum {k in I} exp(x[k])));

s.t. C1:
    exp(x[1]) + 2 * exp(x[2]) + 2 * exp(x[3]) + exp(x[6]) + exp(x[10]) == 2;
s.t. C2:
    exp(x[4]) + 2 * exp(x[5]) + exp(x[6]) + exp(x[7]) == 1;
s.t. C3:
    exp(x[3]) + exp(x[7]) + exp(x[8]) + 2 * exp(x[9]) + exp(x[10]) == 1;

data;

param c :=
      1   -6.089      6    -14.986
      2   -17.164     7    -24.1
      3   -34.054     8    -10.708
      4   -5.914      9    -26.662
      5   -24.721    10    -22.179 ;
