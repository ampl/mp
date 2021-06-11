# hs112.mod	OLR2-MY-10-13
# Original AMPL coding by Elena Bobrovnikova (summer 1996 at Bell Labs).

# Ref.: W. Hock and K. Schittkowski, Test Examples for Nonlinear Programming
# Codes.  Lecture Notes in Economics and Mathematical Systems, v. 187,
# Springer-Verlag, New York, 1981, p. 121.

# Number of variables:  10
# Number of constraints:  13
# Objective partially separable
# Objective nonconvex
# Linear constraints

# There is a mistake in the formulation of the problem: constraint C3 can never
# hold with the given bounds on x. Therefore constraint C3 has been modified
# to make the problem feasible.

set I := 1 .. 10;

param c{I};
var x{I} >= 1e-6, := .1;

minimize Obj:
         sum {i in I} x[i] * (c[i] + log(x[i] / (sum {k in I} x[k])));

s.t. C1:
     x[1] + 2 * x[2] + 2 * x[3] + x[6] + x[10] == 2;
s.t. C2:
     x[4] + 2 * x[5] + x[6] + x[7] == 1;
s.t. C3:
     x[3] + x[7] + x[8] + 2 * x[9] + x[10] == 1.5;

data;

param c :=
      1   -6.089      6    -14.986
      2   -17.164     7    -24.1
      3   -34.054     8    -10.708
      4   -5.914      9    -26.662
      5   -24.721    10    -22.179 ;



