# hs62.mod	OLR2-MY-3-7
# Original AMPL coding by Elena Bobrovnikova (summer 1996 at Bell Labs).

# Ref.: W. Hock and K. Schittkowski, Test Examples for Nonlinear Programming
# Codes.  Lecture Notes in Economics and Mathematical Systems, v. 187,
# Springer-Verlag, New York, 1981, p. 64.

# Number of variables:  3
# Number of constraints:  7
# Objective nonseparable
# Objective nonconvex
# Linear constraints

var x{1..3} >= 0, <= 1;

minimize Obj:
         -32.174 * (255 * log((x[1] + x[2] + x[3] + .03) /
         (.09 * x[1] + x[2] + x[3] + .03))
         + 280 * log((x[2] + x[3] + .03) / (.07 * x[2] + x[3] + .03))
         + 290 * log((x[3] + .03) / (.13 * x[3] + .03)));

s.t. simplex:
     x[1] + x[2] + x[3] == 1;

data;
var x :=
    1   .7
    2   .2
    3   .1;

