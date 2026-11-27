# -------------------------------------------------------------
# inline_lin_rand_02.mod
# Generated: nested NOT in equivalences, sqrt(1 + sin(.)^2) of a
# linear expression (cvt:pre:unnest bit 32).
# Optimum -11.899009165927 by brute force.
# -------------------------------------------------------------
var x{1..4} integer >= -2, <= 2;
minimize o: (2*x[1] + 2*x[3] + 3*x[4] + -2*x[2] + 2) + 0.1*sqrt(1 + (sin((-1*x[4] + 3*x[2] + 3*x[1] + 1)))^2);
s.t. c0: (not ((-1*x[1] + 1*x[3] + -2*x[4]) <= -1)) <==> not (not ((3*x[4] + 1*x[3]) <= 1));
s.t. c1: ((-2*x[3] + 3*x[2] + 2) <= 1) <==> not ((3*x[4] + -2*x[1] + 3*x[2] + 2) <= 2);
s.t. c2: (not ((2*x[3] + -2*x[4] + 1) <= 1)) <==> not ((1*x[1] + -1*x[2]) <= 2);
