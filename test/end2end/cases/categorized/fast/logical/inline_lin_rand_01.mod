# -------------------------------------------------------------
# inline_lin_rand_01.mod
# Generated: equivalences of linear comparisons with NOT; their
# linearizations nest linear expressions (cvt:pre:unnest bit 32).
# Optimum -11 by brute force (quadratic objective term of the
# generated model dropped: MIQP is not for all solvers).
# -------------------------------------------------------------
var x{1..4} integer >= -2, <= 2;
minimize o: (2*x[1] + 1*x[4] + -2*x[3] + -1*x[2]);
s.t. c0: ((3*x[4] + 1*x[1] + -1) >= 1) <==> not ((-2*x[3] + 2*x[2] + 2) >= -1);
s.t. c1: ((-2*x[2] + -1*x[1] + 2) >= 0) <==> not ((2*x[4] + 1*x[2] + 2) <= 0);
