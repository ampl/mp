# -------------------------------------------------------------
# cmp_tol_eq_feastol_01.mod
# x == 0.9999995 with x <= 1: the right-hand side is within
# cvt:pre:feastol of the bound, but not within rounding error,
# so it is kept (not moved to the bound 1).
# Optimum 0 at x = 0.9999995, w = 0, z = 0.
# -------------------------------------------------------------
var x >= 0 <= 1;
var w binary;
var z >= 0;

minimize o: z + 100*w;

s.t. c: x == 0.9999995 or w = 1;
s.t. d: z >= 1e7*x - 9999995;      # z = 5 at x = 1
