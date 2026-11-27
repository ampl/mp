# -------------------------------------------------------------
# round_bnd_fp_02.mod
# Result bounds of floor() from floating-point argument bounds:
# the upper bound of (2*x3 - x2)/3 is computed as 3.9999999999999996,
# but the argument is exactly 4 at x3 = 5, x2 = -2 (feasible).
# Optimum -7 at x2 = -2, x3 = 5.
# -------------------------------------------------------------
var x2 integer >= -2 <= 0;
var x3 integer >= 1 <= 5;

minimize obj: x2 - x3;

s.t. c1: floor((2*x3 - x2) / 3) >= 4;
