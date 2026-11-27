# -------------------------------------------------------------
# round_bnd_fp_05.mod
# As round_bnd_fp_02, with round() at a half-integer: the
# argument is 3.5 at x3 = 5, x2 = -2, rounded to 4, while its
# computed upper bound is 3.4999999999999996.
# Optimum -7 at x2 = -2, x3 = 5.
# -------------------------------------------------------------
var x2 integer >= -2 <= 0;
var x3 integer >= 1 <= 5;

minimize obj: x2 - x3;

s.t. c1: round((2*x3 - x2) / 3 - 0.5) >= 4;
