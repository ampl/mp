# -------------------------------------------------------------
# round_bnd_fp_04.mod
# As round_bnd_fp_02, with div (truncated division): the
# quotient bound 3.9999999999999996 must not truncate to 3.
# ((2*x3 - x2)/3) div 1 = 4 at x3 = 5, x2 = -2.
# Optimum -7 at x2 = -2, x3 = 5.
# -------------------------------------------------------------
var x2 integer >= -2 <= 0;
var x3 integer >= 1 <= 5;

minimize obj: x2 - x3;

s.t. c1: ((2*x3 - x2) / 3) div 1 >= 4;
