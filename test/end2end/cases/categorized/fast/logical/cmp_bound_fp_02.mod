# -------------------------------------------------------------
# cmp_bound_fp_02.mod
# As cmp_bound_fp_01, with an equality: (x + 1)/3 == 2 at x = 5.
# Optimum -8 at x = 4, y = 0.
# -------------------------------------------------------------
var x integer >= 0 <= 5;
var y integer >= 0 <= 5;

minimize o: y - 2*x;

s.t. c: (x + 1) / 3 == 2 ==> y >= 3;
