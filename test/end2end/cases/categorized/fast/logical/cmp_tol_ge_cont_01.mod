# -------------------------------------------------------------
# cmp_tol_ge_cont_01.mod
# As cmp_tol_le_cont_01 at the upper bound: x >= 0.99995 with
# x in [0, 1] is not x == 1.
# Optimum 0.99995 at x = 0.99995, y = 0.
# -------------------------------------------------------------
var x >= 0 <= 1;
var y binary;

minimize o: x + 100*y;

s.t. c: x >= 0.99995 or y >= 1;
