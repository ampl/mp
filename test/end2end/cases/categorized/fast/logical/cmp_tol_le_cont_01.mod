# -------------------------------------------------------------
# cmp_tol_le_cont_01.mod
# A non-strict comparison near a bound of a continuous variable
# is kept as written: x <= 5e-5 with x in [0, 1e-4] is not
# x == 0 (cvt:mip:eps applies to strict comparisons only).
# Optimum 0.5 at x = 5e-5, y = 0.
# -------------------------------------------------------------
var x >= 0 <= 1e-4;
var y binary;

maximize o: 1e4*x - 100*y;

s.t. c: x <= 5e-5 or y >= 1;
