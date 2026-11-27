# -------------------------------------------------------------
# cmp_tol_uenc_01.mod
# Two != comparisons of an integer variable with 100001 values:
# no unary encoding (one binary per value) above
# cvt:uenc:maxdom = 10000; reified inequalities instead.
# Optimum 0 at x = 0, y = 0.
# -------------------------------------------------------------
var x integer >= 0 <= 100000;
var y binary;

minimize o: x + y;

s.t. c2: x != 4 or y >= 1;
s.t. c3: x != 1 or y >= 1;
