# -------------------------------------------------------------
# cmp_tol_big_bnd_02.mod
# As cmp_tol_big_bnd_01: x in [0, 1e12], b <==> x >= 0.001 is not
# decided true at x = 0.
# Optimum 0.001 - 1 = -0.999 at x = 0.001, b = 1.
# -------------------------------------------------------------
var x >= 0 <= 1e12;
var b binary;

minimize o: x - b;

s.t. def: b <==> x >= 1e-3;
