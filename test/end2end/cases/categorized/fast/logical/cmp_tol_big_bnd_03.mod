# -------------------------------------------------------------
# cmp_tol_big_bnd_03.mod
# x + y <= 0.001 with the body's bounds [0, 1e12 + 1]: not
# rewritten as x + y == 0.
# Optimum 1: x + y in [0.0005, 0.001] allows b = 1.
# -------------------------------------------------------------
var x >= 0 <= 1e12;
var y >= 0 <= 1;
var b binary;

maximize o: b;

s.t. def: b <==> x + y <= 1e-3;
s.t. cap: x + y >= 5e-4;
