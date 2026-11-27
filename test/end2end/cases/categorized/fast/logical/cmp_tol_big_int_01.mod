# -------------------------------------------------------------
# cmp_tol_big_int_01.mod
# Integer x in [0, 1e15]: b <==> x >= 2 is not decided true at
# x = 0 by a margin from the large bound.
# Optimum 2 - 10 = -8 at x = 2, b = 1.
# -------------------------------------------------------------
var x integer >= 0 <= 1e15;
var b binary;

minimize o: x - 10*b;

s.t. def: b <==> x >= 2;
