# -------------------------------------------------------------
# cmp_tol_big_bnd_01.mod
# A large data bound (1e20) does not widen the converter's
# decision margin: x >= 5 is not decided true at x = 3.
# The margin comes from the compared constant 5 only.
# Optimum 3 at x = 3, y = 0.
# -------------------------------------------------------------
var x >= 0 <= 1e20;
var y binary;

minimize o: x + 100*y;

s.t. setup: x >= 5 ==> y = 1;      # fixed charge above 5 units
s.t. need: y = 0 ==> x <= 4;
s.t. dem: x >= 3;
