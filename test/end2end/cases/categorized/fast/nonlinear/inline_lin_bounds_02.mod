# -------------------------------------------------------------
# inline_lin_bounds_02.mod
# Same without the bounding constraint: x - 2*y has no stronger bounds and is
# inlined. Optimum v = ln 10 (attainable, |v| <= 6): -ln 10 + 1.
# -------------------------------------------------------------
var x >= -2, <= 2;
var y >= -2, <= 2;
minimize o: -(x - 2*y) + 0.1*exp(x - 2*y);
s.t. c1: x + y <= 3;
