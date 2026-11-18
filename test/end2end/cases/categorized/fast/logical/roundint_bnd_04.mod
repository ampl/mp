# -------------------------------------------------------------
# roundint_bnd_04.mod
# 1400000000/0.7 = 2000000000.0000002 in doubles: the derived bound
# must round to 2000000000, otherwise x == 2000000000 is decided false
# (cvt:pre:roundintbnd). Optimum 0 at x = 2000000000, y = 0.
# AMPL presolve off: MP derives the bound.
# -------------------------------------------------------------
var x integer >= 1999999950 <= 2000000050;
var y binary;
s.t. dl: 0.7*x >= 1400000000;
s.t. c: x == 2000000000 or y >= 1;
minimize o: (x - 2000000000) + 10*y;
