# -------------------------------------------------------------
# roundint_bnd_01.mod
# Rounding of a large integer bound (cvt:pre:roundintbnd):
# 16800000.7/0.1 = 168000006.99999997 in doubles must still round
# to 168000007, the optimum. AMPL presolve off: MP derives the bound.
# -------------------------------------------------------------
var x integer >= 0 <= 1e9;
s.t. c: 0.1*x <= 16800000.7;
s.t. d: x >= 168000007;
minimize obj: x;
