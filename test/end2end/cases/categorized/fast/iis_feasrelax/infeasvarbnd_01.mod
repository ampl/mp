# -------------------------------------------------------------
# infeasvarbnd_01.mod
# Contradictory bounds from singleton constraints (x >= 3, x <= 2).
# The converter declares the model infeasible (cvt:pre:infeasvarbnd),
# but passes the bounds on when alg:feasrelax or iis:find is requested.
# AMPL presolve off: MP derives the bounds.
# -------------------------------------------------------------
var x >= 0 <= 10;
var y >= 0 <= 10;
s.t. lo: x >= 3;
s.t. up: x <= 2;
s.t. c: x + y >= 4;
minimize o: x + y;
