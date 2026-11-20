# -------------------------------------------------------------
# infeasvarbnd_01.mod
# Contradictory bounds from singleton constraints (x >= 3, x <= 2).
# With cvt:pre:infeasvarbnd=1 the converter declares the model
# infeasible; by default (0) it passes the bounds on to the solver.
# AMPL presolve off: MP derives the bounds.
# -------------------------------------------------------------
var x >= 0 <= 10;
var y >= 0 <= 10;
s.t. lo: x >= 3;
s.t. up: x <= 2;
s.t. c: x + y >= 4;
minimize o: x + y;
