# -------------------------------------------------------------
# quadfill_nl_01.mod
# Shared quadratic subexpression (3 terms) used in 4 nonlinear
# functions: with cvt:pre:aggquadfill=0 it stays a variable defined by
# a quadratic equality; by default (no limit) it is inlined.
# -------------------------------------------------------------
var x{1..3} >= 0 <= 1;
minimize o: exp(x[1]*x[2]+x[2]*x[3]+x[1]*x[3]) + sin(x[1]*x[2]+x[2]*x[3]+x[1]*x[3])
          + cos(x[1]*x[2]+x[2]*x[3]+x[1]*x[3]) + log(1+x[1]*x[2]+x[2]*x[3]+x[1]*x[3]);
s.t. c: sum{i in 1..3} x[i] >= 1;
