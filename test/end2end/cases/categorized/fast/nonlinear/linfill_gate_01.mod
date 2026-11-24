# -------------------------------------------------------------
# linfill_gate_01.mod
# Shared linear subexpression (6 terms) used in 4 nonlinear
# functions: its fill exceeds cvt:pre:agglinfill, so by default it stays
# a variable with a linear equality. With cvt:pre:all=0 the limit
# has no effect: the expression is written as before.
# -------------------------------------------------------------
var x{1..6} >= 0 <= 1;
minimize o: exp(sum{i in 1..6} i*x[i]) + sin(sum{i in 1..6} i*x[i])
          + cos(sum{i in 1..6} i*x[i]) + atan(sum{i in 1..6} i*x[i]);
s.t. c: sum{i in 1..6} x[i] >= 1;
