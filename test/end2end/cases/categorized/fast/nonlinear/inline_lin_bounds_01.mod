# -------------------------------------------------------------
# inline_lin_bounds_01.mod
# x - 2*y is a linear subexpression used twice (exp argument, abs argument).
# abs(x - 2*y) <= 1 gives its result variable bounds [-1, 1], stronger than
# the expression's [-6, 6]: it must stay a variable (cvt:pre:unnest bit 32).
# Objective -v + 0.1*exp(v) decreases on v <= ln 10: optimum v = 1, -1 + 0.1*e.
# -------------------------------------------------------------
var x >= -2, <= 2;
var y >= -2, <= 2;
minimize o: -(x - 2*y) + 0.1*exp(x - 2*y);
s.t. c1: abs(x - 2*y) <= 1;
