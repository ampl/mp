# A linear MIP solved with mip:basis=1, for the AMPLS repeated-solve
# test: CPLEX fixes it in the native problem, which must be the user's
# MIP again after each report.
#
# Maximum: b = 1, x = 2, y = 1, objective 6. With the sense reversed
# (minimize) the optimum is 0, with b = 0.

var b binary;
var x >= 0, <= 2;
var y >= 0, <= 10;

maximize obj: b + 2*x + y;

subject to Cap: x <= 1 + b;

# Binding linear row: dual 1 in the fixed LP
subject to Lin: x + y <= 3;
