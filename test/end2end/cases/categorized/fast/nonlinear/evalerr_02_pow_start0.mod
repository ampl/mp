# -------------------------------------------------------------
# evalerr_02_pow_start0.mod
# Derivative of x^0.6 undefined at the starting point (0,0)
# (coin-or/Ipopt#205): the solver must handle the evaluation
# error. Concave terms: global optimum 1.1 at (1,0).
# -------------------------------------------------------------
var x >= 0, := 0;
var y >= 0, := 0;
minimize Cost: x^0.6 + y^0.6 + 0.1*(x - 2)^2;
s.t. Demand: x + y >= 1;
