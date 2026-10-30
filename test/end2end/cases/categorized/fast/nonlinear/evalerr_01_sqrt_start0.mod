# -------------------------------------------------------------
# evalerr_01_sqrt_start0.mod
# Derivative of sqrt undefined at the starting point (0,0):
# the solver must handle the evaluation error (ASL returns it
# to the caller) and move away from the start.
# Optimum: sqrt(5) - 1/4 = 1.9860679775 at r = sqrt(5) - 1/2.
# -------------------------------------------------------------
var x := 0;
var y := 0;
minimize Obj: sqrt(x^2 + y^2) + (x - 1)^2 + (y - 2)^2;
s.t. Ring: sqrt(x^2 + y^2) >= 0.5;
