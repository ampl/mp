# -------------------------------------------------------------
# evalerr_03_objsqrt_lincons.mod
# sqrt'(0) in the objective at the starting point, while all
# constraints (one nonlinear, several linear) are differentiable
# there: the objective-gradient error must not be reported for
# the constraint Jacobian, and the solver must move away from
# the start. Unique optimum in the interior (constraints inactive).
# -------------------------------------------------------------
var x >= 0, <= 10, := 0;
var y >= 0, <= 10, := 1;
minimize Obj: sqrt(x) + (x - 4)^2 + (y - 2)^2;
s.t. NL: x*y >= 1;
s.t. Lin {i in 1..4}: x + i*y >= 1;
