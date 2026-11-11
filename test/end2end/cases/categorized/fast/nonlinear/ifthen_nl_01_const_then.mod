# -------------------------------------------------------------
# ifthen_nl_01_const_then.mod
# if-then-else with a constant "then" branch that is taken at
# the optimum, added to other nonlinear terms: the gradient of
# those terms must not be lost. The penalty max(0, 1-y)^2 is
# continuously differentiable; the problem is convex.
# Optimum: x = 2, y = 3, objective 1.
# -------------------------------------------------------------
var x >= -5, <= 5, := 0;
var y >= -5, <= 5, := 2;
minimize Obj: (x - 2)^2 + exp(y - 3) - (y - 3)
              + (if y >= 1 then 0 else (y - 1)^2);
