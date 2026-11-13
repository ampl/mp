# -------------------------------------------------------------
# IIS with a nonlinear constraint
# infeas_nl_01.mod
# The nonlinear constraint N1 is infeasible within the bounds:
# exp(x) >= 1 for x >= 0. The other constraints are irrelevant.
# Solvers keeping nonlinear constraints in a separate list
# must report N1 as an IIS member.
# -------------------------------------------------------------

var x >= 0, <= 3;
var y >= 0, <= 3;
var b binary;

minimize Obj: x + y;

s.t. A1: abs(x-1) + y >= 0.5;
s.t. N1: exp(x) + y <= 0.5;
s.t. L1: x + y >= 1;
s.t. I1: b = 1 ==> x <= 2;
s.t. N2: sin(x) * y >= -5;
