# -------------------------------------------------------------
# ite_nested_quad_01.mod
# Generated: nested if-then-else with a product in a branch, max and
# sin(exp(.)). Gurobi 13-14 returned a wrong "optimal" -6.29987 with
# InfUnbdInfo=1, which the driver set by default (rays); optimum -10.3
# at x = (2, 2, -2, -2), confirmed by SCIP and a fixed-point check.
# -------------------------------------------------------------
var x{1..4}  >= -2, <= 2;
minimize o: (-1*x[2] + -1*x[1] + 2*x[3] + 1*x[4]) + 0.1*(if (3*x[3] + 3*x[1] + 1*x[4] + 1) <= 1 then (if (-2*x[1] + 1*x[4] + 1) <= -1 then (x[1]*x[3] + (-1*x[3] + -1)) else (1*x[3] + -1*x[1] + 2)) else (3*x[4] + -1*x[1] + -1));
s.t. c0: (2*x[1] + 2*x[4] + 2*x[2] + 1) + max(((1*x[2] + -1*x[3] + 3*x[4]) + (3*x[2] + -1*x[4] + 1*x[3])), (-1*x[3] + -2*x[2] + 2)) >= -4;
s.t. c1: (-2*x[1] + 2*x[2] + 1) + (if (1*x[2] + -1*x[3]) <= 0 then (2*(x[1]*x[3] + (-2*x[3] + 2)) - (-2*x[2] + 2*x[1] + -2*x[4] + 2)) else (-1*x[1] + 1*x[4] + 1*x[2])) >= -6;
s.t. c2: (2*x[3] + -1*x[1] + 1*x[2]) + sin(exp((-2*x[2] + -1*x[3] + -1)/8)) >= -6;
