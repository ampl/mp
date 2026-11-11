# -------------------------------------------------------------
# defvar_nl_02_linear_alias.mod
# Defined variables that are aliases of variables, used in a
# nonlinear defined variable with several product terms
# (as hs54_inspired_3_qp, but convex and not quadratic).
# -------------------------------------------------------------
var xx {1..4} >= -10, <= 10;
var x = xx[1];
var y = xx[2];
var z = xx[3];
var t = xx[4];
var h = ((5*t - 2)^2) / 4
        + ((2*x - 8)^2) + (y - z)*(y - z) + (x - 3)*(x - 3)
        + (2*x - 2*y - 5*z + 12) * (2*x - 2*y - 5*z + 12) / 10
        + exp(0.1*(x + t));
minimize Obj: h;
