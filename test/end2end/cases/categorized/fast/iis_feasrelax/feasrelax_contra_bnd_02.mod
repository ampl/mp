# -------------------------------------------------------------
# feasrelax_contra_bnd_02.mod
# The auxiliary variable t = abs(x) has the deduced lower bound 0
# and gets the upper bound -1 from AbsNeg. Auxiliary variables are
# never relaxed, so the feasibility relaxation cannot repair the
# model. Expected: infeasible, not an error.
# -------------------------------------------------------------

var x >= -4, <= 4;
var y >= 0, <= 10;

minimize Total: x + y;

subject to AbsNeg: abs(x) <= -1;

subject to Cover: x + y >= 1;
