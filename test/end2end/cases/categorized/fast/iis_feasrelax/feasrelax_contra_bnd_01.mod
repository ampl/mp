# -------------------------------------------------------------
# feasrelax_contra_bnd_01.mod
# The single-variable rows Lo and Up become the bounds [5, 3]
# of x (write the model with presolve 0). With lbpen = ubpen = -1
# no variable bound may be relaxed, so the feasibility relaxation
# cannot repair the model. Expected: infeasible, not an error.
# -------------------------------------------------------------

var x >= 0, <= 10;
var y >= 0, <= 10;

minimize Total: x + y;

subject to Lo: x >= 5;
subject to Up: x <= 3;

subject to Cover: x + y >= 1;
