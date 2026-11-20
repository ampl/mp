# -------------------------------------------------------------
# funcbnd_expa_02_neg.mod
# a^x with a negative base: not monotone and not always positive,
# so no bounds from the argument. The bounds [0, -2] were
# contradictory: a feasible model was declared infeasible.
# -------------------------------------------------------------
var n integer >= 1 <= 3;
var y;
minimize o: y;
s.t. c: y = (-2)^n;
