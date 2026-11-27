# -------------------------------------------------------------
# feasrelax_pen_03.mod
# Nothing may be relaxed (negative .rhspen and alg:lbpen/ubpen):
# the model stays infeasible. HiGHS returned a point violating
# c1 as "solved".
# -------------------------------------------------------------
suffix rhspen IN;
var x >= 0 <= 10;
var y >= 0 <= 10;
minimize o: x + y;
s.t. c1: x + y >= 2;
s.t. c2: x + y <= 1;
let c1.rhspen := -1;
let c2.rhspen := -1;
