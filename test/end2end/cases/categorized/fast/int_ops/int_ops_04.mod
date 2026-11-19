# -------------------------------------------------------------
# int_ops_04.mod
# x div y under a condition, y != 0; a = |y| via two rows.
# The quotient's result bounds are hard bounds: SCIP 10.0.2
# reported 22 without them (its presolve); optimum 10 at y = -1.
# -------------------------------------------------------------
var x integer >= 0 <= 9;
var y integer >= -2 <= 2;
var q integer >= -20 <= 20;
var a integer >= 0 <= 2;
minimize o: 10*a + q + x;
s.t. c: x >= 3;
s.t. a1: a >= y;
s.t. a2: a >= -y;
s.t. g: y != 0 ==> q = x div y else q = 0;
