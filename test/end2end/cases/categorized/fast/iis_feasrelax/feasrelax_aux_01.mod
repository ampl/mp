# -------------------------------------------------------------
# feasrelax_aux_01.mod
# abs() is converted with auxiliary variables and constraints,
# which must not be relaxed: before, relaxing them at no cost gave
# x = 0, violating c1 by 4. Expected: c2 relaxed by 2 (cost 2),
# x = 4, c1 satisfied.
# -------------------------------------------------------------
suffix rhspen IN;
var x >= 0 <= 10;
minimize o: x;
s.t. c1: abs(x-5) <= 1;
s.t. c2: x <= 2;
let c1.rhspen := 100;
let c2.rhspen := 1;
