# -------------------------------------------------------------
# bounds_contradict_complement_01.mod
# abs(x-x) == 3 makes an auxiliary variable's bounds contradictory
# ([3, 0]); the if-then-else then complemented it and failed
# instead of reporting infeasibility. Infeasible model.
# -------------------------------------------------------------

var x integer >= 0 <= 3;
var a integer >= 0 <= 3;

minimize o: a;

s.t. c2: abs(x - x) == 3;
s.t. c3: (a < 0) ==> (x <= 3) else (a >= 1);
