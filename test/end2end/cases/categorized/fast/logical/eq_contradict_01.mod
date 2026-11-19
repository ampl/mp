# -------------------------------------------------------------
# eq_contradict_01.mod
# x == 1 and x == -1: the unary encoding of x got an empty domain
# and stopped with the exception "vector" (result 500).
# -------------------------------------------------------------
var x integer >= -2 <= 3;
var y integer >= 0 <= 4;
minimize o: x + y;
s.t. c0: x == 1 and x == -1;
