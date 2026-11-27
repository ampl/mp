# -------------------------------------------------------------
# cmp_tol_strict_false_01.mod
# Integer x, b <==> x < 1e-7 with b = 1: the right-hand side
# within cvt:pre:inteps of 0 is 0, and x < 0 is false at x >= 0,
# so the model is infeasible, while x = 0 satisfies x < 1e-7
# exactly. A warning (strict_cmp_tol) says so.
# With cvt:pre:inteps=-1: optimum 0.
# -------------------------------------------------------------
var x integer >= 0 <= 6;
var b binary;

minimize o: x;

s.t. d: b <==> (x < 1e-7);
s.t. f: b >= 1;
