# -------------------------------------------------------------
# cmp_tol_decide_01.mod
# x >= 4 with x <= 3.9999995 (data): within cvt:pre:feastol, the
# converter does not decide the comparison false; the solver
# accepts y = 1 with x = 3.9999995, objective 10.
# With cvt:pre:feastol=0 the comparison is false: objective 0.
# -------------------------------------------------------------
var x >= 0 <= 3.9999995;
var y binary;

maximize o: 10*y;

s.t. c: y = 1 ==> x >= 4;
