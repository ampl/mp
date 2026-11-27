# -------------------------------------------------------------
# numberof_ref_sort_02.mod
# As numberof_ref_sort_01: with abs() in the objective, AMPL orders
# x2, x3 before x1 in the NL file, so the reference x1 gets the
# largest index. Optimum -1 at x1 = 2, x2 = 1, x3 = 2.
# -------------------------------------------------------------
var x1 integer >= -1 <= 3;
var x2 integer >= 1 <= 2;
var x3 integer >= -3 <= 2;

minimize Obj: -x1 + 2*x2 - x3 + abs(x3 - x2);

s.t. C1: numberof x1 in (x2, x3, 1) >= 1;
