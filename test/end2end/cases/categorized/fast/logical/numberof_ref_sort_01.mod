# -------------------------------------------------------------
# numberof_ref_sort_01.mod
# numberof r in (x, y): preprocessing may reorder the compared
# arguments x, y, but not the reference value r. The reference
# variable is declared last, so it has the largest index.
# numberof r in (x, y) >= 1 means r == x or r == y; r >= 4 > y.
# Optimum 4 at r = x = 4, y = 0.
# -------------------------------------------------------------
var x integer >= 0 <= 5;
var y integer >= 0 <= 3;
var r integer >= 4 <= 9;

minimize o: 2*x + y - r;

s.t. c: numberof r in (x, y) >= 1;
