# -------------------------------------------------------------
# cmp_tol_big_rhs_gt_01.mod
# x > U never holds for x in [0, U]; U = 1e13 as in
# cmp_tol_big_rhs_le_01. Decided false by the converter.
# Optimum 0 at b = 0.
# -------------------------------------------------------------
param U default 1e13;

var x >= 0 <= U;
var b binary;

maximize o: b;

s.t. c: b <==> x > U;
