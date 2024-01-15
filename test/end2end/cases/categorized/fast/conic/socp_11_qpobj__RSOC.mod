# This is actually convex quadratic, not conic,
# but we'd first conify, then QC-fy.
# Test that this works for Mosek 10:
# the linear term must be restored.

var x {1..7} >= 0;
minimize qobj: sum {j in 1..7} x[j]^2;
subj to qconstr: sum {j in 1..6} x[j]^2 <= 6*x[7] - 1000;
