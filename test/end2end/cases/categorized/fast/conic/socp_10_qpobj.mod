# From #229.
# This is actually convex quadratic, not conic.

var x {1..7} >= 0;
minimize qobj: sum {j in 1..7} x[j]^2;
subj to qconstr: sum {j in 1..7} x[j]^2 <= 1000;
