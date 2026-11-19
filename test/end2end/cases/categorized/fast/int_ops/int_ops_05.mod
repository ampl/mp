# -------------------------------------------------------------
# int_ops_05.mod
# div and trunc results inside min and an implication: their
# bounds are hard bounds, so that a big-M solver behind MP2NL
# receives them ("Set bounds on variables participating in
# logical expressions" before). Optimum 1.
# -------------------------------------------------------------
var x{0..2} integer >= 0 <= 4;
maximize o: x[0] - x[1] - x[2];
s.t. c1: min(trunc(x[1]/3-1,0), x[0] + x[2] - 2) <= -1;
s.t. c2: (x[0] >= 3) ==> (trunc(x[1]/3-1,0) >= 1);
s.t. c3: min((-x[1]) div 3, x[0] + x[2] - 2) <= -1;
