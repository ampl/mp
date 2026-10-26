# -------------------------------------------------------------
# ws_hint_01.mod
# MIP start on an integer model with an auxiliary variable
# (count). The start (x[1], x[2]) is feasible, value 25;
# the optimum is 62. CP-SAT cases fix the variables to the
# hints to see that they arrive.
# -------------------------------------------------------------

set I := 1..12;
param w {i in I} := 3 + (7*i) mod 11;
param v {i in I} := 5 + (5*i) mod 13;
var x {I} binary;
s.t. Few: count {i in I} (x[i] = 1) <= 4;   # auxiliary result variable
maximize Value: sum {i in I} v[i]*x[i];
s.t. Cap: sum {i in I} w[i]*x[i] <= 30;
let {i in I} x[i] := if i <= 2 then 1 else 0;   # feasible, not optimal
