# -------------------------------------------------------------
# knap_40_pstart.mod: knapsack with a partial MIP start (12 of 40 items)
# -------------------------------------------------------------

param n := 40;
set I := 1..n;
param w {i in I} := 10 + ((i*7919) mod 53);
param p {i in I} := 5 + ((i*104729) mod 71);

var x {I} binary;

maximize profit: sum {i in I} p[i]*x[i];

s.t. cap: sum {i in I} w[i]*x[i] <= 0.4*sum{i in I} w[i];

## Partial MIP start: items 1..12 taken, the rest unspecified
let {i in 1..12} x[i] := 1;
