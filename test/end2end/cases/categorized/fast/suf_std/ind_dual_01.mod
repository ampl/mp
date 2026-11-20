# Duals of linear rows next to indicator constraints (fixed MIP).
# A driver that stores indicators as solver rows must map the duals of
# the algebraic constraints, not take the first rows. Duals of Demand:
# 1, 2, 3 (x is cheaper than y; b = 0, so the indicators are inactive).
param c {i in 1..3} := i;
param d {i in 1..3} := 2 * i;
var x {1..3} >= 0, <= 10;
var y {1..3} >= 0, <= 10;
var b {1..3} binary;
minimize Cost: sum {i in 1..3} (c[i] * x[i] + 10 * y[i] + 100 * b[i]);
subject to Ind {i in 1..3}: b[i] = 1 ==> x[i] <= 0;
subject to Demand {i in 1..3}: x[i] + y[i] >= d[i];
