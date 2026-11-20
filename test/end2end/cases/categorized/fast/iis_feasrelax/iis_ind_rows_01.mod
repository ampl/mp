# IIS with linear rows and an indicator. A driver that stores indicators
# as solver rows must map IIS rows back to the algebraic constraints.
# IIS: FixB forces b[2] = 1, Ind[2] then forces x[2] = 0, Demand cannot hold.
var x {1..3} >= 0, <= 10;
var y {1..3} >= 0, <= 10;
var z >= 0, <= 1;
var b {1..3} binary;
minimize Cost: sum {i in 1..3} (x[i] + y[i] + b[i]) + z;
subject to Ind {i in 1..3}: b[i] = 1 ==> x[i] <= 0;
subject to Pair {i in 1..3}: x[i] + y[i] >= 1;
subject to FixB: b[2] + z >= 1.5;
subject to Demand: x[2] + y[2] >= 15;
option presolve 0;
