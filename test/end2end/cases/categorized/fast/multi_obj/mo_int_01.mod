# -------------------------------------------------------------
# mo_int_01.mod
# Multiple objectives on an integer-only model
# (CP solvers: CP-SAT): priorities, weights, tolerances.
# -------------------------------------------------------------

set I := 1..6;
param w {I};
param p {I};
param r {I};
param infeas default 0;        # 1: infeasible (mo_int_01_inf.dat)

var x {I} integer in 0..3;

s.t. Capacity: sum {i in I} w[i]*x[i] <= 20;
s.t. MinItems: sum {i in I} x[i] >= 4 + 100*infeas;

maximize Profit: sum {i in I} p[i]*x[i];
minimize Risk:   sum {i in I} r[i]*x[i];
maximize Spread: sum {i in I} (if x[i] >= 1 then 1);   # logical

data;
param:  w   p   r :=
   1    3   7   4
   2    4   9   6
   3    2   4   1
   4    5  11   7
   5    1   2   2
   6    6  12   3 ;
model;

suffix objpriority IN, integer;
suffix objweight IN;
suffix objabstol IN;
suffix objreltol IN;
let Profit.objpriority := 3;
let Risk.objpriority := 2;
let Spread.objpriority := 1;
