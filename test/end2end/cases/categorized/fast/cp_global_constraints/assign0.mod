
# -------------------------------------------------------------
# ASSIGNMENT
# assign0.mod: IP using n^2 zero-one variables
# -------------------------------------------------------------

param n integer > 0;

set JOBS := 1..n;
set MACHINES := 1..n;

param cost {JOBS,MACHINES} > 0;
var Assign {JOBS,MACHINES} binary;

minimize TotalCost:
   sum {j in JOBS, k in MACHINES} cost[j,k] * Assign[j,k];

subj to OneMachinePerJob {j in JOBS}:
   sum {k in MACHINES} Assign[j,k] = 1;

subj to OneJobPerMachine {k in MACHINES}:
   sum {j in JOBS} Assign[j,k] = 1;
