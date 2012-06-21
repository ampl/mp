
# -------------------------------------------------------------
# ASSIGNMENT
# jobassign2.mod: using n integer variables
# with "alldiff" operator and variable as index in objective
# -------------------------------------------------------------

param n integer > 0;

set JOBS := 1..n;
set MACHINES := 1..n;

param cost {JOBS,MACHINES} > 0;

var MachineForJob {JOBS} integer >= 1, <= n;

minimize TotalCost:
   sum {j in JOBS, k in MACHINES} cost[j,MachineForJob[j]];

subj to OneJobPerMachine: 
   alldiff {j in JOBS} MachineForJob[j];
