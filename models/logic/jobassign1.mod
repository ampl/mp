
# -------------------------------------------------------------
# ASSIGNMENT
# jobassign1.mod: using n integer variables
# with "alldiff" operator
# -------------------------------------------------------------

param n integer > 0;

set JOBS := 1..n;
set MACHINES := 1..n;

param cost {JOBS,MACHINES} > 0;

var MachineForJob {JOBS} integer >= 1, <= n;

minimize TotalCost:
   sum {j in JOBS, k in MACHINES} 
      if MachineForJob[j] = k then cost[j,k];

subj to OneJobPerMachine: 
   alldiff {j in JOBS} MachineForJob[j];
