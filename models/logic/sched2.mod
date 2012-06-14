
# -------------------------------------------------------------
# SCHEDULING
# sched2.mod: using n integer variables 
# with "countof" operator
# -------------------------------------------------------------

param n integer > 0;

set JOBS := 1..n;
set MACHINES := 1..n;

param cap {MACHINES} integer >= 0;
param cost {JOBS,MACHINES} > 0;

var MachineForJob {JOBS} integer >= 1, <= n;

minimize TotalCost:
   sum {j in JOBS, k in MACHINES} 
      if MachineForJob[j] = k then cost[j,k];

subj to CapacityOfMachine {k in MACHINES}:
   numberof k in ({j in JOBS} MachineForJob[j]) <= cap[k];
