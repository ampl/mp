
# -------------------------------------------------------------
# ASSIGNMENT
# assign3.mod: using n object-valued variables
# with "alldiff" operator and variable as index in objective
# -------------------------------------------------------------

set JOBS;
set MACHINES;
   check card (JOBS) = card (MACHINES);

param cost {JOBS,MACHINES} > 0;

var MachineForJob {JOBS} in MACHINES;

minimize TotalCost:
   sum {j in JOBS, k in MACHINES} cost[j,MachineForJob[j]];

subj to OneJobPerMachine: 
   alldiff {j in JOBS} MachineForJob[j];
