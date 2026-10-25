
# -------------------------------------------------------------
# JOB-SHOP SCHEDULING
# jobshop.mod: each job is a sequence of operations, each on a
# given machine; a machine processes one operation at a time.
# Minimize the makespan.
# Disjunctive formulation with "or" (no overlap of two operations
# on the same machine).
# -------------------------------------------------------------

set JOBS ordered;
param nops integer > 0;           # operations per job
set OPS := 1..nops;

param mach {JOBS, OPS} integer > 0;       # machine of operation
param dur  {JOBS, OPS} integer > 0;       # duration of operation

param horizon := sum {j in JOBS, k in OPS} dur[j,k];

var Start {j in JOBS, k in OPS} integer >= 0, <= horizon;
var Makespan integer >= 0, <= horizon;

minimize TotalTime: Makespan;

subj to JobOrder {j in JOBS, k in 1..nops-1}:
   Start[j,k+1] >= Start[j,k] + dur[j,k];

subj to EndOfJob {j in JOBS}:
   Makespan >= Start[j,nops] + dur[j,nops];

# Two operations on the same machine do not overlap
subj to NoOverlap {j1 in JOBS, k1 in OPS, j2 in JOBS, k2 in OPS:
      ord(j1) < ord(j2) and mach[j1,k1] = mach[j2,k2]}:
   Start[j1,k1] + dur[j1,k1] <= Start[j2,k2]
   or
   Start[j2,k2] + dur[j2,k2] <= Start[j1,k1];
