
# -------------------------------------------------------------
# SEQUENCING
# seq1a.mod: using n integer and n object-valued variables
# with variables as indices in objective & constraints
# -------------------------------------------------------------

param n integer > 0;

set JOBS;
set SLOTS := 0..n;

   check: card (JOBS) = n;
   check: "start" in JOBS;

param procTime {JOBS} > 0;
param dueTime {JOBS} > 0;
param duePenalty {JOBS} > 0;

set JOBPAIRS := {j1 in JOBS, j2 in JOBS: j1 != j2};
set PREC within JOBPAIRS;

param setupCost {JOBPAIRS} >= 0;
param setupTime {JOBPAIRS} >= 0;


var SlotForJob {JOBS} integer >= 0, <= n;
var JobForSlot {SLOTS} in JOBS;

var FinishTime {JOBS} >= 0;


minimize CostPlusPenalty:
   sum {k in SLOTS} setupCost[JobForSlot[k-1],JobForSlot[k]] +
   sum {j in JOBS} duePenalty[j] * (dueTime[j] - FinishTime[j]);


subj to Definition {j in JOBS}:
   JobForSlot[SlotForJob[j]] = j;

subj to InitialSlot: JobForSlot[0] = "start";
subj to InitialTime: FinishTime["start"] = 0;


subj to Deadline {j in JOBS}:
   FinishTime[j] <= dueTime[j];

subj to TimeNeeded {k in SLOTS}:
   FinishTime[JobForSlot[k-1]] 
      + setupTime[JobForSlot[k-1],JobForSlot[k]]
      + procTime[JobForSlot[k]] 
           <= FinishTime[JobForSlot[k]];

subj to Precedence {(j1,j2) in PREC}:
   SlotForJob[j1] < SlotForJob[j2];
