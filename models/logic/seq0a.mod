
# -------------------------------------------------------------
# SEQUENCING
# seq0a.mod: IP using n^2 zero-one variables,
# with arbitrary job set
# -------------------------------------------------------------

set JOBS;
   check: "s" not in JOBS;  # start of sequence
   check: "t" not in JOBS;  # end of sequence

param procTime {JOBS} >= 0;
param dueTime {JOBS} >= 0;
param duePenalty {JOBS} >= 0;

set JOBPAIRS 
   := {j1 in {"s"} union JOBS, 
       j2 in {"t"} union JOBS: j1 <> j2} diff {("s","t")};
set PREC within JOBPAIRS;

param setupCost {JOBPAIRS} >= 0;
param setupTime {JOBPAIRS} >= 0;

param BIG := max {j in JOBS} dueTime[j];


var FinishTime {JOBS} >= 0 integer;

var Seq {JOBPAIRS} binary;


minimize CostPlusPenalty:
   sum {(j1,j2) in JOBPAIRS} setupCost[j1,j2] * Seq[j1,j2] +
   sum {j in JOBS} duePenalty[j] * (dueTime[j] - FinishTime[j]);


subj to OneBefore {j in {"t"} union JOBS}:
   sum {(j1,j) in JOBPAIRS} Seq[j1,j] = 1;

subj to OneAfter {j in {"s"} union JOBS}:
   sum {(j,j2) in JOBPAIRS} Seq[j,j2] = 1;


subj to Deadline {j in JOBS}:
   FinishTime[j] <= dueTime[j];

subj to TimeNeeded {(j1,j2) in JOBPAIRS: j2 <> "t"}:
   (if j1 <> "s" then FinishTime[j1]) + setupTime[j1,j2] + procTime[j2] 
      <= FinishTime[j2] + BIG * (1 - Seq[j1,j2]);

subj to Precedence {(j1,j2) in PREC}: Seq[j2,j1] = 0;


data;

set JOBS :=  1 2 3 4 ;

param: procTime dueTime duePenalty :=
   1      6       21      9
   2      6       33      9
   3      3       21      2
   4      9       33      6 ;

set PREC :=  1 2  3 4 ;

param setupTime:  1   2   3   4   t :=
              s   3   3   2   2   .
              1   .   0   1   1   0
              2   0   .   1   1   0
              3   3   3   .   0   0
              4   3   3   0   .   0 ;

param setupCost:  1   2   3   4   t :=
              s  90  90  60  60   .
              1   .   0  30  30   0
              2   0   .  30  30   0
              3  90  90   .   0   0
              4  90  90   0   .   0 ;
