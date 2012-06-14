
# -------------------------------------------------------------
# FLOWSHOP
# flowshop1.mod: using disjunctive constraints
# -------------------------------------------------------------

set JOBS ordered;
set ALL_MACH ordered;

set MACH {JOBS} ordered within ALL_MACH;

param t_proc {i in JOBS, MACH[i]} > 0;

param t_cum {i in JOBS, j in MACH[i]} :=
   sum {jj in MACH[i]: ord(jj) <= ord(j)} t_proc[i,jj];

param t_offset {i1 in JOBS, i2 in JOBS: i1 <> i2} :=
   max {j in MACH[i1] inter MACH[i2]} 
      (t_cum[i1,j] - t_cum[i2,j] + t_proc[i2,j]);

var End >= 0;
var Start {JOBS} >= 0;

minimize Makespan: End;

subj to Makespan_Defn {i in JOBS}:
   End >= Start[i] + sum {j in MACH[i]} t_proc[i,j];

subj to No_Conflict {i1 in JOBS, i2 in JOBS: ord(i1) < ord(i2)}:
   Start[i2] >= Start[i1] + t_offset[i1,i2]  or
   Start[i1] >= Start[i2] + t_offset[i2,i1];
