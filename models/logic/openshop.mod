
# Adapted from:
# www.g12.csse.unimelb.edu.au/minizinc/downloads/examples-latest/oss.mzn

param endTime integer > 0;
param nMach integer > 0;
param nJobs integer > 0;

param duration {1..nMach, 1..nJobs};

var Start {1..nMach, 1..nJobs} integer >= 0, <= endTime;
var Makespan integer >= 0, <= endTime;

minimize Objective: Makespan;

subject to NoJobConflicts 
      {m in 1..nMach, j1 in 1..nJobs, j2 in j1+1..nJobs}:
   Start[m,j1] + duration[m,j1] <= Start[m,j2]  or
   Start[m,j2] + duration[m,j2] <= Start[m,j1];

subject to NoMachineConflicts 
      {m1 in 1..nMach, m2 in m1+1..nMach, j in 1..nJobs}:
   Start[m1,j] + duration[m1,j] <= Start[m2,j]  or
   Start[m2,j] + duration[m2,j] <= Start[m1,j];

subj to MakespanDefn {m in 1..nMach, j in 1..nJobs}:
   Start[m,j] + duration[m,j] <= Makespan;
