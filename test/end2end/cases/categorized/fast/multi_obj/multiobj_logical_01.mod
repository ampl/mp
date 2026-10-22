# -------------------------------------------------------------
# multiobj_logical_01.mod
# Lexicographic objectives on a model with logical expressions.
# The objective suffixes must stay on the objectives:
# MP2NL used to write them for the auxiliary variables and
# constraints too, and Knitro 16 (use_asl=0) misreads binary
# NL files with suffixes.
# -------------------------------------------------------------

set U := 1..4;
param cost {u in U} := 10 + 3*u;
param emis {u in U} := 9 - 2*u;
param fix {u in U} := 20 - 2*u;

var P {U} integer >= 0, <= 10;

minimize Cost:
  sum {u in U} (cost[u]*P[u] + (if P[u] >= 1 then fix[u]));
minimize Emissions:
  sum {u in U} emis[u]*P[u];

s.t. Demand: sum {u in U} P[u] >= 15;
s.t. MinRun {u in U}: P[u] >= 1 ==> P[u] >= 3;

suffix objpriority IN, integer;
suffix objreltol IN;
let Cost.objpriority := 2;
let Emissions.objpriority := 1;
let Cost.objreltol := 0.05;
