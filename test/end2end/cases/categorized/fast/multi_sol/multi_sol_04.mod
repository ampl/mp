
# -------------------------------------------------------------
# multi_sol_04.mod: a feasibility problem, i.e., NO OBJECTIVE.
# Requesting multiple solutions (countsolutions, solutionstub)
# must not write objective suffix values: there is no objective
# to hold them.
# -------------------------------------------------------------

var x integer >= 0, <= 3;
var y integer >= 0, <= 3;

subj to C1:  x + y >= 3;
subj to C2:  x <> y;
