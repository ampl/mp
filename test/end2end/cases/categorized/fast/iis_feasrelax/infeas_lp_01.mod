# -------------------------------------------------------------
# IIS LP 01
# infeas_lp_01.mod
# The IIS consists of the upper bound of x and rows C1, C2
# -------------------------------------------------------------

var x <= 0;
var y;

minimize Total:
    x + y;

subj to C1:
     x + y >= 3;

subj to C2:
     y <= 2;

subj to C3:
     x - y <= 10;
