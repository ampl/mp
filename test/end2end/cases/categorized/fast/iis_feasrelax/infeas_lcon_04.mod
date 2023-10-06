# -------------------------------------------------------------
# IIS LogicalCon 04
# infeas_lcon_04.mod
# Test that AMPL can receive suffixes on logical constraints
# -------------------------------------------------------------

var x >= -30, <= 30;
var y >= -40, <= 40;
var z in interval [3, 8];

minimize TotalSum:
    x - 2*y;

subj to C1:
    x + y >= 7;

subj to C2_lcon:
    (x + z <= 6 and 5*y-z <= 5)
      or x + 2*y <= -150;

subj to C3:
    20*x +    y <= 200;

