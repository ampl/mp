
# -------------------------------------------------------------
# mip_round_01.mod
# -------------------------------------------------------------

var x;
var b binary;

maximize TotalSum:
    x - 2e7*b;

subj to C1:
       x <= 1e6*b;

