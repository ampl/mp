
# -------------------------------------------------------------
# IMPLICATION
# impl3_vev_lhs.mod: x==y ==> ....
# -------------------------------------------------------------

param ub integer := 10;

var b logical;
var x integer >= -41, <= ub;
var y integer >= -41, <= ub;

maximize TotalSum:
    10*b+x + y;

subj to IMPL:
    x==y ==> b<=0;
