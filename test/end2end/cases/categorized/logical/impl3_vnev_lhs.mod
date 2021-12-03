
# -------------------------------------------------------------
# IMPLICATION
# impl3_vnev_lhs.mod: x!=y ==> ....
# -------------------------------------------------------------

param ub integer := 10;

var b logical;
var x integer >= -41, <= ub-1;
var y integer >= -41, <= ub;

maximize TotalSum:
    100*b+x + y;

subj to IMPL:
    x!=y ==> b<=0;
