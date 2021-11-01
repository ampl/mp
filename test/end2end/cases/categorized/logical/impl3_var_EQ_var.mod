
# -------------------------------------------------------------
# IMPLICATION
# impl3_var_EQ_var.mod: imply x<=y
# -------------------------------------------------------------

param ub integer := 10;

var b logical;
var x integer >= -41, <= ub;
var y integer >= -41, <= ub;

maximize TotalSum:
    100*b+x -y;

subj to IMPL:
    b==1 ==> x==y;
