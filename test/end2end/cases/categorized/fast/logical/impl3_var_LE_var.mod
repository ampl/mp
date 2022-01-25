
# -------------------------------------------------------------
# IMPLICATION
# impl3_var_GE_var.mod: imply y>=x
# -------------------------------------------------------------

param ub integer := 10;

var b logical;
var x integer >= -41, <= ub;
var y integer >= -41, <= ub;

maximize TotalSum:
    b+x -y;

subj to IMPL:
    b==1 ==> y>=x;
