
# -------------------------------------------------------------
# IMPLICATION
# impl1.mod: binary==1 to integer variable <= 0
# -------------------------------------------------------------

param ub integer := 10;

var b logical;
var x integer >= -41, <= ub;

maximize TotalSum:
    b+x;

subj to IMPL: 
    b==1 ==> x<=0;
