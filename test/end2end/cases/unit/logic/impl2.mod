
# -------------------------------------------------------------
# IMPLICATION
# impl2.mod: integer variable <= 0 to binary==1
# -------------------------------------------------------------

param ub integer := 10;

var b logical;
var x integer >= -41, <= ub;

minimize TotalSum:
    b+x;

subj to RIMPL: 
    x<=0 ==> b==1;
