
# -------------------------------------------------------------
# IMPLICATION
# impl3.mod: two implications reusing a common subexpression
# -------------------------------------------------------------

param ub integer := 10;

var b logical;
var x integer >= -41, <= ub;
var y integer >= -41, <= ub;

var ce = y;

minimize TotalSum:
    b+x -y;

subj to RIMPL: 
    x<=0 ==> b==1;
    
subj to IMPL:
    b==1 ==> ce<=0;
