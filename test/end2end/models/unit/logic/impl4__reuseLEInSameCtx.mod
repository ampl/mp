
# -------------------------------------------------------------
# IMPLICATION
# impl4.mod: two implications reusing a common subexpression
# -------------------------------------------------------------

param ub integer := 10;

var a logical;
var b logical;
var x integer >= -41, <= ub;

maximize TotalSum:
    a+b;

subj to IMPL1: 
    b==1 ==> x<=0;
    
subj to IMPL2:
    a==1 ==> x<=0;
