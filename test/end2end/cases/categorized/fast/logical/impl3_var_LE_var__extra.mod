
# -------------------------------------------------------------
# IMPLICATION
# impl3_var_LE_var__extra.mod: imply x<=y
# -------------------------------------------------------------

param ub integer := 10;

var b logical;
var x integer >= -41, <= ub;
var y integer >= -41, <= ub;

maximize TotalSum:
    b+x -y;

subj to IMPL:
    b==1 ==> x<=y;

subj to RIMPL:
    x>=y ==> b==0;

subj to IMPLEQ:
    b==1 ==> x==y;

subj to RIMPLNEQ:
    x!=y ==> b==0;
