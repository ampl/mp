
# -------------------------------------------------------------
# Reified equality comparison var==const
# eqVarConst01.mod
# -------------------------------------------------------------

var y integer >=-3, <=-1;
var x integer >= -41, <= -38;

minimize TotalSum:
    y+x;

subj to RIMPL: 
    y==-3 ==> x==-39;
