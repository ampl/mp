# -------------------------------------------------------------
# exp_02_check_2con.mod.
# -------------------------------------------------------------

var x >=-100, <=100;
var y >=-100, <=100;

minimize Obj1:
       y;

subj to C1:
       y >= exp(x);

subj to C2:
       2*y <= x+2;
