# -------------------------------------------------------------
# exp_01_check_obj.mod.
# -------------------------------------------------------------

var x >=-100, <=100;
var y >=-100, <=100;

minimize Obj1:
       2*y - x;

subj to C1:
       y >= exp(x);

