# -------------------------------------------------------------
# funcpieces_02.mod.
# Tests that the Gurobi interface correctly filters
# functional nonlinear constraint arguments:
# the FuncPieces attribute only applies to such general constraints
# and not to all general constraints.
# -------------------------------------------------------------

var x >=-100, <=100;
var y >=-100, <=100;
var z >=-100, <=100;

minimize Obj1:
    y - 2*x - 3*z;

subj to C1:
       x+y >= 1;

subj to C2:
       y + max(3, z + abs(x+3)) <= 1.83;

subj to C3:
       z + log(y + 3.8*exp(x+3)) >= -14.265;


suffix funcpieces IN;

let C1.funcpieces := 1;
let C2.funcpieces := 60;
let C3.funcpieces := 50;
