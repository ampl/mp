# -------------------------------------------------------------
# funcpieces_01_01_obj.mod.
# Tests propagation of suffixes through expression trees
# into flat constraints, both from constraints and objective.
# In this example, .funcpieces for exp(x+3) receives
# different values through Obj1, C2 and C3.
# With debug=1, x-gurobi should return
# Initial.test_funcpieces_presolved = 58 (maximum)
# -------------------------------------------------------------

var x;
var y;
var z;

minimize Obj1:
    y - 2*x - 3*exp(x+3);

subj to C1:
       x+y >= 1;

subj to C2:
       y + log(z + exp(x+3)) <= 1.83;

subj to C3:
       z + log(y + 3.8*exp(x+3)) >= -14.265;


suffix funcpieces IN;

let Obj1.funcpieces := 58;

let C1.funcpieces := 12;
let C2.funcpieces := 23;
let C3.funcpieces := 38;
