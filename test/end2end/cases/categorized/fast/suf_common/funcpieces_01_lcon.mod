# -------------------------------------------------------------
# funcpieces_01.mod.
# Test propagation of suffixes through expression trees
# into flat constraints.
# Moreover, test IN suffixes (AMPL -> solver)
# from logical constraints.
# In this example, .funcpieces for exp(x+3) receives
# different values through C2 and C3.
# With debug=1, x-gurobi should return
# Initial.test_funcpieces_presolved = 38 (maximum)
# -------------------------------------------------------------

var x;
var y;
var z;

minimize Obj1:
       y - 2*x - 3*z;

subj to C1:
       x+y >= 1;

subj to C2:
       y + log(z + exp(x+3)) <= 1.83;

subj to C3_lcon:
       z + log(y + 3.8*exp(x+3)) >= -14.265
       or
         y - 5*z >= 12;


suffix funcpieces IN;

let C1.funcpieces := 12;
let C2.funcpieces := 23;
let C3_lcon.funcpieces := 38;
