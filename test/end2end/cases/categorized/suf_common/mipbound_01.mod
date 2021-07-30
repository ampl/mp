
# -------------------------------------------------------------
# multi_sol_01.mod
# -------------------------------------------------------------

var x binary;
var y binary;
var z binary;

minimize TotalSum:
    z + 1;

subj to C1:
       x+y >= 1;

