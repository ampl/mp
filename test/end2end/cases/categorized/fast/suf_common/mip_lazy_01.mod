
# -------------------------------------------------------------
# mip_lazy_01.mod
# -------------------------------------------------------------

var x binary;
var y binary;
var z binary;

minimize TotalSum:
    z + 1;

subj to C1:
       x+y >= 1;

subj to C2:
       x^2+y^2+(z-0.7)^2 <= 1.83;

subj to C3:
       z==1 ==> x-y <= 2;


suffix lazy IN;

let C1.lazy := 1;
let C2.lazy := 2;
let C3.lazy := 3;
