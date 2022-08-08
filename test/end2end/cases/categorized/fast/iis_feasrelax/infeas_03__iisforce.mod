
# -------------------------------------------------------------
# IIS 03 + IISForce
# infeas_03__iisforce.mod
# -------------------------------------------------------------

var x <= 1;
var y >= 1;
var z <= 1;

minimize TotalSum:
    x - 2*y;


subj to Cxy0:
       y + 0.0001*x >= 1;

subj to Cxy1:
       x + y <=  2;

subj to Cxy2:
       x + y <=  1.8;

subj to Cx_y1:
       x - y >= -0.5;

subj to Cx_y2:
       x - y >= -0.1;

subj to Cx_y3:
       1 >= x - y >=  0.5;

subj to Cx_y4:
       1.52 >= x - y >=  1.23;


subj to Cxz1:
       x + z <=  1;

subj to Cx_z1:
       x - z >=  0.5;


suffix iisubforce IN;
let x.iisubforce := 1; 

suffix iislbforce IN;
let Cx_y4.iislbforce := -1; 

suffix iisforce IN;
let Cxy0.iisforce := -1; 

