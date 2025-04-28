
# -------------------------------------------------------------
# Natural logarithm
# log_1_ccv.mod
# Concave usage
# -------------------------------------------------------------

param ubx integer := 10;
param uby integer := 20;

var x  >= 0, <= ubx;
var y  >= -41, <= uby;

minimize TotalSum:
    x+y;

subj to LOG: 
    y >= -log(x);
    
