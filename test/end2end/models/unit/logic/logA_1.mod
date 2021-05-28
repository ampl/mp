
# -------------------------------------------------------------
# Logarithm with given base
# logA_1.mod
# -------------------------------------------------------------

param ubx integer := 10;
param uby integer := 20;

var x  >= 0, <= ubx;
var y  >= -41, <= uby;

minimize TotalSum:
    x+y;

subj to LOGA: 
    y == -log10(x);
    
