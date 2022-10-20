
# -------------------------------------------------------------
# Exponent with given base
# expA_1.mod
# -------------------------------------------------------------

param ubx integer := 10;
param uby integer := 20;

var x  >= 0, <= ubx;
var y  >= -41, <= uby;

maximize TotalSum:
    x-y;

subj to ExpA: 
    y >= 1.28 ^ x;
    
