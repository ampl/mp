
# -------------------------------------------------------------
# Exponent with given base
# expA_1.mod
# -------------------------------------------------------------

param ubx integer := 100;
param uby integer := 200;

var x  >= -100, <= ubx;
var y  >= -205, <= uby;

maximize TotalSum:
    -x-y;

subj to ExpA: 
    y >= 0.1 ^ x;
    
