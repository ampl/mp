
# -------------------------------------------------------------
# obj_suf_01.mod
# -------------------------------------------------------------

var x >= 0.0, <= 1.0;
var y >= 0.0, <= 1.0;
var z >= 0.0, <= 1.0;

maximize XY_Z:
    x + y - z;
    
minimize Y:
    y;

subj to C1:
       x + y <= 1;

suffix objpriority IN, integer;
suffix objweight IN, >=-1e20, <=1e20;

let XY_Z.objpriority := 1;

let XY_Z.objweight := 1;
let Y.objweight := 1.5;
