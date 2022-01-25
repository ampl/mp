
# -------------------------------------------------------------
# obj_none_01.mod: no objective specified
# -------------------------------------------------------------

var x >= 0.0, <= 1.0;
var y >= 0.0, <= 1.0;
var t >= 0.0, <= 1.0;  ## Need t, otherwise solved by presolve
var z >= 1.0, <= 1.0;

subj to C1:
       x + y <= 1;

subj to C2:            ## TODO without t and C2, z==0 after 1st solve
       t + z <= 1;
