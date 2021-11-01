
# -------------------------------------------------------------
# x!=const for integer variable
# -------------------------------------------------------------

param ub integer := 10;

var x integer >= -41, <= ub;
var y integer >= -41, <= ub;

minimize XMinus5:
    y+3;

subj to NE:
    x!=5;

subj to YUp:
    y >= x-5;

subj to YDown:
    y >= 5-x;
