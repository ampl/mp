
# -------------------------------------------------------------
# x!=const for integer variable
# -------------------------------------------------------------

var x integer >= -2e2, <= 1e2;
var y integer;

minimize XMinus5:
    y+3;

subj to NE:
    x!=(-12);

subj to YUp:
    y >= x-(-12);

subj to YDown:
    y >= (-12)-x;
