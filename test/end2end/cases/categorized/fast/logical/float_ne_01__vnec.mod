
# -------------------------------------------------------------
# x!=const for float variable
# so we use big-M or indicators
# -------------------------------------------------------------

var x >= -2e9, <= -1e9;
var y;

minimize XMinus5:
    y+3;

subj to NE:
    x!=(-1e9-2);

subj to YUp:
    y >= x-(-1e9-2);

subj to YDown:
    y >= (-1e9-2)-x;
