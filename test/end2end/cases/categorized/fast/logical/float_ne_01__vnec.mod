
# -------------------------------------------------------------
# x!=const for float variable
# so we use big-M or indicators
# -------------------------------------------------------------

var x >= -2e3, <= 1e3;
var y;

minimize XMinus5:
    y;

subj to NE:
    x!=(-1e2-2);

subj to YUp:
    y >= x-(-1e2-2);

subj to YDown:
    y >= (-1e2-2)-x;
