
# -------------------------------------------------------------
# disj2ind_03.mod: implication with a quadratic implied part.
#   Converted into a quadratic indicator constraint only with
#   bit 2 of option cvt:pre:disj2ind.
# -------------------------------------------------------------

var b {1..2} binary;
var x >= -10, <= 10;
var y >= -10, <= 10;

maximize Obj:
    10*b[1] + b[2] + x + y;

subj to C1: b[1]==1 ==> x*y <= 6;
subj to C2: b[2]==0 ==> x + y <= 3;
subj to C3: b[1] + b[2] >= 1;
subj to C4: x <= y;
