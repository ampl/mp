
# -------------------------------------------------------------
# disj2ind_02.mod: disjunctions which are *not* convertible
#   into indicator constraints, mixed with those which are.
#   See option cvt:pre:disj2ind.
# -------------------------------------------------------------

var b {1..4} binary;
var x >= -10, <= 10 integer;
var y >= -10, <= 10;

maximize Obj:
    b[1] + b[2] + b[3] + b[4] + x + y;

subj to C1: b[1]==1 ==> x <= 2;          # indicator
subj to C2: b[2]==1 ==> x != 3;          # body not an indicator's
subj to C3: x >= 5 ==> y <= 1.5;         # condition not a binary variable
subj to C4: b[3]==1 ==> y < 4.5;         # strict inequality in the body
subj to C5:                              # only the 2nd conjunct is an indicator
    (b[4]==1 ==> x != -3) and (b[4]==1 ==> y >= -1);

subj to C6: b[1] + b[2] >= 1;
subj to C7: b[3] + b[4] >= 1;
