
# -------------------------------------------------------------
# disj2ind_01.mod: implications [bvar==0/1 ==> linear relation],
#   which AMPL sends as disjunctions [bvar!=0/1 || relation],
#   to be recognized as indicator constraints.
#   See option cvt:pre:disj2ind.
# -------------------------------------------------------------

var b {1..5} binary;
var x >= -10, <= 10;
var y >= -10, <= 10 integer;

maximize Obj:
    3*b[1] + 2*b[2] + b[3] + b[4] + b[5] - x + y;

subj to C1: b[1]==1 ==> x >= 4;                # NE(b, 1) || GE
subj to C2: b[2]==0 ==> y <= 3;                # NE(b, 0) || LE
subj to C3: b[3]==1 ==> x + 2*y == 6;          # NE(b, 1) || EQ
subj to C4: b[4]==0 or y >= 8;                 # EQ(b, 0) || GE
subj to C5: b[5]==1 ==> (x <= 1 and y >= -2);  # conjunction in the body

subj to C6:                                    # conjunction of implications
    (b[1]==1 ==> y <= 5) and (b[2]==1 ==> x >= -3);

subj to C7: b[1] + b[2] + b[3] >= 2;
