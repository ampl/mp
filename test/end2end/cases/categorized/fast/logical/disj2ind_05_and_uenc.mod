
# -------------------------------------------------------------
# disj2ind_05_and_uenc.mod: implication with a conjunction
#   whose first conjunct x==3 awaits the unary encoding
#   decision for x. The whole chain is delayed: absorbed
#   into indicators on b if x is not encoded
#   (cvt:uenc:negctx:max=2), converted as usual otherwise.
#   See option cvt:pre:disj2ind.
# -------------------------------------------------------------

var x integer >= 1, <= 9;
var y >= 0, <= 10;
var z >= 1, <= 200;
var b binary;

maximize Obj: 20*b + x + y - 0.01*z;

subj to C1: b==1 ==> (x==3 and y<=2);
subj to C2: (x!=2 and x!=7) or z>=4;
