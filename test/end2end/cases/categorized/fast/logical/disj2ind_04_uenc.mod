
# -------------------------------------------------------------
# disj2ind_04_uenc.mod: implication with consequent x==k,
#   where integer x may receive the unary encoding.
#   With the encoding (default here, because 2 values
#   of x are compared in negative context), the disjunction
#   stays a row on the encoding flag;
#   without it (cvt:uenc:negctx:max=2), it is absorbed
#   into the indicator b==1 ==> x==5.
#   See option cvt:pre:disj2ind.
# -------------------------------------------------------------

var x integer >= 1, <= 9;
var y >= 1, <= 200;
var b binary;

maximize Obj: 20*b + x - 0.01*y;

subj to C1: b==1 ==> x==5;
subj to C2: (x!=2 and x!=7) or y>=4;
