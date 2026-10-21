
# -------------------------------------------------------------
# disj2ind_06_compl.mod: complementarity with a single
#   finite bound is converted into the disjunction
#   x<=0 || y-1<=0 of two conditional comparisons.
#   With cvt:pre:disj2ind bit 4, one of them is absorbed
#   into an indicator on the other one's result:
#   a single binary.
# -------------------------------------------------------------

var x >= 0;
var y >= -10, <= 10;

maximize Obj: x + y;

subj to Lim: x <= 10;
subj to C: x >= 0 complements y - 1 >= 0;
