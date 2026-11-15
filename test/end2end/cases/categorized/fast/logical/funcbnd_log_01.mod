# -------------------------------------------------------------
# funcbnd_log_01.mod: implication with log(x)
#   in a logical constraint whose argument is bounded: the result of
#   the function needs bounds derived from the argument's bounds
#   (for big-M). See option cvt:pre:funcbndtol.
# -------------------------------------------------------------

var x >= 1, <= 5;
var b binary;

minimize Obj: -x + 8*b;

subj to C1: b == 0 ==> log(x) <= 1;
