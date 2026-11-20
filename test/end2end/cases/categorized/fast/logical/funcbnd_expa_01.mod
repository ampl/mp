# -------------------------------------------------------------
# funcbnd_expa_01.mod: implication with a^x, base below 1
#   in a logical constraint whose argument is bounded: the result of
#   the function needs bounds derived from the argument's bounds
#   (for big-M). See option cvt:pre:funcbndtol.
# -------------------------------------------------------------

var x >= 0, <= 5;
var b binary;

minimize Obj: -x + 8*b;

subj to C1: b == 0 ==> 0.5^x >= 0.25;
