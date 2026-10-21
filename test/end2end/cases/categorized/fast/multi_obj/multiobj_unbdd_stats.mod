# The 1st pass is unbounded, the 2nd pass is not run.
# Test .objpass* suffixes for the skipped objective.
var x >= 0;
var y >= 0, <= 5;

maximize o1: x + y;
minimize o2: y;

suffix objpriority IN;

let o1.objpriority := 2;
let o2.objpriority := 1;
