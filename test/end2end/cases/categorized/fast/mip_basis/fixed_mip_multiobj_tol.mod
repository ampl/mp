# mip:basis with two objectives (obj:multi=1) where the first may lose
# up to 5 (objabstol) for the second.  The reported point must keep the
# values of both objectives: a fixed LP of one objective could move x.

suffix objpriority IN;
suffix objabstol IN;

var b binary;
var x {1..2} >= 0, <= 10;

maximize o1: x[1] + b;
maximize o2: x[2];

subject to C: x[1] + x[2] <= 12;
subject to D: x[2] <= 3;

let o1.objpriority := 2;
let o2.objpriority := 1;
let o1.objabstol := 5;
