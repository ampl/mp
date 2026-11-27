# mip:basis with two objectives (obj:multi=1): the second picks among the
# alternative optima of the first.  The reported point must keep the
# value of the second objective (the fixed model of the first objective
# alone could move x and lose it).

var b {1..3} binary;
var x {1..3} >= 0, <= 10;

maximize o1: sum {i in 1..3} b[i];
maximize o2: x[1] - x[2];

subject to lnk {i in 1..3}: x[i] <= 10*b[i];
subject to dem: sum {i in 1..3} x[i] >= 12;
subject to cap: sum {i in 1..3} x[i] <= 15;
