# mip:basis on an infeasible MIP: no incumbent, so no fixed problem.
# The driver says that the fixed solve was skipped, not that it failed.

var b {1..2} binary;
var x >= 0, <= 1;

maximize o: b[1] + b[2] + x;

subject to odd: b[1] + b[2] = 1.5;
