# Redundant rows with finite bounds beyond the 64-bit integer range
# ("big-M" right-hand sides). Such bounds cannot bind: the rows are
# accepted and the optimum is that of Cap alone, x = 0, y = 6.

param bigM := 1e19;

var x integer >= 0, <= 10;
var y integer >= 0, <= 10;

maximize Total: x + y;

subject to Cap: 3*x + 2*y <= 12;

subject to BigUpper: x - y <= bigM;
subject to BigLower: x + y >= -bigM;
