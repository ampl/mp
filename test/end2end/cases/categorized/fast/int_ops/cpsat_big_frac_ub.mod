# A fractional upper row bound of large magnitude is a real fraction,
# not rounding noise: 3e14 + 0.5 is rounded inward to 3e14, so that
# x = 1, y = 1 (row value 3e14 + 1) stays infeasible.

param big := 3e14;

var x integer >= 0, <= 10;
var y binary;

maximize Total: x + big*y;

subject to AtMost: x + big*y <= big + 0.5;
