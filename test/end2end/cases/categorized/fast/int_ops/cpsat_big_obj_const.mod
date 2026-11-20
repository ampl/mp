# An objective constant above the default tech:varbound (1e12) is data,
# not a bound to cap: the model is feasible with optimum 2e12 + 10.

var x integer >= 0, <= 100;
var y integer >= 0, <= 100;

minimize Total: x + y + 2e12;

subject to AtLeast: x + y >= 10;
