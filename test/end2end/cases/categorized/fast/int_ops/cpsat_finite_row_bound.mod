var x binary;
var y binary;

maximize obj: x + y;

subject to finite_limit: 6e17*x + 6e17*y <= 1e18;
