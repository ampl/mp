var x binary;
var y binary;

maximize obj: x;

subject to fractional: 1000000000.5*x + 2*y <= 1000000000.75;
