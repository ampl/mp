# A coefficient that is an integer up to a rounding error:
# 0.1*3*10 evaluates to 3.0000000000000004 and is taken as 3.
# With 3*x + 2*y <= 12 the optimum is x = 4, y = 0.

param a := 0.1*3*10;

var x integer >= 0, <= 10;
var y integer >= 0, <= 10;

maximize Total: 2*x + y;

subject to Cap: a*x + 2*y <= 12;
