# An equality whose right-hand side is an integer up to a rounding
# error: 0.07*100 evaluates to 7.000000000000001. Rounding both bounds
# inward without a tolerance would leave an empty range [8, 7].
# The only solution of 2*x + 3*y = 7 in the bounds is x = 2, y = 1.

var x integer >= 0, <= 10;
var y integer >= 0, <= 10;

minimize Total: x + y;

subject to Exact: 2*x + 3*y = 0.07*100;
