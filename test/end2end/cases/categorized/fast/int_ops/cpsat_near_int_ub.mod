# An upper row bound that is an integer up to a rounding error:
# 4.35*100 evaluates to 434.99999999999994. The row means
# x + 2*y <= 435; rounding the bound down without a tolerance
# would give 434.

var x integer >= 0, <= 1000;
var y integer >= 0, <= 1000;

maximize Total: x + 2*y;

subject to AtMost: x + 2*y <= 4.35*100;
