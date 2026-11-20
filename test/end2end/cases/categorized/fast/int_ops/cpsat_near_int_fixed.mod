# A continuous variable fixed at an integer up to a rounding error:
# 2.3 - 0.3 evaluates to 1.9999999999999998 and is taken as 2,
# so that y can be used in a constraint. Then x <= 3.

var x integer >= 0, <= 5;
var y >= 2.3 - 0.3, <= 2.3 - 0.3;

maximize Total: x;

subject to Cap: x + y <= 5;
