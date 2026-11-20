# A lower row bound that is an integer up to a rounding error:
# 0.1*3*10 evaluates to 3.0000000000000004. The row means x + y >= 3;
# rounding the bound up without a tolerance would give x + y >= 4.

var x integer >= 0, <= 10;
var y integer >= 0, <= 10;

minimize Total: x + y;

subject to AtLeast: x + y >= 0.1*3*10;
