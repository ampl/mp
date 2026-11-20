# A fractional coefficient of large magnitude is an error, not rounding
# noise: 1e13 + 1/64 differs from 1e13 by 8 units in the last place.

var x integer >= 0, <= 10;
var y binary;

maximize Total: x + y;

subject to Fractional: x + (1e13 + 1/64)*y <= 2e13;
