# Finite variable bounds beyond the default cvt:mip_max_bound (1e12) are
# passed as given: the optimum needs x = 2e14 and y = 2e12. An infinite
# bound opposite such a bound must not make the domain empty.

var x integer >= 0, <= 2e14;
var y integer >= 2e12;           # upper bound infinite

maximize Total: x - y;

subject to Link: x - 100*y >= 0;
