# x has no upper bound: cvt:mip_max_bound (default 1e12) sets it, and the
# maximum lies at that bound. Both cvt:mip_max_bound warnings are expected:
# the bound was replaced, and the solution is at it.

var x integer >= 0;
var y integer >= 0, <= 10;

maximize Total: x + y;

subject to Link: y <= x;
