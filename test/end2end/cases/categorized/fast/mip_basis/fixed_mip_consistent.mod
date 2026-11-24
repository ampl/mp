# mip:basis on a MIP with several optimal points: the u split the demand
# 25 among the open facilities in many ways.  The reported primal values
# must be those of the fixed problem whose basis and duals are reported:
# a nonbasic variable lies at its bound, a nonzero dual sits on an active
# row.  (Reporting the MIP incumbent with the fixed LP's basis gave, e.g.,
# u = 5 with status upp and upper bound 10.)
#
# Minimum: 3 open facilities, objective 3.

param n := 6;
set I := 1..n;

var b {I} binary;
var u {I} >= 0, <= 10;

minimize o: sum {i in I} b[i];

subject to link {i in I}: u[i] <= 10*b[i];
subject to dem: sum {i in I} u[i] >= 25;

# Neighbors share a capacity
subject to cap {i in I}: u[i] + u[if i < n then i+1 else 1] <= 14;
