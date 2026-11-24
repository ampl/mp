# mip:basis after a MIP stopped at its first incumbent (lim:sol=1): the
# fixed LP improves the continuous part.  The reported objective value
# must be that of the reported point, the fixed LP's (before: the point
# of the fixed LP with the incumbent's objective value, and the solution
# check warned).

param n := 30;

var b {1..n} binary;
var x {1..n} >= 0, <= 10;

maximize o: sum {i in 1..n} b[i] + 0.01*sum {i in 1..n} (i mod 7)*x[i];

subject to lnk {i in 1..n}: x[i] <= 10*b[i];

# Neighbors open together only with enough flow
subject to pair {i in 1..n-1}:
    b[i] + b[i+1] <= 1 + 0.5*(x[i] + x[i+1])/10;

subject to dem: sum {i in 1..n} x[i] >= 12;
subject to cap: sum {i in 1..n} x[i] <= 150;
