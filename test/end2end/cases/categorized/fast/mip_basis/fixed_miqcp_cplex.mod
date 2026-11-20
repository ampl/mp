# A mixed-integer problem with a quadratic constraint, solved with
# mip:basis=1: the driver fixes b at its optimal value and solves the
# continuous problem for the duals.
#
# Optimum: b = 1, x = sqrt(2), y = 3 - sqrt(2), objective 4 + sqrt(2).
#
# Dual of Lin, by hand. With b = 1 the fixed problem is
#     maximize 2*x + y  subject to  x*x <= 2 (mu),  x + y <= 3 (lambda),
# and y is strictly inside its bounds, so stationarity in y gives
#     1 - lambda = 0,  lambda = 1,
# and in x:  2 - 2*mu*x - lambda = 0,  mu = 1/(2*sqrt(2)) > 0.
# Both rows are active. Lin.dual = 1: one more unit of the right-hand
# side is used by y and adds 1 to the objective.

var b binary;
var x >= 0, <= 2;
var y >= 0, <= 10;

maximize obj: b + 2*x + y;

subject to qc: x*x <= 1 + b;

# Binding linear row: its dual comes from the fixed continuous problem
subject to Lin: x + y <= 3;
