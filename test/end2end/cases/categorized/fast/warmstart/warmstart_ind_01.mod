# Primal and dual initial values on a model with indicator constraints.
# A driver that adds indicators as solver rows must still pass the dual
# start in the solver's row count (more rows than algebraic constraints).
# Each b[i] must be 1: with b[i] = 0, x[i] <= 2 cannot meet Demand[i].
# Optimum: x = 0, b = 1, objective 0.5 * n.
# A driver that reads past the dual vector fails reliably only under a
# memory checker (AddressSanitizer, Guard Malloc), not in a plain run.

param n := 3000;

var x {1..n} >= 0, <= 10;
var b {1..n} binary;

minimize Cost: sum {i in 1..n} (x[i] + 0.5 * b[i]);

subject to Demand {i in 1..n}: x[i] + 5 * b[i] >= 3;

subject to Cap {i in 1..n}: b[i] = 0 ==> x[i] <= 2;

# Primal and dual initial values, no basis statuses: the driver
# receives a primal-dual start.
let {i in 1..n} x[i] := 1;
let {i in 1..n} Demand[i] := 1;
