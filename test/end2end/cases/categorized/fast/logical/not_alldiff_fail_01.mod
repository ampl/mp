# -------------------------------------------------------------
# not_alldiff_fail_01.mod
# The negation of alldiff is not supported: the driver must report
# a failure (solve_result_num 500), not a solved status.
# -------------------------------------------------------------
var x integer >= 0 <= 4;
var y integer >= 0 <= 4;
var z integer >= 0 <= 1;
var b binary;

maximize Obj: x + y + z - 10*b;

s.t. C: not alldiff(x, y, z) or b = 1;
