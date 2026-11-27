# -------------------------------------------------------------
# chk_ctx_max_01.mod
# As chk_ctx_count_01 for max(x, y) used below an upper bound
# (negative context): its result only needs to be >= x and >= y,
# and a solver may return 1 while the recomputed value is 0.
# Optimum 4 at x = y = 0, u = 4.
# -------------------------------------------------------------
var x >= 0 <= 3;
var y >= 0 <= 3;
var u >= 0 <= 9;

maximize o: u - x - y;

s.t. c: max(x, y) + u <= 5;
s.t. d: u <= 4;
