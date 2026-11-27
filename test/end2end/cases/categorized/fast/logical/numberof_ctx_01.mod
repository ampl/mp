# -------------------------------------------------------------
# numberof_ctx_01.mod
# The arguments of numberof are compared for equality, so they
# are needed in both directions (mixed context), whatever the
# context of numberof itself. Here the inner numberof must be 1
# when x = 2, so c requires x != 2 (and y = 0).
# Optimum -1 at x = 1 or x = 3, y = 0.
# -------------------------------------------------------------
var x integer >= 0 <= 3;
var y integer >= -1 <= 1;

maximize o: -abs(x - 2);

# both the inner numberof and y must be 0
s.t. c: numberof 0 in (numberof 2 in (x, y), y) >= 2;
