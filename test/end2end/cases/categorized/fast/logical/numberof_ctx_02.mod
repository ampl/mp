# -------------------------------------------------------------
# numberof_ctx_02.mod
# As numberof_ctx_01, with max() as a compared argument: it needs
# both directions, which failed a conversion before
# ("Part of expression redefinition could be lost").
# Optimum 3, e.g., at x = 3, y = 0.
# -------------------------------------------------------------
var x integer >= 0 <= 5;
var y integer >= 0 <= 5;

minimize o: x + y;

s.t. c: numberof 3 in (max(x, y), x + y) >= 1;
