# -------------------------------------------------------------
# chk_ctx_count_01.mod
# sol:chk:mode=64 compares the solver's values of expressions
# with recomputed ones on the side their context requires.
# x >= 1 becomes not (x == 0); the binary of x == 0 is bounded
# from one side only, and a solver may return 1 at x = 4.
# The count is then 2 from the solver's values and 3 recomputed:
# in positive context (count >= 2) that is no violation.
# Optimum -4 at x = 4.
# -------------------------------------------------------------
var x integer >= 0 <= 4;

minimize o: -x;

s.t. c: atleast 2 (x >= 2, x >= 1, x >= 3);
