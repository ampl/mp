# -------------------------------------------------------------
# funcbnd_log_03_zero.mod
# log over an argument whose bounds become [0, 0]: the result range
# is degenerate (both ends -infinity). Such a range is not used:
# HiGHS stopped with an error in Highs_addRows; now infeasible.
# -------------------------------------------------------------
var x >= -1, <= 0;
var y >= -5, <= 5;
minimize o: y;
s.t. c: y >= log(x) + 1;
