# -------------------------------------------------------------
# round_bnd_fp_01.mod
# Result bounds of ceil() from floating-point argument bounds:
# the lower bound of (-x - 1)/3 is computed as -1.9999999999999998,
# but ceil(-2) = -2 is attained at x = 5. Optimum -2 at x = 5.
# -------------------------------------------------------------
var x integer >= 0 <= 5;

minimize o: ceil((-x - 1) / 3);
