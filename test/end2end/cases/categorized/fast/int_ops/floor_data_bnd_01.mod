# -------------------------------------------------------------
# floor_data_bnd_01.mod
# floor(x) with the data bound x <= 3.9999995: the result bound is
# floor(3.9999995) = 3; only a rounding error of a computed bound
# (cvt:pre:ulps) may move it to 4, not a feasibility tolerance.
# Optimum 3.
# -------------------------------------------------------------
var x >= 0 <= 3.9999995;

maximize o: floor(x);
