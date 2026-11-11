# -------------------------------------------------------------
# asinh_02_neg_arg.mod
# asinh with a negative argument at the optimum (x = -1.7456):
# asinh is odd, asinh(x) = -asinh(|x|).
# Unique stationary point of (x+1)^2 + 3*asinh(x) on [-3, 4].
# -------------------------------------------------------------
var x >= -3, <= 4, := -1;
var y;
minimize Obj: (x + 1)^2 + 3*y;
s.t. Asinh: y == asinh(x);
