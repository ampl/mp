# -------------------------------------------------------------
# pow_08_frac_obj_con.mod
# Powers with constant non-integer exponents (x^1.5, y^2.5,
# y^1.25) in the objective and in an active constraint.
# Convex; checks values, gradients and Hessians of x^c, c != 2.
# -------------------------------------------------------------
var x >= 0.1, <= 10, := 1;
var y >= 0.1, <= 10, := 1;
minimize Obj: x^1.5 + y^2.5 - 3*x - 5*y;
s.t. Cap: x^1.5 + y^1.25 <= 9;
