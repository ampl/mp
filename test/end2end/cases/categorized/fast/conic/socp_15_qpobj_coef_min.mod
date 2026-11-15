# As socp_14_qpobj_coef_max, minimization (convex quadratic objective).
# Optimum: x = 1.875, y = 0, objective -2.8125.

var x >= 0 <= 10;
var y >= 0 <= 10;
minimize qobj: 0.8*x^2 - 3*x + y;
subj to qconstr: x^2 + y^2 <= 100;
