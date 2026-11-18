# -------------------------------------------------------------
# exp_square_ctx_01.mod
# Generated: log(1 + exp(.)^2) in a >= constraint. The exp result
# was relaxed to v <= exp(.) using its best-known bound v >= 0, which
# was not submitted: the solver took v very negative (Gurobi 14:
# infeasible point, -11.5). Optimum -8.7 by brute force.
# -------------------------------------------------------------
var x{1..4} integer >= -2, <= 2;
minimize o: (3*x[3] + 1*x[4] + 2*x[2] + -2*x[1]) + 0.1*(2*(3*x[3] + 3*x[1] + -1*x[2] + 1) - (x[1]*x[4] + (2*x[2] + 1)));
s.t. c0: (3*x[3] + -2*x[1] + 3*x[2]) + log(1 + (exp((3*x[4] + 2*x[3] + 1*x[1] + -1)/8))^2) >= -5;
s.t. c1: (not ((3*x[2] + -2*x[4] + 2) <= 1)) ==> (exp((2*x[2] + 2*x[1] + 1*x[4] + 2)/8) <= 3);
s.t. c2: ((1*x[3] + 3*x[4] + -1) <= 1) <==> not ((3*x[2] + -1*x[4] + 2) <= 2);
