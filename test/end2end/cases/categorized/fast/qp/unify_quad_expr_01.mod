######################################
## unify_quad_expr_01.mod
## Test that equibvalent QP expressions are unified
######################################

var x {1..3} >=-1 <=1;

var dx = x[2]*x[3] + 1.5*x[2];

maximize Obj: exp(x[2]*x[1] + x[1] + dx);

s.t. C1: cos(1.5*x[2] + x[1]*x[2] + x[1] + x[3]*x[2])
         - sin(dx + x[1]*x[2] + x[1]) <= 0.5;
s.t. C2: x[1]*x[2] + x[1] + x[3] <= 1.5;
