##########################################
## unify_lin_expr_01.mod
## Test that equivalent linear expressions are unified
##########################################

var x {1..4} >=-3, <=5.27 integer;

var dx = 3*x[2] + 2*x[3];

minimize Obj: count (-5*x[4] + dx == 3,
                      2*x[3] - 5*x[4] + 3*x[2] == 8,
                      6*x[1] + 2*x[2] == 0);

s.t. C1: alldiff (dx - 5*x[4], 3*x[4] - 7*x[1], 17*x[3]);

s.t. C2: 2 <= numberof (2*x[3] - 5*x[4] + 3*x[2])
                in (-5*x[4] + dx, 3*x[4] - 7*x[1], 17*x[3]);