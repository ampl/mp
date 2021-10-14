# ellipse_nonconvex.mod
# Adding a DStart value

# Number of variables: 2
# Number of constraints:  2
# Objective linear
# Quadratic constraints

param N integer, := 2;
set I := 1..N;

param a >= 0;
param b >= 0;
param c >= 0;
param d >= 0;
param e >= 0;
param f >= 0;
param g >= 0;
param h >= 0;

var x{I};

maximize Sum:
     -5 * (abs(x[1])-0.7)^2 + x[2];

s.t. c1: 2 * x[1] + x[2] <= 2;
s.t. c2: (1.0 / 9) * x[1]^2 + x[2]^2 >= 1;
s.t. c3: 8 * x[1]^2 + x[2] >= 0.5;


data;
param a := 0.0025;
param b := 0.01;
param c := 833.3325;
param d := 100;
param e := 83333.33;
param f := 1250;
param g := 1250000;
param h := 2500;

var x :=
    1  0.5   2  1;

    
let _con[1] := 0;
