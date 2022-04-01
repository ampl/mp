# div01.mod

param N integer, := 2;
set I := 1..N;

var x{I} >= -2.8;

maximize Div:
     -5 * (x[1]-0.7) / ((x[2]+3)^2);

s.t. c1: 2 * x[1] + x[2] <= 10.2;
s.t. c3: 8 * x[1]^2 + x[2] >= 0.5;

