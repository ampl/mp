# power5.mod

param N integer, := 2;
set I := 1..N;

var x{I} >= 2.8;

maximize Sum:
     -5 * (x[1]-0.7)^2 + x[2]^7;

s.t. c1: 2 * x[1] + x[2] <= 10.2;
s.t. c2: (1.0 / 9) *
           (2.3*x[1] + 1.57*x[2] - 3.4)^5 +
            x[2]^2 >= 1;
s.t. c3: 8 * x[1]^2 + x[2] >= 0.5;

