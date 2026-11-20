var x{1..6} integer >= -2 <= 2;

minimize o: 4*x[1] + 3*x[2] - 4*x[3] + 5*x[4] - 3*x[5] + 0*x[6];

s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: atleast 2 {i in 1..6} (x[i] >= 0);
