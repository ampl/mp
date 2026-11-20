var x{1..6} integer >= -2 <= 2;

minimize o: x[1] - 2*x[2] + 2*x[3] - 5*x[4] + 5*x[5] + 6*x[6];

s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: forall {i in 1..5} x[i] != x[i+1];
