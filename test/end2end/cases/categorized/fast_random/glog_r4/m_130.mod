var x{1..6} integer >= 0 <= 4;

minimize o: 5*x[1] + 0*x[2] + 0*x[3] - 4*x[4] + 6*x[5] + 2*x[6];

s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: forall {i in 1..5} x[i] != x[i+1];
