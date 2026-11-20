var x{1..6} integer >= -2 <= 2;

maximize o: 0*x[1] + x[2] + 3*x[3] - 5*x[4] - 5*x[5] + 6*x[6];

s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: exists {i in 1..6} x[i] == -1;
