var x{1..6} integer >= -2 <= 1;

maximize o: x[1] - 4*x[2] - x[3] + 4*x[4] - 4*x[5] - 5*x[6];

s.t. c1: forall {i in 1..5} x[i] != x[i+1];
