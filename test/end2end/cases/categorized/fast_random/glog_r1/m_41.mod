var x{1..5} integer >= 0 <= 2;

minimize o: -4*x[1] + 5*x[2] + 6*x[3] - 3*x[4] - 4*x[5];

s.t. c1: forall {i in 1..4} x[i] != x[i+1];
