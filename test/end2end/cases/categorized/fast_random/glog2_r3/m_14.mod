var x{1..5} integer >= -2 <= 3;

minimize o: 6*x[1] + 2*x[2] + 4*x[3] + x[4] + x[5];

s.t. c1: forall {i in 1..4} x[i] != x[i+1];
