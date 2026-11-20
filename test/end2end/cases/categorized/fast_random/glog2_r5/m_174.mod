var x{1..4} integer >= -3 <= 1;

minimize o: -3*x[1] + 2*x[2] + 0*x[3] - 3*x[4];

s.t. c1: forall {i in 1..3} x[i] != x[i+1];
