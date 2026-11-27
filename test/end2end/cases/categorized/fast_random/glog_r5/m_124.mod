var x{1..4} integer >= -2 <= 0;

minimize o: -5*x[1] + 0*x[2] + 4*x[3] - 2*x[4];

s.t. c1: exists {i in 1..4} x[i] == -1;
s.t. c2: forall {i in 1..3} x[i] != x[i+1];
