var x{1..4} integer >= 0 <= 1;

minimize o: 3*x[1] + 3*x[2] + 6*x[3] + 5*x[4];

s.t. c1: forall {i in 1..3} x[i] != x[i+1];
s.t. c2: (x[1] <= 0) or (x[2] < 1) or (x[3] == 1);
