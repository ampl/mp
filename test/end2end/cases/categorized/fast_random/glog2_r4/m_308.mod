var x{1..4} integer >= 0 <= 1;

minimize o: -3*x[1] + 6*x[2] + 4*x[3] + 4*x[4];

s.t. c1: exists {i in 1..4} x[i] == 1;
s.t. c2: (x[4] > 0) and (x[2] <= 0) ==> (x[4] >= 0);
