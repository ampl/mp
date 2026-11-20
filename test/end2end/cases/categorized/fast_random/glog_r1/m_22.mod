var x{1..6} integer >= 0 <= 4;

minimize o: -x[1] + 3*x[2] + 5*x[3] + 6*x[4] + 5*x[5] - 4*x[6];

s.t. c1: exists {i in 1..6} x[i] == 1;
s.t. c2: (x[4] <= 4) and (x[2] >= 1) ==> (x[3] >= 4);
