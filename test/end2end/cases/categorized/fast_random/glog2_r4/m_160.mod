var x{1..6} integer >= 0 <= 2;

minimize o: 2*x[1] + 4*x[2] + 6*x[3] - 2*x[4] + 5*x[5] - x[6];

s.t. c1: (x[4] < 2) and (x[6] != 0) ==> (x[1] >= 2);
s.t. c2: exists {i in 1..6} x[i] == 1;
