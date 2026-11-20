var x{1..6} integer >= 0 <= 2;

minimize o: -x[1] + 5*x[2] + 2*x[3] + 3*x[4] + 2*x[5] - 5*x[6];

s.t. c1: exists {i in 1..6} x[i] == 1;
s.t. c2: (x[6] < 2) and (x[5] < 0) ==> (x[1] >= 2);
