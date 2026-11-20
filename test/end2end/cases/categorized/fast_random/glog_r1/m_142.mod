var x{1..6} integer >= -2 <= 0;

minimize o: 4*x[1] - 5*x[2] - 5*x[3] + 0*x[4] + 2*x[5] + 3*x[6];

s.t. c1: exists {i in 1..6} x[i] == -2;
s.t. c2: max(x[4], x[6]) - min(x[4], x[6]) > -2;
s.t. c3: (x[2] > 0) ==> (x[4] > 0);
s.t. c4: (x[1] == -2) ==> (x[6] <= 0);
