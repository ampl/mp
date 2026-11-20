var x{1..6} integer >= 0 <= 4;

maximize o: -4*x[1] - 4*x[2] - 2*x[3] - 4*x[4] - 5*x[5] - 5*x[6];

s.t. c1: exists {i in 1..6} x[i] == 3;
s.t. c2: (x[3] == 4) ==> (x[1] >= 0) else (x[4] > 2);
