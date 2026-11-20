var x{1..6} integer >= 0 <= 3;

minimize o: -3*x[1] - 2*x[2] - 2*x[3] - 4*x[4] + 4*x[5] + 5*x[6];

s.t. c1: exists {i in 1..6} x[i] == 2;
s.t. c2: (x[6] >= 3) ==> (x[4] == 1) else (x[4] == 2);
