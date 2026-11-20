var x{1..6} integer >= -2 <= 1;

maximize o: 3*x[1] - 2*x[2] + 3*x[3] + 0*x[4] + 5*x[5] - x[6];

s.t. c1: exists {i in 1..6} x[i] == -1;
s.t. c2: (x[3] < -2) ==> (x[4] <= -2) else (x[3] != 0);
s.t. c3: (x[6] >= -1) and (x[2] != 0) ==> (x[2] >= -1);
