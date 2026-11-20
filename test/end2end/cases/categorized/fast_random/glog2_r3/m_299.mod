var x{1..6} integer >= 0 <= 3;

maximize o: -5*x[1] - 5*x[2] + 5*x[3] + 5*x[4] - 3*x[5] + x[6];

s.t. c1: (x[3] <= 1) ==> (x[4] >= 0) else (x[4] < 1);
s.t. c2: exists {i in 1..6} x[i] == 3;
s.t. c3: (x[3] <= 1) ==> (x[4] >= 0) else (x[4] < 1);
s.t. c4: exists {i in 1..6} x[i] == 3;
