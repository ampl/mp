var x{1..6} integer >= 0 <= 4;

maximize o: 5*x[1] + 4*x[2] + 3*x[3] - x[4] + 4*x[5] - 4*x[6];

s.t. c1: exists {i in 1..6} x[i] == 2;
s.t. c2: (x[3] <= 1) ==> (x[2] < 0) else (x[6] == 4);
