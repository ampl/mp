var x{1..6} integer >= 0 <= 4;

maximize o: 4*x[1] + 4*x[2] - 4*x[3] + 5*x[4] + 2*x[5] - 3*x[6];

s.t. c1: exists {i in 1..6} x[i] == 2;
s.t. c2: (x[2] < 1) ==> (x[6] == 2) else (x[2] > 2);
