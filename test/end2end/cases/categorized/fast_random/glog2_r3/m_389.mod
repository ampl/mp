var x{1..6} integer >= 0 <= 2;

maximize o: 6*x[1] + 2*x[2] - x[3] + 4*x[4] - 2*x[5] + 4*x[6];

s.t. c1: exists {i in 1..6} x[i] == 1;
s.t. c2: (if x[3] > 1 then x[1] else x[4]) < 1;
s.t. c3: exists {i in 1..6} x[i] == 1;
s.t. c4: (if x[3] > 1 then x[1] else x[4]) < 1;
