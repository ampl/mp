var x{1..6} integer >= 0 <= 3;

maximize o: -3*x[1] + 6*x[2] - 2*x[3] - 3*x[4] - 3*x[5] - 5*x[6];

s.t. c1: (if x[2] <= 2 then x[3] else x[3]) > 0;
s.t. c2: exists {i in 1..6} x[i] == 3;
