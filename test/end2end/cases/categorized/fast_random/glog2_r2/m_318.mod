var x{1..6} integer >= 0 <= 1;

minimize o: 5*x[1] + 0*x[2] - 2*x[3] + 5*x[4] + 2*x[5] + 5*x[6];

s.t. c1: (if x[2] >= 1 then x[4] else x[4]) > 1;
s.t. c2: exists {i in 1..6} x[i] == 1;
