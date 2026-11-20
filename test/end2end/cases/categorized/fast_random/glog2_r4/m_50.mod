var x{1..6} integer >= -2 <= 2;
minimize o: (5)*x[1] + (6)*x[2] + (2)*x[3] + (-2)*x[4] + (3)*x[5] + (1)*x[6];
s.t. c1: exists {i in 1..6} x[i] == -1;
s.t. c2: max(x[5], x[4]) - min(x[1], x[1]) > 1;
