var x{1..6} integer >= 0 <= 3;
maximize o: (2)*x[1] + (6)*x[2] + (4)*x[3] + (2)*x[4] + (-2)*x[5] + (3)*x[6];
s.t. c1: exists {i in 1..6} x[i] == 1;
s.t. c2: not ((x[3] > 3) and (x[4] >= 1));
