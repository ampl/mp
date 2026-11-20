var x{1..6} integer >= -2 <= 2;
maximize o: (5)*x[1] + (6)*x[2] + (6)*x[3] + (-3)*x[4] + (-3)*x[5] + (6)*x[6];
s.t. c1: exists {i in 1..6} x[i] == -2;
s.t. c2: max(x[6], x[6]) - min(x[6], x[6]) >= -2;
s.t. c3: not ((x[5] <= 2) and (x[5] >= 1));
