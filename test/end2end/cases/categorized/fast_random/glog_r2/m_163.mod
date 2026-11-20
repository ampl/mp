var x{1..5} integer >= 0 <= 3;
maximize o: (6)*x[1] + (0)*x[2] + (2)*x[3] + (-3)*x[4] + (4)*x[5];
s.t. c1: (x[2] >= 0) or (x[4] != 3) or (x[3] <= 0);
s.t. c2: alldiff(x[3], x[2], x[1], x[4]);
s.t. c3: abs(x[4] - x[4]) == 3;
s.t. c4: exists {i in 1..5} x[i] == 2;
