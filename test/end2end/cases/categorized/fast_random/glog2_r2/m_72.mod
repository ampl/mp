var x{1..5} integer >= 0 <= 2;
maximize o: (2)*x[1] + (4)*x[2] + (6)*x[3] + (1)*x[4] + (0)*x[5];
s.t. c1: numberof 2 in (x[1], x[2], x[3], x[4], x[5]) <= 3;
s.t. c2: (x[4] >= 0) or (x[4] != 1) or (x[4] <= 1);
s.t. c3: exists {i in 1..5} x[i] == 2;
