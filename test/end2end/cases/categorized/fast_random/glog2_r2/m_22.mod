var x{1..5} integer >= 0 <= 4;
maximize o: (1)*x[1] + (-3)*x[2] + (6)*x[3] + (3)*x[4] + (2)*x[5];
s.t. c1: alldiff(x[4], x[3], x[2], x[5]);
s.t. c2: max(x[5], x[1]) - min(x[1], x[4]) != 3;
s.t. c3: numberof 0 in (x[1], x[2], x[3], x[4], x[5]) >= 3;
s.t. c4: exists {i in 1..5} x[i] == 1;
