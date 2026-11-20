var x{1..5} integer >= -2 <= 1;
maximize o: (4)*x[1] + (2)*x[2] + (6)*x[3] + (-2)*x[4] + (-3)*x[5];
s.t. c1: numberof 1 in (x[1], x[2], x[3], x[4], x[5]) <= 3;
s.t. c2: alldiff(x[2], x[5], x[3], x[4]);
s.t. c3: exists {i in 1..5} x[i] == -1;
