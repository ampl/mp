var x{1..5} integer >= 0 <= 4;
maximize o: (-1)*x[1] + (2)*x[2] + (2)*x[3] + (3)*x[4] + (6)*x[5];
s.t. c1: exists {i in 1..5} x[i] == 3;
s.t. c2: max(x[2], x[5]) - min(x[5], x[3]) == 2;
s.t. c3: alldiff(x[3], x[5]);
