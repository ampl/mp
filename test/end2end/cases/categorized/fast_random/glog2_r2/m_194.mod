var x{1..5} integer >= -2 <= 2;
maximize o: (5)*x[1] + (5)*x[2] + (0)*x[3] + (3)*x[4] + (-2)*x[5];
s.t. c1: not ((x[5] <= -1) and (x[2] > 0));
s.t. c2: alldiff(x[5], x[4], x[2], x[3]);
s.t. c3: exists {i in 1..5} x[i] == 2;
