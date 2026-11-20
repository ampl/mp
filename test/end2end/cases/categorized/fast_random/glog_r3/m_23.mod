var x{1..4} integer >= 0 <= 2;
maximize o: (3)*x[1] + (1)*x[2] + (4)*x[3] + (5)*x[4];
s.t. c1: exists {i in 1..4} x[i] == 0;
s.t. c2: alldiff(x[2], x[4]);
s.t. c3: alldiff(x[1], x[2], x[3], x[4]);
s.t. c4: max(x[2], x[2]) - min(x[2], x[3]) < 0;
