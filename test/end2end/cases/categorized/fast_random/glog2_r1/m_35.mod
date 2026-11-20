var x{1..3} integer >= 0 <= 2;
maximize o: (5)*x[1] + (2)*x[2] + (6)*x[3];
s.t. c1: max(x[2], x[1]) - min(x[1], x[3]) < 1;
s.t. c2: exists {i in 1..3} x[i] == 1;
s.t. c3: (x[1] >= 0) ==> (x[2] == 2);
s.t. c4: not ((x[1] > 0) and (x[3] != 0));
s.t. c5: alldiff(x[1], x[2]);
s.t. c6: abs(x[1] - x[3]) < 1;
s.t. c7: abs(x[1] - x[3]) <= 1;
