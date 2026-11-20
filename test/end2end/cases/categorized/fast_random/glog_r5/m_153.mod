var x{1..4} integer >= -2 <= 2;
maximize o: (5)*x[1] + (1)*x[2] + (5)*x[3] + (0)*x[4];
s.t. c1: alldiff(x[4], x[1], x[2], x[3]);
s.t. c2: (x[4] != 0) ==> (x[4] == 1) else (x[4] < 0);
s.t. c3: alldiff(x[4], x[2], x[1]);
s.t. c4: not ((x[3] == 2) and (x[1] > -1));
