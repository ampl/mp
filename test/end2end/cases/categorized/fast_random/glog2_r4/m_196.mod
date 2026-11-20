var x{1..3} integer >= 0 <= 4;
maximize o: (0)*x[1] + (0)*x[2] + (0)*x[3];
s.t. c1: (x[3] <= 3) and (x[3] > 3) ==> (x[3] >= 3);
s.t. c2: not ((x[3] == 0) and (x[1] != 2));
s.t. c3: alldiff(x[1], x[3]);
