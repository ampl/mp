var x{1..4} integer >= 0 <= 3;
maximize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4];
s.t. c1: alldiff(x[3], x[1]);
s.t. c2: (x[4] != 0) ==> (x[2] < 1);
s.t. c3: (x[2] < 0) and (x[2] >= 1) ==> (x[2] >= 0);
