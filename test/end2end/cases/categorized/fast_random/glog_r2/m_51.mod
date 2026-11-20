var x{1..3} integer >= 0 <= 3;
maximize o: (1)*x[1] + (0)*x[2] + (-3)*x[3];
s.t. c1: (x[1] < 0) and (x[3] == 3) ==> (x[2] >= 0);
s.t. c2: alldiff(x[2], x[3]);
s.t. c3: alldiff(x[1], x[3]);
