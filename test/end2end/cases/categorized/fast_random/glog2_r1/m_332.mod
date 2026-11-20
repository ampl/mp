var x{1..4} integer >= 0 <= 2;
maximize o: (2)*x[1] + (-1)*x[2] + (-1)*x[3] + (-5)*x[4];
s.t. c1: (x[2] < 1) ==> (x[4] == 0) else (x[4] != 2);
s.t. c2: alldiff(x[4], x[1], x[2], x[3]);
