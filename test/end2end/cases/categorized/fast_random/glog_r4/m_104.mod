var x{1..5} integer >= 0 <= 4;
maximize o: (0)*x[1] + (3)*x[2] + (6)*x[3] + (-1)*x[4] + (-5)*x[5];
s.t. c1: (x[1] < 0) ==> (x[1] < 3) else (x[5] <= 1);
s.t. c2: max(x[1], x[3]) - min(x[5], x[4]) <= 4;
s.t. c3: alldiff(x[4], x[3], x[2]);
