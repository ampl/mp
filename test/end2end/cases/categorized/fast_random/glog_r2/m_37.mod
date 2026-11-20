var x{1..4} integer >= 0 <= 4;
maximize o: (-3)*x[1] + (4)*x[2] + (5)*x[3] + (6)*x[4];
s.t. c1: (x[3] != 1) ==> (x[2] < 2) else (x[1] <= 0);
s.t. c2: alldiff(x[1], x[2]);
