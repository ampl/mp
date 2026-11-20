var x{1..4} integer >= 0 <= 4;
maximize o: (6)*x[1] + (-4)*x[2] + (-4)*x[3] + (-3)*x[4];
s.t. c1: alldiff(x[4], x[2]);
s.t. c2: (x[4] <= 3) ==> (x[1] <= 2) else (x[2] < 3);
s.t. c3: alldiff(x[4], x[3], x[1], x[2]);
