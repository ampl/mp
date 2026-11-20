var x{1..3} integer >= 0 <= 3;
maximize o: (-3)*x[1] + (6)*x[2] + (6)*x[3];
s.t. c1: ((if x[2] < 2 then 1 else 0) + (if x[3] >= 2 then 1 else 0)) >= 1;
s.t. c2: alldiff(x[2], x[1], x[3]);
