var x{1..5} integer >= 0 <= 2;
maximize o: (5)*x[1] + (3)*x[2] + (1)*x[3] + (4)*x[4] + (6)*x[5];
s.t. c1: (if x[2] <= 0 then x[1] else x[2]) == 3;
s.t. c2: alldiff(x[5], x[3], x[1]);
