var x{1..5} integer >= 0 <= 1;
maximize o: (-3)*x[1] + (2)*x[2] + (-3)*x[3] + (2)*x[4] + (6)*x[5];
s.t. c1: ((if x[3] == 1 then 1 else 0) + (if x[2] >= 1 then 1 else 0)) >= 1;
s.t. c2: max(x[1], x[2]) - min(x[1], x[4]) <= 1;
s.t. c3: max(x[2], x[4]) - min(x[4], x[3]) > 0;
s.t. c4: alldiff(x[4], x[1]);
