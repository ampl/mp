var x{1..3} integer >= 0 <= 2;
maximize o: (5)*x[1] + (-1)*x[2] + (-1)*x[3];
s.t. c1: alldiff(x[2], x[3], x[1]);
s.t. c2: ((if x[1] < 0 then 1 else 0) + (if x[3] >= 0 then 1 else 0)) >= 1;
s.t. c3: max(x[1], x[3]) - min(x[1], x[2]) == 2;
