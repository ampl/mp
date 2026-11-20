var x{1..3} integer >= -2 <= 3;
maximize o: (5)*x[1] + (2)*x[2] + (-1)*x[3];
s.t. c1: ((if x[2] > 1 then 1 else 0) + (if x[2] >= 1 then 1 else 0)) >= 1;
s.t. c2: alldiff(x[2], x[1], x[3]);
s.t. c3: max(x[1], x[2]) - min(x[2], x[1]) >= 3;
s.t. c4: numberof 0 in (x[1], x[2], x[3]) == 3;
