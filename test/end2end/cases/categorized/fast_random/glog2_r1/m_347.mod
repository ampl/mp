var x{1..3} integer >= 0 <= 2;
maximize o: (-2000)*x[1] + (3000)*x[2] + (-999)*x[3];
s.t. c1: alldiff(x[3], x[1], x[2]);
s.t. c2: alldiff(x[3], x[2]);
s.t. c3: ((if x[1] != 1 then 1 else 0) + (if x[1] >= 1 then 1 else 0)) >= 1;
