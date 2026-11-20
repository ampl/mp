var x{1..4} integer >= 0 <= 3;
maximize o: (-2)*x[1] + (3)*x[2] + (0)*x[3] + (-1)*x[4];
s.t. c1: ((if x[3] == 0 then 1 else 0) + (if x[4] >= 0 then 1 else 0)) >= 1;
s.t. c2: alldiff(x[4], x[2]);
s.t. c3: ((if x[3] == 0 then 1 else 0) + (if x[4] >= 0 then 1 else 0)) >= 1;
s.t. c4: alldiff(x[4], x[2]);
