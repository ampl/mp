var x{1..5} integer >= 0 <= 2;
maximize o: (-666)*x[1] + (-999)*x[2] + (-666)*x[3] + (-999)*x[4] + (-1000)*x[5];
s.t. c1: (x[1] < 1) <==> (x[4] + x[2] == 0);
s.t. c2: ((if x[5] < 0 then 1 else 0) + (if x[3] >= 0 then 1 else 0)) >= 1;
s.t. c3: alldiff(x[4], x[3]);
