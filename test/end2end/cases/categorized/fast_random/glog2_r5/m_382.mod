var x{1..4} integer >= 0 <= 4;
maximize o: (6)*x[1] + (2)*x[2] + (6)*x[3] + (1)*x[4];
s.t. c1: ((if x[3] == 3 then 1 else 0) + (if x[4] >= 3 then 1 else 0)) >= 1;
s.t. c2: max(x[1], x[4]) - min(x[1], x[1]) != 0;
s.t. c3: not ((x[4] < 1) and (x[4] >= 3));
s.t. c4: (x[4] > 2) <==> (x[4] + x[1] <= 2);
