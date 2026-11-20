var x{1..3} integer >= -2 <= 1;
maximize o: (5)*x[1] + (1)*x[2] + (5)*x[3];
s.t. c1: alldiff(x[3], x[1]);
s.t. c2: (x[3] == -2) ==> (x[3] < -3) else (x[1] <= -3);
s.t. c3: max(x[3], x[1]) - min(x[2], x[3]) == -2;
s.t. c4: ((if x[2] == -2 then 1 else 0) + (if x[3] >= -2 then 1 else 0)) >= 1;
