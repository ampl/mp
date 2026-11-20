var x{1..5} integer >= -1 <= 2;
maximize o: (4)*x[1] + (-2)*x[2] + (3)*x[3] + (0)*x[4] + (-4)*x[5] + (2)*(if x[5] <= 1 then 1 else 0);
s.t. c1: (min(x[5], x[1]) < -1) <==> (x[5] <= 1);
s.t. c2: not ((x[5] <= 1) and (x[5] <= 1));
s.t. c3: (if x[5] <= 1 then 1 else 0) + (if min(x[5], x[1]) < -1 then 1 else 0) + (if x[5] <= 1 then 1 else 0) >= 2;
