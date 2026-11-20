var x{1..4} integer >= -1 <= 2;
maximize o: (6)*x[1] + (5)*x[2] + (4)*x[3] + (5)*x[4] + (-5)*(if min(x[4], x[1]) <= -1 then 1 else 0);
s.t. c1: (if x[2] >= 0 then 1 else 0) + (if x[2] >= 0 then 1 else 0) + (if min(x[4], x[1]) <= -1 then 1 else 0) >= 2;
s.t. c2: not ((min(x[4], x[1]) <= -1) and (x[2] >= 0));
