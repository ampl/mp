var x{1..5} integer >= 0 <= 2;
maximize o: (6)*x[1] + (6)*x[2] + (4)*x[3] + (-3)*x[4] + (6)*x[5] + (-1)*(if x[4] > 2 then 1 else 0);
s.t. c1: not (x[5] >= 2);
s.t. c2: (if x[5] >= 2 then x[4] else 0) + (if x[4] > 2 then x[3] else 0) <= 2;
