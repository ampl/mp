var x{1..4} integer >= -2 <= 0;
minimize o: (3)*x[1] + (4)*x[2] + (5)*x[3] + (6)*x[4];
s.t. c1: ((if x[1] >= -2 then 1 else 0) + (if x[1] >= -2 then 1 else 0)) >= 1;
s.t. c2: max(x[3], x[1]) - min(x[3], x[3]) <= -2;
