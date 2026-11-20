var x{1..5} integer >= 0 <= 3;
minimize o: (0)*x[1] + (-2)*x[2] + (-4)*x[3] + (6)*x[4] + (4)*x[5];
s.t. c1: not ((max(x[3], x[3]) >= 2) and (min(x[1], x[5]) > 3));
s.t. c2: (max(x[3], x[3]) >= 2) and ((min(x[1], x[5]) > 3) or not (max(x[3], x[3]) >= 2));
s.t. c3: (if max(x[3], x[3]) >= 2 then 1 else 0) + (if max(x[3], x[3]) >= 2 then 1 else 0) + (if max(x[3], x[3]) >= 2 then 1 else 0) >= 1;
