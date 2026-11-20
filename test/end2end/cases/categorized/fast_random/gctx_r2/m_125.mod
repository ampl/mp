var x{1..4} integer >= -1 <= 2;
minimize o: (6)*x[1] + (6)*x[2] + (0)*x[3] + (-1)*x[4] + (-1)*(if min(x[4], x[2]) <= 1 then 1 else 0);
s.t. c1: not (min(x[4], x[2]) <= 1);
s.t. c2: (if x[4] == 0 then 1 else 0) + (if x[4] == 0 then 1 else 0) + (if min(x[4], x[2]) <= 1 then 1 else 0) >= 2;
