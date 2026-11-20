var x{1..5} integer >= 0 <= 2;
minimize o: (0)*x[1] + (6)*x[2] + (-2)*x[3] + (-3)*x[4] + (3)*x[5];
s.t. c1: (x[5] != x[4]) <==> (min(x[2], x[4]) > 2);
s.t. c2: (if x[4] + x[2] > 2 then 1 else 0) + (if x[5] != x[4] then 1 else 0) + (if min(x[2], x[4]) > 2 then 1 else 0) >= 2;
