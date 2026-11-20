var x{1..5} integer >= -1 <= 2;
minimize o: (2)*x[1] + (-1)*x[2] + (1)*x[3] + (6)*x[4] + (-2)*x[5] + (-3)*(if max(x[4], x[2]) == 0 then 1 else 0);
s.t. c1: (max(x[4], x[2]) == 0) <==> (x[1] <= -1);
s.t. c2: (if min(x[1], x[3]) >= -1 then x[3] else 0) + (if min(x[1], x[3]) >= -1 then x[4] else 0) < 0;
s.t. c3: (min(x[1], x[3]) >= -1) or (min(x[1], x[3]) >= -1);
s.t. c4: (if x[1] <= -1 then 1 else 0) + (if min(x[3], x[3]) != 1 then 1 else 0) + (if x[1] <= -1 then 1 else 0) >= 1;
