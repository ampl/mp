var x{1..4} integer >= -1 <= 1;
minimize o: (4)*x[1] + (0)*x[2] + (6)*x[3] + (-3)*x[4] + (-4)*(if x[1] <= x[4] then 1 else 0);
s.t. c1: (abs(x[4] - x[2]) == 2) <==> (max(x[2], x[4]) <= 1);
s.t. c2: (max(x[2], x[4]) <= 1) ==> (abs(x[2] - x[3]) == 1) else (max(x[2], x[4]) <= 1);
s.t. c3: (if abs(x[4] - x[2]) == 2 then x[3] else 0) + (if abs(x[2] - x[3]) == 1 then x[3] else 0) > -1;
