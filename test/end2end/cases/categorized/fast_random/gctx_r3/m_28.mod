var x{1..3} integer >= 0 <= 3;
minimize o: (-4)*x[1] + (-3)*x[2] + (6)*x[3] + (1)*(if abs(x[2] - x[3]) >= 2 then 1 else 0);
s.t. c1: (max(x[3], x[1]) == 2) <==> (max(x[3], x[1]) == 2);
s.t. c2: (if x[1] > x[3] then x[1] else 0) + (if x[1] > x[3] then x[3] else 0) >= 0;
