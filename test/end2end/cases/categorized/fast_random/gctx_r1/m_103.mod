var x{1..4} integer >= 0 <= 2;
minimize o: (-3)*x[1] + (-4)*x[2] + (-1)*x[3] + (6)*x[4];
s.t. c1: (min(x[1], x[2]) == 2) ==> (abs(x[3] - x[3]) < 0) else (min(x[1], x[2]) == 2);
s.t. c2: (abs(x[3] - x[3]) < 0) and ((min(x[1], x[2]) == 2) or not (min(x[1], x[2]) == 2));
s.t. c3: (if x[1] < x[1] then x[3] else 0) + (if abs(x[4] - x[1]) < 0 then x[3] else 0) > 0;
