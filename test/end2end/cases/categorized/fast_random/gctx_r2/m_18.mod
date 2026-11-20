var x{1..4} integer >= 0 <= 3;
minimize o: (-5)*x[1] + (6)*x[2] + (-2)*x[3] + (-2)*x[4] + (0)*(if abs(x[1] - x[2]) <= 3 then 1 else 0);
s.t. c1: (if x[4] >= 0 then x[1] else 0) + (if x[4] >= 0 then x[4] else 0) >= 1;
s.t. c2: not ((abs(x[1] - x[2]) <= 3) and (x[4] >= 0));
s.t. c3: not ((abs(x[1] - x[2]) <= 3) and (x[4] >= 0));
