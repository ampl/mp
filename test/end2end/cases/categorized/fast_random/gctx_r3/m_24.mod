var x{1..4} integer >= 0 <= 2;
maximize o: (6)*x[1] + (0)*x[2] + (-5)*x[3] + (2)*x[4] + (-3)*(if x[2] == x[4] then 1 else 0);
s.t. c1: (if min(x[2], x[1]) <= 1 then x[1] else 0) + (if min(x[2], x[1]) <= 1 then x[1] else 0) <= 1;
s.t. c2: (if abs(x[3] - x[2]) < 2 then x[3] else 0) + (if abs(x[3] - x[2]) < 2 then x[1] else 0) < 0;
