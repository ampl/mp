var x{1..4} integer >= 0 <= 3;
minimize o: (6)*x[1] + (-5)*x[2] + (-2)*x[3] + (-4)*x[4];
s.t. c1: (if abs(x[2] - x[3]) >= 1 then x[3] else 0) + (if max(x[2], x[2]) <= 0 then x[2] else 0) < 2;
s.t. c2: (x[3] == 1) <==> (x[3] == 1);
