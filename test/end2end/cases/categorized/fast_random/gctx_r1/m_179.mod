var x{1..4} integer >= -1 <= 2;
minimize o: (-2)*x[1] + (5)*x[2] + (1)*x[3] + (6)*x[4];
s.t. c1: (if max(x[1], x[2]) >= 2 then x[4] else 0) + (if x[2] + x[1] != 0 then x[4] else 0) > 0;
s.t. c2: (x[2] + x[1] != 0) <==> (max(x[1], x[2]) >= 2);
