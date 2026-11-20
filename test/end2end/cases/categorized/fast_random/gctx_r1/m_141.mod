var x{1..3} integer >= 0 <= 2;
minimize o: (-3)*x[1] + (0)*x[2] + (-1)*x[3];
s.t. c1: (if x[3] > x[1] then 1 else 0) + (if max(x[1], x[3]) > 0 then 1 else 0) + (if max(x[1], x[3]) > 0 then 1 else 0) >= 2;
s.t. c2: (if max(x[1], x[3]) > 0 then x[1] else 0) + (if max(x[1], x[3]) > 0 then x[2] else 0) < 0;
s.t. c3: not ((x[3] > x[1]) and (x[3] > x[1]));
s.t. c4: (max(x[1], x[3]) > 0) <==> (max(x[1], x[3]) > 0);
