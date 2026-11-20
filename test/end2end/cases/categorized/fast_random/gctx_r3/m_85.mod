var x{1..4} integer >= 0 <= 2;
maximize o: (0)*x[1] + (1)*x[2] + (6)*x[3] + (-3)*x[4];
s.t. c1: (if x[4] + x[3] <= 2 then x[3] else 0) + (if x[2] <= x[2] then x[3] else 0) > 2;
s.t. c2: (x[4] + x[3] <= 2) or (x[2] <= x[2]);
s.t. c3: (x[4] + x[3] <= 2) or (x[2] <= x[2]);
