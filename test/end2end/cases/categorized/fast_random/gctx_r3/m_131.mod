var x{1..5} integer >= 0 <= 3;
maximize o: (1)*x[1] + (2)*x[2] + (6)*x[3] + (6)*x[4] + (5)*x[5];
s.t. c1: not (x[5] >= 2);
s.t. c2: (if x[5] >= 2 then x[1] else 0) + (if x[1] >= x[3] then x[4] else 0) >= 0;
