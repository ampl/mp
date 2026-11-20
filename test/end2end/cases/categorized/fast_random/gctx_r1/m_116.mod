var x{1..3} integer >= -1 <= 2;
maximize o: (-1)*x[1] + (4)*x[2] + (6)*x[3];
s.t. c1: not ((x[3] + x[2] <= 2) and (x[1] <= x[3]));
s.t. c2: (if x[1] <= x[3] then 1 else 0) + (if x[3] + x[2] <= 2 then 1 else 0) + (if x[1] <= x[3] then 1 else 0) >= 2;
