var x{1..3} integer >= -1 <= 1;
maximize o: (2)*x[1] + (6)*x[2] + (6)*x[3];
s.t. c1: not ((x[3] > x[3]) and (x[3] > x[3]));
s.t. c2: (if x[3] <= 1 then x[2] else 0) + (if x[3] > x[3] then x[2] else 0) < -1;
s.t. c3: not (x[3] > x[3]);
