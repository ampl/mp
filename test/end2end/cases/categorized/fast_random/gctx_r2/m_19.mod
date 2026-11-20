var x{1..5} integer >= -1 <= 1;
maximize o: (-4)*x[1] + (2)*x[2] + (6)*x[3] + (5)*x[4] + (6)*x[5];
s.t. c1: (if x[1] != x[3] then x[1] else 0) + (if x[1] != x[3] then x[3] else 0) >= 1;
s.t. c2: (x[1] != x[3]) and ((x[1] != x[3]) or not (x[1] != x[3]));
s.t. c3: not ((x[1] != x[3]) and (x[1] != x[3]));
