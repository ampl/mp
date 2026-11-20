var x{1..3} integer >= -1 <= 2;
maximize o: (3)*x[1] + (6)*x[2] + (2)*x[3];
s.t. c1: not ((min(x[2], x[3]) != 1) and (max(x[3], x[3]) != 2));
s.t. c2: not (x[3] != x[3]);
s.t. c3: (if x[3] != x[3] then x[2] else 0) + (if min(x[2], x[3]) != 1 then x[3] else 0) <= 1;
s.t. c4: (min(x[2], x[3]) != 1) <==> (x[3] != x[3]);
