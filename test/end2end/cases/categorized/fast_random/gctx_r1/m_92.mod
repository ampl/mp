var x{1..4} integer >= -1 <= 1;
maximize o: (4)*x[1] + (0)*x[2] + (-4)*x[3] + (1)*x[4];
s.t. c1: (if min(x[3], x[3]) != 1 then x[3] else 0) + (if min(x[3], x[3]) != 1 then x[1] else 0) < 0;
s.t. c2: (min(x[3], x[3]) != 1) and ((x[4] != x[2]) or not (x[4] != x[2]));
s.t. c3: not ((min(x[3], x[3]) != 1) and (min(x[3], x[3]) != 1));
s.t. c4: (if x[4] != x[2] then 1 else 0) + (if min(x[3], x[3]) != 1 then 1 else 0) + (if x[4] != x[2] then 1 else 0) >= 3;
