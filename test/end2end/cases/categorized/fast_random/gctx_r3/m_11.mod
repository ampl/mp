var x{1..3} integer >= 0 <= 2;
maximize o: (6)*x[1] + (-2)*x[2] + (-4)*x[3] + (1)*(if abs(x[1] - x[3]) != 0 then 1 else 0);
s.t. c1: (if x[3] + x[2] == 1 then x[2] else x[1]) == 2;
s.t. c2: (x[3] + x[2] == 1) <==> (x[3] >= x[3]);
s.t. c3: (if x[3] >= x[3] then x[3] else 0) + (if x[3] + x[2] == 1 then x[2] else 0) < 1;
