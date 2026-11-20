var x{1..4} integer >= 0 <= 3;
minimize o: (-2)*x[1] + (-4)*x[2] + (3)*x[3] + (-1)*x[4];
s.t. c1: not ((x[1] > 3) and (x[4] + x[4] == 1));
s.t. c2: (x[1] > 3) <==> (max(x[3], x[4]) < 1);
s.t. c3: (x[4] + x[4] == 1) <==> (x[2] == x[3]);
s.t. c4: (if x[2] == x[3] then x[4] else x[3]) <= 0;
