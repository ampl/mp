var x{1..3} integer >= 0 <= 2;
minimize o: (4)*x[1] + (-3)*x[2] + (-1)*x[3];
s.t. c1: (x[3] <= x[1]) or (x[3] <= x[1]);
s.t. c2: (max(x[2], x[3]) == 0) <==> (max(x[2], x[3]) == 0);
s.t. c3: (if x[3] <= x[1] then x[3] else x[2]) > 1;
