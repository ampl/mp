var x{1..4} integer >= 0 <= 3;
minimize o: (3)*x[1] + (5)*x[2] + (4)*x[3] + (-3)*x[4];
s.t. c1: not (x[2] >= x[2]);
s.t. c2: (if x[2] >= x[2] then x[2] else x[2]) <= 0;
s.t. c3: (if min(x[2], x[4]) < 0 then x[4] else x[3]) < 0;
