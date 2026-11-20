var x{1..3} integer >= 0 <= 3;
minimize o: (-1)*x[1] + (5)*x[2] + (4)*x[3];
s.t. c1: not (max(x[3], x[2]) <= 3);
s.t. c2: (if min(x[1], x[3]) <= 3 then x[2] else 0) + (if min(x[1], x[3]) <= 3 then x[1] else 0) <= 2;
