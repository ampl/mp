var x{1..3} integer >= 0 <= 2;
minimize o: (-1)*x[1] + (3)*x[2] + (6)*x[3];
s.t. c1: not (min(x[2], x[3]) <= 1);
s.t. c2: (min(x[2], x[3]) <= 1) ==> (max(x[3], x[2]) >= 0);
