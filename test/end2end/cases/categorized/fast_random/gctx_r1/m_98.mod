var x{1..3} integer >= -1 <= 1;
minimize o: (6)*x[1] + (5)*x[2] + (-3)*x[3];
s.t. c1: (max(x[2], x[3]) > 0) ==> (min(x[3], x[3]) >= 1);
s.t. c2: (max(x[2], x[3]) > 0) ==> (min(x[3], x[3]) >= 1) else (max(x[2], x[3]) > 0);
