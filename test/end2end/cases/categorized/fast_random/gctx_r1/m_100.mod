var x{1..4} integer >= 0 <= 2;
minimize o: (-5)*x[1] + (-3)*x[2] + (6)*x[3] + (4)*x[4];
s.t. c1: (if x[4] <= 1 then x[2] else x[4]) <= 0;
s.t. c2: (x[4] < 2) ==> (max(x[2], x[4]) != 0) else (x[4] <= 1);
s.t. c3: (min(x[1], x[2]) >= 0) <==> (max(x[2], x[4]) != 0);
