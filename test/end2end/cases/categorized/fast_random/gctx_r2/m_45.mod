var x{1..3} integer >= 0 <= 2;
minimize o: (-4)*x[1] + (2)*x[2] + (6)*x[3];
s.t. c1: (if x[2] == 1 then x[1] else x[2]) < 0;
s.t. c2: not ((x[2] == 1) and (x[2] == 1));
s.t. c3: (max(x[2], x[2]) <= 2) ==> (x[2] == 1) else (x[2] == 1);
