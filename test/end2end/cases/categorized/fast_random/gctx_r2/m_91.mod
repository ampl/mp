var x{1..3} integer >= 0 <= 3;
minimize o: (3)*x[1] + (6)*x[2] + (5)*x[3] + (5)*(if x[2] > 3 then 1 else 0);
s.t. c1: (if x[2] > 3 then x[2] else 0) + (if x[2] > 3 then x[3] else 0) < 1;
s.t. c2: (x[1] != 2) <==> (x[2] > 3);
