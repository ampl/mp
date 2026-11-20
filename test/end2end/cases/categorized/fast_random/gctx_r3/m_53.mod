var x{1..5} integer >= 0 <= 2;
minimize o: (3)*x[1] + (5)*x[2] + (5)*x[3] + (6)*x[4] + (5)*x[5];
s.t. c1: (if x[5] + x[2] < 2 then 1 else 0) + (if x[5] + x[2] < 2 then 1 else 0) + (if min(x[1], x[4]) < 2 then 1 else 0) >= 3;
s.t. c2: (x[4] > 0) <==> (x[4] > 0);
s.t. c3: (if min(x[1], x[4]) < 2 then x[3] else x[2]) > 2;
