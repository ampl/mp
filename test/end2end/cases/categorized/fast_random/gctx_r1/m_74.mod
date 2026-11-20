var x{1..4} integer >= -1 <= 1;
minimize o: (6)*x[1] + (5)*x[2] + (3)*x[3] + (3)*x[4] + (1)*(if x[2] != x[1] then 1 else 0);
s.t. c1: (abs(x[3] - x[2]) != 1) <==> (abs(x[4] - x[1]) < 1);
s.t. c2: (abs(x[3] - x[2]) != 1) ==> (x[2] != x[1]) else (abs(x[4] - x[1]) < 1);
s.t. c3: (abs(x[3] - x[2]) != 1) <==> (x[2] != x[1]);
s.t. c4: (if x[2] != x[1] then x[4] else 0) + (if x[2] != x[1] then x[1] else 0) >= -1;
