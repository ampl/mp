var x{1..4} integer >= -1 <= 1;
minimize o: (6)*x[1] + (-3)*x[2] + (-2)*x[3] + (5)*x[4];
s.t. c1: (x[4] == -1) ==> (x[4] == -1);
s.t. c2: (x[4] > 1) <==> (abs(x[2] - x[4]) <= 1);
s.t. c3: (x[4] == -1) or (x[4] > 1);
s.t. c4: (if abs(x[2] - x[4]) <= 1 then 1 else 0) + (if x[2] <= 1 then 1 else 0) + (if x[4] == -1 then 1 else 0) >= 2;
