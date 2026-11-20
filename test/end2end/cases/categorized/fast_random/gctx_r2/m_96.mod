var x{1..4} integer >= 0 <= 3;
minimize o: (-4)*x[1] + (5)*x[2] + (6)*x[3] + (-3)*x[4];
s.t. c1: not (x[1] < x[1]);
s.t. c2: (x[3] == x[4]) ==> (x[1] < x[1]) else (x[1] >= 1);
s.t. c3: (if x[1] < x[1] then 1 else 0) + (if x[1] >= 1 then 1 else 0) + (if x[3] == x[4] then 1 else 0) >= 2;
