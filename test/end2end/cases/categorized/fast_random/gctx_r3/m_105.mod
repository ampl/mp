var x{1..4} integer >= -1 <= 2;
maximize o: (6)*x[1] + (4)*x[2] + (-3)*x[3] + (-3)*x[4];
s.t. c1: (x[4] < x[1]) ==> (x[4] == x[1]) else (x[4] < x[1]);
s.t. c2: (if x[4] < x[1] then 1 else 0) + (if x[4] == x[1] then 1 else 0) + (if x[4] == x[1] then 1 else 0) >= 2;
s.t. c3: (x[4] == x[1]) <==> (x[4] < x[1]);
s.t. c4: not (x[4] == x[1]);
