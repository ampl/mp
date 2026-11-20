var x{1..3} integer >= 0 <= 2;
maximize o: (-3)*x[1] + (3)*x[2] + (-4)*x[3];
s.t. c1: (abs(x[2] - x[3]) >= 2) <==> (x[3] == x[1]);
s.t. c2: (if abs(x[2] - x[3]) >= 2 then 1 else 0) + (if x[3] == x[1] then 1 else 0) + (if x[3] == x[1] then 1 else 0) >= 1;
s.t. c3: (if x[3] == x[1] then 1 else 0) + (if x[3] == x[1] then 1 else 0) + (if x[3] == x[1] then 1 else 0) >= 3;
s.t. c4: (abs(x[2] - x[3]) >= 2) ==> (x[3] == x[1]) else (x[3] == x[1]);
