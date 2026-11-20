var x{1..3} integer >= 0 <= 3;
maximize o: (5)*x[1] + (1)*x[2] + (-2)*x[3];
s.t. c1: (x[1] == x[2]) <==> (x[3] >= x[2]);
s.t. c2: not (abs(x[3] - x[3]) <= 3);
s.t. c3: (if x[3] >= x[2] then 1 else 0) + (if x[3] >= x[2] then 1 else 0) + (if x[3] >= x[2] then 1 else 0) >= 1;
s.t. c4: (if x[1] == x[2] then 1 else 0) + (if abs(x[3] - x[3]) <= 3 then 1 else 0) + (if x[3] >= x[2] then 1 else 0) >= 3;
