var x{1..3} integer >= 0 <= 3;
maximize o: (6)*x[1] + (2)*x[2] + (-4)*x[3] + (0)*(if x[1] == 3 then 1 else 0);
s.t. c1: (abs(x[3] - x[3]) >= 3) ==> (abs(x[3] - x[3]) >= 3) else (x[1] == 3);
s.t. c2: (x[1] == 3) or (x[1] == 3);
s.t. c3: (if x[1] == 3 then 1 else 0) + (if x[1] == 3 then 1 else 0) + (if abs(x[3] - x[3]) >= 3 then 1 else 0) >= 2;
