var x{1..4} integer >= -2 <= 0;
maximize o: (-2)*x[1] + (6)*x[2] + (1)*x[3] + (-4)*x[4];
s.t. c1: ((if x[2] <= -2 then 1 else 0) + (if x[2] >= -2 then 1 else 0)) >= 1;
s.t. c2: (x[1] == 0) and (x[4] < 0) ==> (x[3] >= 0);
s.t. c3: (x[4] >= -2) or (x[3] != -2) or (x[2] < -2);
s.t. c4: (x[4] == -2) ==> (x[2] >= -2) else (x[2] != -1);
