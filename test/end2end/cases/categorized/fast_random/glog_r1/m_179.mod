var x{1..4} integer >= -2 <= 2;
maximize o: (4)*x[1] + (-5)*x[2] + (0)*x[3] + (3)*x[4];
s.t. c1: (x[4] >= 2) and (x[2] >= 2) ==> (x[2] >= 2);
s.t. c2: (x[2] != 2) ==> (x[1] <= -2) else (x[3] >= -1);
s.t. c3: max(x[1], x[3]) - min(x[1], x[2]) <= -2;
s.t. c4: ((if x[1] != 2 then 1 else 0) + (if x[2] >= 2 then 1 else 0)) >= 1;
