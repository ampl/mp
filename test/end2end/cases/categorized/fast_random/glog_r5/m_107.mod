var x{1..4} integer >= 0 <= 3;
maximize o: (-2)*x[1] + (6)*x[2] + (-2)*x[3] + (-2)*x[4];
s.t. c1: (x[3] != 1) or (x[2] < 1) or (x[3] <= 0);
s.t. c2: ((if x[3] == 0 then 1 else 0) + (if x[4] >= 0 then 1 else 0)) >= 1;
