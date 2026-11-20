var x{1..4} integer >= -3 <= 1;
maximize o: (2)*x[1] + (-4)*x[2] + (5)*x[3] + (0)*x[4];
s.t. c1: numberof -3 in (x[1], x[2], x[3], x[4]) >= 1;
s.t. c2: ((if x[2] < 1 then 1 else 0) + (if x[2] >= 1 then 1 else 0)) >= 1;
s.t. c3: (x[4] != -1) ==> (x[4] > -3);
