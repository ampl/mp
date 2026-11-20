var x{1..4} integer >= 0 <= 4;
maximize o: (1)*x[1] + (-3)*x[2] + (1)*x[3] + (3)*x[4];
s.t. c1: (x[2] <= 3) <==> (x[3] + x[3] == 1);
s.t. c2: numberof 3 in (x[1], x[2], x[3], x[4]) <= 1;
s.t. c3: ((if x[2] < 2 then 1 else 0) + (if x[3] >= 2 then 1 else 0)) >= 1;
