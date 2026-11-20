var x{1..3} integer >= 0 <= 2;
maximize o: (2)*x[1] + (3)*x[2] + (-2)*x[3];
s.t. c1: (x[2] <= 2) <==> (x[1] + x[3] == 1);
s.t. c2: ((if x[1] != 2 then 1 else 0) + (if x[2] >= 2 then 1 else 0)) >= 1;
s.t. c3: numberof 0 in (x[1], x[2], x[3]) > 0;
s.t. c4: abs(x[1] - x[3]) != 0;
