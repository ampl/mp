var x{1..4} integer >= 0 <= 2;
maximize o: (5)*x[1] + (6)*x[2] + (-2)*x[3] + (-2)*x[4];
s.t. c1: numberof 1 in (x[1], x[2], x[3], x[4]) >= 3;
s.t. c2: (x[1] >= 1) <==> (x[1] + x[1] >= 1);
