var x{1..3} integer >= 0 <= 2;
maximize o: (5)*x[1] + (-3)*x[2] + (-1)*x[3];
s.t. c1: (x[3] > 1) <==> (x[1] + x[2] == 2);
s.t. c2: not ((x[3] == 1) and (x[1] >= 1));
s.t. c3: numberof 1 in (x[1], x[2], x[3]) == 1;
