var x{1..3} integer >= 0 <= 3;
maximize o: (4)*x[1] + (-4)*x[2] + (3)*x[3];
s.t. c1: numberof 0 in (x[1], x[2], x[3]) >= 0;
s.t. c2: abs(x[3] - x[1]) == 2;
s.t. c3: (x[3] <= 2) <==> (x[3] + x[3] > 2);
